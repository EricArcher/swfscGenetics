#' @title Accession NGS sequence metadata to SQL Server table
#' @description Accession NGS sequence metadata to SQL Server table.
#' 
#' @param df formatted indexing data frame resulting 
#'   from \code{\link{ngsFormatDF}}.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{ngsFormatDF}}, \code{\link{ngsRename}}
#' 
#' @export
#' 
ngsAccession <- function(df) {
  .valOrNull <- function(x) {
    if(is.na(x)) "NULL" else {
      if(is.character(x)) paste0("'", x, "'")
      else as.character(x)
    }
  }
  
  .makeNewFilename <- function(labid.num, df, i, id) {
    paste0(
      "z", sprintf("%07d", labid.num), 
      "_", df$species[i], 
      "_", df$run.library[i], 
      "_n", sprintf("%07d", id), 
      "_", df$read.direction[i],
      ".fastq.gz"
    ) 
  }
  
  .lookupQuery <- function(df, i) {
    paste0(
      # SET NOCOUNT ON added to prevent extra result sets from 
      # interfering with SELECT statements.
      "SET NOCOUNT ON; ", 
      
      #0 is the default "not found" value
      "DECLARE @Return int = 0, @Count int, @MaxID int; ", 
      
      # This returns the count of matching records, and the Max of their FASTQ_ID's. 
      # If exactly one record was found, that Max will be the one FASTQ_ID. "
      "SELECT @Count = COUNT(*), @MaxID = MAX(FASTQ_ID) ",
      "FROM tbl_NGS_FASTQ ",
      "WHERE ISNULL(Run_Library, '') = ",
      "ISNULL(", .valOrNull(df$run.library[i]), ", '') ",
      "AND ISNULL(Original_FIlename, '') = ",
      "ISNULL(", .valOrNull(df$original.filename[i]), ", ''); ",
      
      # If there was exactly 1 matching record, return its FASTQ_ID
      "IF @Count = 1 SET @Return = @MaxID ELSE IF @Count > 1 SET @Return = -1; ",
      
      "SELECT @Return;"
    )
  }
  
  .insertQuery <- function(df, i, labid.num) {
    paste0(
      # SET NOCOUNT ON added to prevent extra result sets from 
      # interfering with SELECT statements.
      "SET NOCOUNT ON; ", 
      
      "IF NOT EXISTS ( ",
      # check if record has already been entered
      "SELECT * FROM tbl_NGS_FASTQ ",
      "WHERE LABID = ", labid.num,  
      " AND Run_Library = ", .valOrNull(df$run.library[i]),
      " AND i7_index = ", .valOrNull(df$i7.index[i]),
      " AND Library_Directory = ", .valOrNull(df$library.directory[i]),
      " AND Original_Filename = ", .valOrNull(df$original.filename[i]),
      " AND D_id = ", .valOrNull(df$d.id[i]),
      " AND Read_direction = ", .valOrNull(df$read.direction[i]),
      " AND i5_index = ", .valOrNull(df$i5.index[i]),
      ") BEGIN ",
      # enter record data
      "INSERT INTO tbl_NGS_FASTQ (",
      "LABID, Run_Library, i7_index, Library_Directory, Original_Filename, ",
      "D_id, Read_direction, i5_index, Comments) ",
      "VALUES (",
      labid.num, ", ",
      .valOrNull(df$run.library[i]), ", ",
      .valOrNull(df$i7.index[i]), ", ",
      .valOrNull(df$library.directory[i]), ", ",
      .valOrNull(df$original.filename[i]), ", ",
      .valOrNull(df$d.id[i]), ", ",
      .valOrNull(df$read.direction[i]), ", ",
      .valOrNull(df$i5.index[i]), ", ",
      .valOrNull(df$comments[i]),
      "); ",
      "SELECT SCOPE_IDENTITY(); ",
      "END ",
      # return 0 if record has already been entered
      "ELSE SELECT 0;"
    )
  }
  
  RODBC::odbcCloseAll()
  conn <- RODBC::odbcDriverConnect(
    paste(
      "DRIVER=SQL Server Native Client 11.0",
      "DATABASE=Genetics",
      "Trusted_Connection=Yes",
      "SERVER=161.55.235.187",
      sep = ";"
    )
  )
  
  result <- do.call(rbind, lapply(1:nrow(df), function(i) {
    labid.num <- as.numeric(
      regmatches(df$labid[i], regexpr("[[:digit:]]*", df$labid[i]))
    )
  
    # Check if record exists
    qry.result <- RODBC::sqlQuery(conn, .lookupQuery(df, i))
    if(is.character(qry.result)) stop(qry.result)
    id <- as.numeric(unlist(qry.result))
    
    fname <- if(id < 0) {
      message(
        format(Sys.time()), " : WARNING! - Found more than one record matching: ",
        "Run Library = '", df$run.library[i], "' and ",
        "Original Filename = '", df$original.filename[i], "'"
      )
      id <- NA
      NA
    } else if(id > 0) {     
      message(
        format(Sys.time()), " : Found ngs.id = ", id, " matching: ",
        "Run Library = '", df$run.library[i], "' and ",
        "Original Filename = '", df$original.filename[i], "'"
      )
      .makeNewFilename(labid.num, df, i, id)
    } else {
      # Insert row and get new ID
      qry.result <- RODBC::sqlQuery(conn, .insertQuery(df, i, labid.num))
      if(is.character(qry.result)) stop(qry.result)
      id <- as.numeric(unlist(qry.result))
      message(
        format(Sys.time()), 
        " : Inserted ID: ", id, ", LABID: ", df$labid[i]
      )
      
      # Create file name and update database
      fname <- .makeNewFilename(labid.num, df, i, id)
      qryStr <- paste0(
        "SET NOCOUNT ON ",
        "UPDATE tbl_NGS_FASTQ ", 
        "SET New_Filename = '", fname, "' WHERE FASTQ_ID = ", id
      )
      RODBC::sqlQuery(conn, qryStr)
      fname
    }
    
    data.frame(
      index.id = df$index.id[i],
      ngs.id = id, 
      new.filename = fname, 
      stringsAsFactors = FALSE
    )
  }))
  RODBC::odbcCloseAll()
  
  message(
    format(Sys.time()), " : Accessioned ", 
    sum(result$ngs.id > 0), " of ", nrow(result), " records."
  )
  acc.df <- merge(df, result, by = "index.id", all.x = TRUE)
  acc.fname <- paste0(
    paste(unique(sort(acc.df$run.library)), collapse = "."), 
    "_library_accession_report_", 
    format(Sys.time(), "%Y%m%d_%H%M"),
    ".csv"
  )
  utils::write.csv(acc.df, file = acc.fname, row.names = FALSE)
  
  acc.df
}