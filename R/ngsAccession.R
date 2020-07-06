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
  ts <- format(Sys.time(), "%Y%m%d_%H%M")
  
  connStr <- paste(
    "DRIVER=SQL Server Native Client 11.0",
    "DATABASE=Genetics",
    "Trusted_Connection=Yes",
    "SERVER=161.55.235.187",
    sep = ";"
  )
  RODBC::odbcCloseAll()
  conn <- RODBC::odbcDriverConnect(connection = connStr)
  
  df$labid.num <- as.numeric(
    regmatches(df$labid, regexpr("[[:digit:]]*", df$labid))
  )
  
  result <- do.call(rbind, lapply(1:nrow(df), function(i) {
    # Insert row
    qryStr <- paste0(
      "EXEC sp_NextGenSequence_Insert ",
      df$labid.num[i], ", ", 
      ifelse(
        is.na(df$run.library[i]), 
        "NULL", 
        paste("'", df$run.library[i], "'", sep = "")
      ), ", ", 
      ifelse(is.na(df$i7.index[i]), "NULL", df$i7.index[i]), ", ", 
      "NULL, ", 
      ifelse(
        is.na(df$original.filename[i]), 
        "NULL", 
        paste("'", df$original.filename[i], "'", sep = "")
      ), ", ",
      ifelse(
        is.na(df$d.id[i]), 
        "NULL", 
        paste("'", df$d.id[i], "'", sep = "")
      ), ", ", 
      ifelse(
        is.na(df$read.direction[i]), 
        "NULL", 
        paste("'", df$read.direction[i], "'", sep = "")
      ), ", ", 
      ifelse(is.na(df$i5.index[i]), "NULL", df$i5.index[i])
    )
    
    # Get ID
    id <- unlist(RODBC::sqlQuery(conn, qryStr))
    message("inserting id:, ", id, ", LABID: ", df$labid[i])
    fname <- NA
    
    if(id > 0) {
      # Create new filename
      fname <- paste(
        "z", sprintf("%07d", df$labid.num[i]), 
        "_", df$species[i], 
        "_", df$run.library[i], 
        "_n", sprintf("%07d", id), 
        "_", df$read.direction[i],
        ".fastq.gz", 
        sep = ""
      )      
      
      # Update filename
      qryStr <- paste(
        "SET NOCOUNT ON ",
        "UPDATE tbl_NextGenSequence ",
        "SET New_Filename = '", fname, "' ",
        "WHERE ID = ", id,
        sep = ""
      )
      RODBC::sqlQuery(conn, qryStr)
    }
    
    data.frame(
      index.id = df$index.id[i],
      ngs.id = id, 
      new.filename = fname, 
      stringsAsFactors = FALSE
    )
  }))
  RODBC::odbcCloseAll()
  
  acc.df <- merge(df, result, by = "index.id", all.x = TRUE)
  library.name <- paste(unique(sort(acc.df$run.library)), collapse = ".")
  acc.fname <- paste0(library.name, "_library_accession_report_", ts, ".csv")
  utils::write.csv(acc.df, file = acc.fname, row.names = FALSE)
  
  message(
    format(Sys.time()), " : Accessioned ", 
    sum(result$ngs.id > 0), " of ", nrow(result), " records."
  )
  
  acc.df
}