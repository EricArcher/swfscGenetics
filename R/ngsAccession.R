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
  print("starting function")
  
  .valOrNull <- function(x) {
    if(is.na(x)) "NULL" else {
      if(is.character(x)) paste0("'", x, "'")
      else as.character(x)
    }
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
    print(i)
    labid.num <- as.numeric(
      regmatches(df$labid[i], regexpr("[[:digit:]]*", df$labid[i]))
    )
  
    print("checking if record exists")
    # Check if record exists
    qry.result <- RODBC::sqlQuery(
      conn, 
      paste0(
        "EXEC sp_NextGenSequence_LookupID ",
        .valOrNull(df$run.library[i]), ", ",
        .valOrNull(df$original.filename[i])
      )
    )
    print(qry.result)
    if(is.character(qry.result)) stop(qry.result)
    id <- as.numeric(unlist(qry.result))
    fname <- NA
    
    print("creating row")
    if(id != 0) {
      x <- if(id < 0) "more than one" else "one"
      message(
        format(Sys.time()), " : Found ", x, " record matching: ",
        "Run Library = '", df$run.library[i], "' and ",
        "Original Filename = '", df$original.filename[i], "'"
      )
    } else {
      # Insert row and get new ID
      qry.result <- RODBC::sqlQuery(
        conn,
        paste0(
          "EXEC sp_NextGenSequence_Insert ",
          labid.num, ", ", 
          .valOrNull(df$run.library[i]), ", ",
          .valOrNull(df$i7.index[i]), ", ",
          .valOrNull(df$original.filename[i]), ", ",
          .valOrNull(df$d.id[i]), ", ",
          .valOrNull(df$read.direction[i]), ", ",
          .valOrNull(df$i5.index[i])
        )
      )
      message(qry.result)
      if(is.character(qry.result)) stop(qry.result)
      id <- as.numeric(unlist(qry.result))
      message(
        format(Sys.time()), 
        " : Inserted ID: ", id, ", LABID: ", df$labid[i]
      )
      
      # Update filename
      fname <- paste0(
        "z", sprintf("%07d", labid.num), 
        "_", df$species[i], 
        "_", df$run.library[i], 
        "_n", sprintf("%07d", id), 
        "_", df$read.direction[i],
        ".fastq.gz"
      ) 
      qryStr <- paste0(
        "SET NOCOUNT ON UPDATE tbl_NextGenSequence ", 
        "SET New_Filename = '", fname, "' WHERE ID = ", id
      )
      message(qryStr)
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