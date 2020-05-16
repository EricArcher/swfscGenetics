#' @title Accession NGS sequence metadata to SQL Server table
#' @description Accession NGS sequence metadata to SQL Server table.
#' 
#' @param df formatted indexing data frame resulting 
#'   from \code{\link{ngsFormatDF}}.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
ngsAccession <- function(df) {
  connStr <- paste(
    "DRIVER=SQL Server Native Client 11.0",
    "DATABASE=Genetics",
    "Trusted_Connection=Yes",
    "SERVER=161.55.235.187",
    sep = ";"
  )
  RODBC::odbcCloseAll()
  conn <- RODBC::odbcDriverConnect(connection = connStr)
  
  result <- do.call(rbind, lapply(1:nrow(df), function(i) {
    # Insert row
    qryStr <- paste0(
      "EXEC sp_NextGenSequence_Insert ",
      df$LABID[i], ", ", 
      ifelse(
        is.na(df$Run.Library[i]), 
        "NULL", 
        paste("'", df$Run.Library[i], "'", sep = "")
      ), ", ", 
      ifelse(is.na(df$F_index[i]), "NULL", df$F_index[i]), ", ", 
      "NULL, ", 
      ifelse(
        is.na(df$Original.Filename[i]), 
        "NULL", 
        paste("'", df$Original.Filename[i], "'", sep = "")
      ), ", ",
      ifelse(
        is.na(df$D_id[i]), 
        "NULL", 
        paste("'", df$D_id[i], "'", sep = "")
      ), ", ", 
      ifelse(
        is.na(df$Read_Direction[i]), 
        "NULL", 
        paste("'", df$Read_Direction[i], "'", sep = "")
      ), ", ", 
      ifelse(is.na(df$R_index[i]), "NULL", df$R_index[i])
    )
    
    # Get ID
    id <- unlist(RODBC::sqlQuery(conn, qryStr))
    message("inserting id:, ", id, ", LABID: ", df$LABID[i])
    fname <- NA
    
    if(id > 0) {
      # Create new filename
      fname <- paste(
        "z", sprintf("%07d", df$LABID[i]), 
        "_", df$Gspp[i], 
        "_", df$Run.Library[i], 
        "_n", sprintf("%07d", id), 
        "_", df$Read_Direction[i], #added March2017
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
      New.Filename = fname, 
      stringsAsFactors = FALSE
    )
  }))
  RODBC::odbcCloseAll()

  message(
    format(Sys.time()), " : Accessioned ", 
    sum(result$ngs.id > 0), " of ", nrow(result), " records."
  )
  result
}