#' @title Format NGS data frame
#' @description Format NGS data frame for accessioning. Assigns correct column 
#'   names, checks that files found.
#' 
#' @param library.name name of run library
#' @param library.filename filename of library spreadsheet
#' @param folder folder where original FASTQ files are stored
#' @param labid SWFSC LABID column name
#' @param species species name column name
#' @param d.id DNA ID column name
#' @param i7.index forward index column name
#' @param i5.index reverse index column name
#' @param read.direction direction of reads column name
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{ngsRename}}, \code{\link{ngsAccession}}
#' 
#' @export
#'
ngsFormatDF <- function(library.name, library.filename, folder,
                        labid = "labid", species = "species", d.id = "d_id",
                        i7.index = "i7_index", i5.index = "i5_index", 
                        read.direction = "read_direction") {
  df <- if(file.exists(library.filename)) {
    utils::read.csv(library.filename, stringsAsFactors = FALSE)
  } else {
    stop("the file, '", library.filename, "' cannot be found.")
  }
  
  ts <- format(Sys.time(), "%Y%m%d_%H%M")
  
  df$index.id <- 1:nrow(df)
  df$run.library <- library.name
  df$species <- df[[species]]
  df$species <- gsub(" ", "", df$species)
  df$i7.index <- df[[i7.index]]
  df$i5.index <- df[[i5.index]]
  df$labid <- df[[labid]]
  df$read.direction <- df[[read.direction]]
  
  old.fnames <- dir(
    folder, 
    pattern = ".fastq", 
    full.names = TRUE, 
    recursive = TRUE
  )
  df$original.filename <- sapply(1:nrow(df), function(i) {
    has.labid <- grep(
      paste0(df$labid[i], "_"), 
      old.fnames, 
      value = TRUE
    )
    has.read.dir <- grep(
      paste0(df$read.direction[i], "_"), 
      old.fnames, 
      value = TRUE
    )
    fname <- intersect(has.labid, has.read.dir)
    if(length(fname) != 1) NA else basename(fname)
  })
  
  library.name <- paste(unique(sort(df$run.library)), collapse = ".")
  fname <- paste0(library.name, "_index_", ts, ".csv")
  utils::write.csv(df, file = fname, row.names = FALSE)
  
  df
}
