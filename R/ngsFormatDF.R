#' @title Format NGS data frame
#' @description Format NGS data frame for accessioning. Assigns correct column 
#'   names, checks that files found.
#' 
#' @param library.name name of run library. This should also be the folder 
#'   name where the files to be archived are located.
#' @param library.filename filename of library spreadsheet. If \code{NULL}, 
#'   the .csv file in the \code{library.name} folder
#' @param labid SWFSC LABID column name
#' @param d.id DNA ID column name
#' @param species species name column name
#' @param i7.index forward index column name
#' @param i5.index reverse index column name
#' @param read.direction direction of reads column name
#' @param received.filename original filename of \code{.fastq.gz} file. This 
#'   column is optional. If it does not exist, or values are empty for a row,
#'   the function will look for a file with '\code{<labid>_}' and 
#'   '\code{<read.direction>_}' in the filename.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{ngsRename}}, \code{\link{ngsAccession}}
#' 
#' @export
#'
ngsFormatDF <- function(library.name, library.filename = NULL,
                        labid = "LabID", d.id = "D_id", species = "species", 
                        i7.index = "i7_index", i5.index = "i5_index", 
                        read.direction = "Read_Direction",
                        received.filename = "received_filename") {
  
  if(!dir.exists(library.name)) {
    stop("the folder '", library.name, "' does not exist.")
  }
  
  if(is.null(library.filename)) {
    library.filename <- dir(library.name, pattern = ".csv$", full.names = TRUE)
    if(length(library.filename) > 1) {
      stop("more than one .csv file found in '", library.name, "'")
    }
    if(length(library.filename) == 0) {
      stop("no .csv file found in '", library.name, "'")
    }
  }
  
  df <- if(file.exists(library.filename)) {
    utils::read.csv(
      library.filename, 
      na.strings = c("", "NA"),
      stringsAsFactors = FALSE
    )
  } else {
    stop("the file, '", library.filename, "' cannot be found.")
  }
  
  for(x in c(labid, d.id, species, i7.index, i5.index, read.direction)) {
    if(!x %in% colnames(df)) {
      stop("can't find column '", x, "' in '", library.filename, "'")
    }
  }
  
  df$index.id <- 1:nrow(df)
  df$run.library <- gsub("[[:space:]|[:punct:]]+", ".", library.name)
  df$labid <- df[[labid]]
  df$d.id <- df[[d.id]]
  df$species <- df[[species]]
  df$species <- gsub(" ", "", df$species)
  df$i7.index <- df[[i7.index]]
  df$i5.index <- df[[i5.index]]
  df$read.direction <- df[[read.direction]]
  df$received.filename <- if(received.filename %in% colnames(df)) {
    df[[received.filename]]
  } else {
    as.character(rep(NA, nrow(df)))
  }
  
  old.fnames <- dir(
    library.name, 
    pattern = ".fastq.gz$", 
    full.names = TRUE, 
    recursive = TRUE
  )
  df$original.filename <- sapply(1:nrow(df), function(i) {
    if(!is.na(df$received.filename[i])) df$received.filename[i] else {
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
    }
  })
  
  library.name <- paste(unique(sort(df$run.library)), collapse = ".")
  fname <- paste0(
    library.name, 
    "_index_", 
    format(Sys.time(), "%Y%m%d_%H%M"), 
    ".csv"
  )
  utils::write.csv(df, file = fname, row.names = FALSE)
  
  df
}
