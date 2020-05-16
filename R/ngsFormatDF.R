#' @title Format NGS data frame
#' @description Format NGS data frame for accessioning. Assigns correct column 
#'   names, checks that files found.
#' 
#' @param library.name name of run library
#' @param df original library data frame
#' @param folder folder where original FASTQ files are stored
#' @param species species name
#' @param f.index forward index
#' @param r.index reverse index
#' @param labid SWFSC LABID
#' @param filename filename of original FASTQ file
#' @param read.direction direction of reads
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'
ngsFormatDF <- function(library.name, df, folder,
                        species = "species", f.index = "F_index", 
                        r.index = "R_index", labid = "LabID", 
                        filename = "filename", 
                        read.direction = "Read_Direction") {
  df$index.id <- 1:nrow(df)
  df$run.library <- library.name
  df$species <- df[[species]]
  df$species <- gsub(" ", "", df$species)
  df$f.index <- df[[f.index]]
  df$r.index <- df[[r.index]]
  df$labid <- df[[labid]]
  df$filename <- df[[filename]]
  df$read.direction <- df[[read.direction]]
  
  old.fnames <- dir(folder)
  filenames <- do.call(rbind, lapply(1:nrow(df), function(i) {
    #added to index file from directory list for Illumina NextSeq files.
    idx <- df$filename[i]
    fnames <- grep(idx, old.fnames, value = TRUE)
    #fixes problem of grabbing numbers that contain shorter numbers (e.g. 371, 11371)
    #  fnames <- grep(paste0('^(',idx,')$'), old.fnames, value = TRUE) 
    fnames <- if(length(fnames) == 0) NA else fnames
    data.frame(
      index.id = rep(i, length(fnames)), 
      original.filename = fnames,
      stringsAsFactors = FALSE
    )
  }))
  
  df <- merge(df, filenames, by = "index.id", all = TRUE)
  df[!is.na(df$original.filename), ]
}
