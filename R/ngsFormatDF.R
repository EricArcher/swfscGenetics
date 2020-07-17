#' @title Format NGS data frame
#' @description Format NGS data frame for accessioning. Assigns correct column 
#'   names, checks that files found.
#' 
#' @param library.name name of run library. This should also be the folder 
#'   name where the files to be archived are located.
#' @param library.filename filename of library spreadsheet. If \code{NULL}, 
#'   the .csv file in the \code{library.name} folder
#' @param labid SWFSC LABID column name
#' @param d.id DNA ID column name (optional)
#' @param species species name column name
#' @param i7.index forward index column name
#' @param i5.index reverse index column name
#' @param read.direction direction of reads column name
#' @param library.directory original folder of file within the run folder 
#'   specified by \code{library.name} (optional)
#' @param original.filename original filename of \code{.fastq.gz} file 
#'   (optional). If the column does not exist, or values are empty for a row,
#'   the function will look for a file beginning with '\code{<labid>_}' and has 
#'   '\code{<read.direction>_}' somewhere in the filename
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{ngsRename}}, \code{\link{ngsAccession}}
#' 
#' @export
#'
ngsFormatDF <- function(library.name, library.filename = NULL,                     
                        labid = "LABID", d.id = "D_id", species = "species", 
                        i7.index = "i7_index",  i5.index = "i5_index", 
                        read.direction = "Read_direction",
                        library.directory = "Library_Directory",
                        original.filename = "Original_Filename") {
  
  # Check that run folder exists
  if(!dir.exists(library.name)) {
    stop("The folder '", library.name, "' does not exist.")
  }
  
  # Get filename
  if(is.null(library.filename)) {
    library.filename <- dir(library.name, pattern = ".csv$", full.names = TRUE)
    if(length(library.filename) > 1) {
      stop("More than one .csv file found in '", library.name, "'.")
    }
    if(length(library.filename) == 0) {
      stop("No .csv file found in '", library.name, "'.")
    }
  }
  
  # Read file
  df <- if(file.exists(library.filename)) {
    utils::read.csv(
      library.filename, 
      na.strings = c("", "NA"),
      stringsAsFactors = FALSE
    )
  } else {
    stop("The file, '", library.filename, "' cannot be found.")
  }
  
  # Check that column names are in data frame and that columns have data
  for(x in c(labid, species, i7.index, i5.index, read.direction)) {
    if(!x %in% colnames(df)) {
      stop("Can't find column '", x, "' in '", library.filename, "'.")
    }
    if(any(is.na(df[[x]]))) {
      stop("Required column '", x, "' is missing data.")
    }
  }
  if(!d.id %in% colnames(df)) df[[d.id]] <- NA
  if(!library.directory %in% colnames(df)) df[[library.directory]] <- NA
  
  # Check that species is 4 characters long
  df[[species]] <- gsub("[[:space:]]|[[:punct:]]+", "", df[[species]])
  num.chars <- nchar(df[[species]])
  if(any(num.chars != 4)) stop("Not all species are 4 characters long.")
  
  # Add and format necessary columns
  df$index.id <- 1:nrow(df)
  df$run.library <- gsub("[[:space:]|[:punct:]]+", ".", library.name)
  df$labid <- df[[labid]]
  df$d.id <- df[[d.id]]
  df$species <- df[[species]]
  df$i7.index <- df[[i7.index]]
  df$i5.index <- df[[i5.index]]
  df$read.direction <- df[[read.direction]]
  df$library.directory <- df[[library.directory]]
  df$original.filename <- if(original.filename %in% colnames(df)) {
    df[[original.filename]]
  } else {
    as.character(rep(NA, nrow(df)))
  }
  
  # Match original filenames
  old.fnames <- dir(
    library.name, 
    pattern = ".fastq.gz$", 
    full.names = TRUE, 
    recursive = TRUE
  )
  df$original.filename <- sapply(1:nrow(df), function(i) {
    if(!is.na(df$original.filename[i])) df$original.filename[i] else {
      has.labid <- grep(
        paste0("^", df$labid[i], "_"), 
        basename(old.fnames), 
        value = TRUE
      )
      has.read.dir <- grep(
        paste0(df$read.direction[i], "_"), 
        basename(old.fnames), 
        value = TRUE
      )
      fname <- intersect(has.labid, has.read.dir)
      if(length(fname) != 1) NA else basename(fname)
    }
  })
  
  # Check that no filenames are duplicated in df
  is.dup <- duplicated(df$original.filename)
  if(any(is.dup & !is.na(df$original.filename))) {
    stop(
      "The following filenames were assigned to more than one record: ",
      paste(df$original.filename[is.dup], collapse = ", ")
    )
  }
  
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
