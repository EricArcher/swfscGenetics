#' @title Rename NGS files
#' @description Rename NGS files using new accessioned names from indexing
#'   data frame.
#' 
#' @param df formatted indexing data frame resulting 
#'   from \code{\link{ngsFormatDF}}.
#' @param old.folder folder where original files reside.
#' @param new.folder folder where renamed files should be placed.
#' @param leave.files leave files in \code{old.folder}? If \code{FALSE} files 
#'   will be moved.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{ngsFormatDF}}, \code{\link{ngsAccession}}
#' 
#' @export
#'
ngsRename <- function(df, old.folder, new.folder, leave.files = TRUE) {
  df <- df[order(df$species, decreasing = T), ]
  
  df$file.written <- sapply(1:nrow(df), function(i) {
    if(is.na(df$original.filename[i]) | is.na(df$new.filename[i])) return(FALSE)
    # get old filename and check that it exists
    old.path <- dir(
      old.folder,
      pattern = paste0("^", df$original.filename[i], "$"),
      full.names = TRUE,
      recursive = TRUE
    )
    if(length(old.path) != 1) {
      message("'", df$orginal.filename[i], "' can't be found in, '", old.folder, "'")
      return(FALSE)
    }
    # create new folder if necessary
    new.folder <- file.path(new.folder, df$species[i], df$run.library[i])
    if(!dir.exists(new.folder)) dir.create(new.folder, recursive = TRUE)
    # create new path and check that it exists
    new.path <- file.path(new.folder, df$new.filename[i])
    if(file.exists(new.path)) {
      message("'", new.path, "' already exists.")
      return(FALSE)
    }
    # rename files
    message(
      format(Sys.time()), 
      " : Renaming ", i, " / ", nrow(df), 
      " '", df$new.filename[i], "'"
    )
    if(leave.files) {
      file.copy(from = old.path, to = new.path)
    } else {
      file.rename(from = old.path, to = new.path)
    }
  })
  
  message(format(Sys.time()), " : Renamed ", sum(df$file.written), " files")
  
  library.name <- paste(unique(sort(df$run.library)), collapse = ".")
  out.fname <- paste0(library.name, "_run.library.csv")
  utils::write.csv(df, out.fname, row.names = FALSE)
}