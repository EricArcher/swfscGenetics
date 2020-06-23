#' @title Rename NGS files
#' @description Rename NGS files using new accessioned names from indexing
#'   data frame.
#' 
#' @param df formatted indexing data frame resulting 
#'   from \code{\link{ngsFormatDF}}.
#' @param old.folder folder where original files reside.
#' @param new.folder folder where renamed files should be placed.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{ngsFormatDF}}, \code{\link{ngsAccession}}
#' 
#' @export
#'
ngsRename <- function(df, old.folder, new.folder) {
  df <- df[order(df$species, decreasing = T), ]
  
  df$file.written <- sapply(1:nrow(df), function(i) {
    species <- df$species[i]
    run.library <- df$run.library[i]
    new.folder <- file.path(new.folder, species, run.library)
    if(!dir.exists(new.folder)) dir.create(new.folder, recursive = TRUE)
    new.path <- file.path(new.folder, df$new.filename[i])
    if(!file.exists(new.path)) {
      old.path <- file.path(old.folder, df$original.filename[i])
      message(
        format(Sys.time()), 
        " : Renaming ", i, " / ", nrow(df), 
        " '", df$new.filename[i], "'"
      )
      file.rename(from = old.path, to = new.path)
    } else {
      cat("'", new.path, "' already exists\n")
      FALSE
    }
  })
  
  message(format(Sys.time()), " : Renamed ", sum(df$file.written), " files")
  
  library.name <- paste(unique(sort(df$run.library)), collapse = ".")
  out.fname <- paste0(library.name, "_run.library.csv")
  utils::write.csv(df, out.fname, row.names = FALSE)
}