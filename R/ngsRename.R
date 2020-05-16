#' @title Rename NGS files
#' @description Rename NGS files using new accessioned names from indexing
#'   data frame.
#' 
#' @param df formatted indexing data frame resulting 
#'   from \code{\link{ngsFormatDF}}.
#' @param folder folder where original files reside.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'
ngsRename <- function(df, folder) {
  df <- df[order(df$Gspp, decreasing = T), ]
  df$file.written <- sapply(1:nrow(df), function(i) {
    new.path <- file.path(folder, df$New.Filename[i])
    if(!file.exists(new.path)) {
      old.path <- file.path(folder, df$Original.Filename[i])
      message(
        format(Sys.time()), 
        " : Renaming ", i, " / ", nrow(df), 
        " '", df$New.Filename[i], "'"
      )
      file.rename(from = old.path, to = new.path)
    } else FALSE
  })
  message(format(Sys.time()), " : Renamed ", sum(df$file.written), " files")
}