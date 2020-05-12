#' @title Read pileup file
#' @description Reads a pileup format file created by samtools mpileup
#' 
#' @param fname filename of pileup file
#' @param all.columns If \code{FALSE} (default), only the first 6 columns are 
#'   returned. If \code{TRUE}, all columns are returned, but only the six are 
#'   informatively named.
#' 
#' @return data frame representing the pileup file with the following columns:
#' \tabular{ll}{ 
#'   \code{chrom} \tab Chromosome name.\cr 
#'   \code{pos} \tab 1-based position on the chromosome.\cr 
#'   \code{ref} \tab Reference base at this position.\cr 
#'   \code{cov} \tab Number of reads covering this position.\cr 
#'   \code{bases} \tab Read bases.\cr 
#'   \code{quals} \tab Base qualities, encoded as ASCII characters.\cr 
#' }
#' 
#' @note The input pileup file should be the result of a call to 
#'   'samtools mpileup' on a single BAM file. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
readPileup <- function(fname, all.columns = FALSE) {
  col.names <- c("chrom", "pos", "ref", "cov", "bases", "quals")
  plp <- data.table::fread(fname, sep = "\t")
  colnames(plp)[1:6] <- col.names
  plp$ref <- toupper(plp$ref)
  plp$bases <- toupper(plp$bases)
  plp <- as.data.frame(plp)
  if(all.columns) plp else plp[, col.names]
}