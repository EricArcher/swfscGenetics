#' @title Create consensus sequences from pileup file
#' @description Reads a pileup formatted file created by \code{samtools mpileup}
#'   and creates a consensus sequence for each chromosome listed.
#'
#' @param fname filename of pileup file
#' @param min.cov minimum coverage for base calling. Sites with coverage below
#'   this are assigned \code{N}'s.
#' @param min.freq minimum frequency of either the reference or alternate base
#'   for calling. If both bases are below this frequency, an \code{N} 
#'   is assigned.
#' @param min.freq.cov minimum coverage above which \code{min.freq} is applied. 
#'   Sites below this value and >= than \code{min.cov} will only be called 
#'   if all reads agree.
#' @param num.cores number of cores to use during processing. If \code{NULL}, 
#'   will default to \code{parallel::detectCores() - 1}.
#'
#' @return list with the following elements:
#' \tabular{ll}{ 
#'   \code{cons.seq} \tab named list of consensus sequences for each chromosome.\cr 
#'   \code{plp} \tab data frame of the pileup data with base frequencies.\cr 
#'   \code{insertions} \tab data frame of insertion locations and frequencies.
#'     If there are no insertions, this will be \code{NULL}.\cr 
#' }
#'
#' @note The input pileup file should be the result of a call to 
#'   \code{samtools mpileup} on a single BAM file.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#' 
pileupConsensus <- function(fname, min.cov, min.freq, min.freq.cov, 
                            num.cores = NULL) {
  # convert parameters and make sure the values are in proper range
  params <- list(
    min.cov = max(as.numeric(min.cov), 1),
    min.freq = min(max(0.5, as.numeric(min.freq)), 1),
    min.freq.cov = max(min.cov, as.numeric(min.freq.cov))
  )
  
  message("Reading file...")
  plp.file <- readPileup(fname)
  plp.file <- split(plp.file, plp.file$chrom)
  
  smry <- lapply(1:length(plp.file), function(i) {
    message(
      "Processing pileup for '", names(plp.file)[[i]], "' ",
      "(", i, "/", length(plp.file), ")..."
    ) 
    
    # call bases and summarize frequencies
    plp.smry <- .plpSummary(plp.file[[i]], params, num.cores)
    
    # compile insertion frequencies
    insertions <- .compileInsertions(plp.smry)
    
    # compile pileup summary data frames
    plp.smry <- do.call(rbind, lapply(plp.smry, function(x) x$plp))
    plp.smry <- plp.smry[order(plp.smry$pos, plp.smry$pos.order), ]
    plp.smry$cons.pos <- 1:nrow(plp.smry)
    rownames(plp.smry) <- NULL
    colnames(plp.smry)[2] <- "ref.pos"
    
    cons.seq <- unlist(strsplit(plp.smry$consensus, ""))
    
    list(cons.seq = cons.seq, plp = plp.smry, insertions = insertions)
  })
  
  .compileSummaries(smry, names(plp.file))
}