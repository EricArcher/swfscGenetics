#' @title Create consensus sequences from pileup file
#' @description Reads a pileup formatted file created by \code{samtools mpileup}
#'   and creates a consensus sequence for each chromosome listed.
#'
#' @param fname filename of pileup file
#' @param min.cov minimum coverage for base calling (sites with coverage below
#'   this are assigned N's).
#' @param min.freq minimum frequency of either the reference or alternate base
#'   for calling. If both bases are below this frequency, an N is assigned.
#' @param min.freq.cov minimum coverage above which min.freq is applied. Sites
#'   below this and >= than min.cov will only be called if all reads agree.
#'
#' @return  list with the following elements:
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
pileupConsensus <- function(fname, min.cov, min.freq, min.freq.cov) {
  # convert parameters and make sure the values are in proper range
  params <- list(
    min.cov = max(as.numeric(min.cov), 1),
    min.freq = min(max(0.5, as.numeric(min.freq)), 1),
    min.freq.cov = max(min.cov, as.numeric(min.freq.cov))
  )
  
  plp.file <- readPileup(fname)
  
  results <- sapply(split(plp.file, plp.file$chrom), function(plp) {
    pos.freqs <- .calcFreqs(plp)
    selected <- .makeSelections(pos.freqs, params)
    cons.seq <- selected$bases
    
    ins <- selected$ins
    insertions <- if(is.null(ins)) NULL else {
      ins.pos <- names(ins)
      # insert insertions at correct position
      i <- match(ins.pos, names(cons.seq))
      cons.seq <- R.utils::insert(cons.seq, i + 1, ins)
      # compile data frame about insertion frequencies
      ins <- sapply(ins, paste, collapse = "")
      ins.freq <- sapply(ins.pos, function(x) max(pos.freqs[[x]]$in.freq))
      ins.df <- data.frame(
        ref.pos = as.numeric(ins.pos),
        insertion = do.call(c, ins),
        freq = ins.freq
      )
      ins.df <- cbind(
        chrom = rep(unique(plp$chrom), nrow(ins.df)),
        ins.df
      )
      rownames(ins.df) <- NULL
      ins.df
    }
    
    plp <- .plpSummary(plp, cons.seq, pos.freqs)
    
    # replace names of inserted bases with position of prior base
    insert.i <- which(is.na(as.numeric(names(cons.seq))))
    if(length(insert.i) > 0) {
      for(i in insert.i) names(cons.seq)[i] <- names(cons.seq)[i - 1]
    }
    
    list(cons.seq = .insertNs(cons.seq), plp = plp, insertions = insertions)
  }, simplify = F, USE.NAMES = T)
  
  cons.seq <- sapply(results, function(x) x$cons.seq, simplify = F)
  plp <- do.call(rbind, lapply(results, function(x) x$plp))
  ins.list <- lapply(results, function(x) x$insertions)
  ins.list <- ins.list[!sapply(ins.list, is.null)]
  insertions <- do.call(rbind, ins.list)
  
  list(cons.seq = cons.seq, plp = plp, insertions = insertions)
}