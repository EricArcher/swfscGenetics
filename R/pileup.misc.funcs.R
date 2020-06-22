#' @noRd
.indelFreq <- function(coverage, bases, flag) {
  pos <- gregexpr(
    paste("[", flag, "][0-9]+[ACGTNacgtn]+", sep = ""), 
    bases
  )[[1]]
  if(pos[1] == -1) return(NULL)
  indels <- sapply(1:length(pos), function(i) {
    # get text initially matched by regexp
    indel.end <- pos[i] + attr(pos, "match.length")[i] - 1
    indel <- substr(bases, pos[i], indel.end)
    
    # get indel length value
    length.start <- regexpr("[0-9]+", indel)
    length.end <- length.start + attr(length.start, "match.length") - 1
    indel.length <- as.numeric(substr(indel, length.start, length.end))
    
    # get bases that are indel.length characters after number
    bases.start <- pos[i] + 1 + nchar(as.character(indel.length))
    bases.end <- bases.start + indel.length - 1
    indel.bases <- substr(bases, bases.start, bases.end)
    
    # get actual indel text
    indel <- substr(bases, pos[i], bases.end)
    c(text = indel, bases = indel.bases)
  })
  indels <- data.frame(t(indels), stringsAsFactors = FALSE)
  freq <- as.data.frame(
    table(bases = indels$bases), 
    responseName = "freq",
    stringsAsFactors = FALSE
  ) 
  freq$freq <- freq$freq / coverage
  indel.unique <- unique(indels)
  merge(freq, unique(indels), by = "bases", all.x = TRUE)
}


#' @noRd
.removeIndels <- function(bases, text) {
  if(!is.null(text)) for(x in text) bases <- gsub(x, "", bases, fixed = T)
  bases
}


#' @noRd
.freqFunc <- function(i, plp) {
  bases <- plp$bases[i]
  coverage <- plp$cov[i]
  # remove insertions
  ins.freq <- .indelFreq(coverage, bases, "+")
  bases <- .removeIndels(bases, ins.freq$text)
  # remove deletions
  del.freq <- .indelFreq(coverage, bases, "-")
  bases <- .removeIndels(bases, del.freq$text)
  # remove read markers
  bases <- gsub("[ ^][[:alnum:]|[:punct:]]|[$]", "", bases)
  # insert reference base
  bases <- gsub("[.|,]", plp$ref[i], bases)
  # insert deletions
  bases <- gsub("[*]", "-", bases)
  # create factor
  bases <- strsplit(bases, "")[[1]]
  bases <- factor(bases, levels = c("A", "C", "G", "T", "N", "-"))
  
  base.freq <- table(bases = bases) / coverage
  base.freq <- stats::setNames(as.vector(base.freq), names(base.freq))
  ins.freq <- if(!is.null(ins.freq)) {
    stats::setNames(ins.freq$freq, ins.freq$bases)
  }
  
  list(bases = base.freq, ins = ins.freq)
}


#' @noRd
.consensusBase <- function(coverage, freqs, params) {
  if(is.null(freqs)) return(NULL)
  if(coverage < params$min.cov) "N" else {
    bases <- names(freqs)
    if(coverage < params$min.freq.cov) {
      if(all(freqs < 1)) "N" else bases[freqs == 1]
    } else {
      if(!any(freqs >= params$min.freq)) "N" else bases[which.max(freqs)]
    }
  }
}


#' @noRd
.compilePileup <- function(i, plp, params) {
  freqs <- .freqFunc(i, plp)
  plp <-  plp[i, c("chrom", "pos", "ref", "cov", "base.quals")]
  plp$consensus <- .consensusBase(plp$cov, freqs$bases, params)
  plp$cons.freq <- freqs$bases[plp$consensus]
  plp <- cbind(plp, rbind(freqs$bases))
  if(!is.null(freqs$ins)) {
    plp.ins <- plp
    plp.ins$consensus <- .consensusBase(plp$cov, freqs$ins, params)
    plp.ins[, names(freqs$bases)] <- 0
    plp.ins$cons.freq <- plp.ins[, "-"] <- freqs$ins[plp.ins$consensus]
    plp <- rbind(plp, plp.ins)
  }
  plp$pos.order <- 1:nrow(plp)
  list(plp = plp, ins.freq = freqs$ins)
}


#' @noRd
.plpSummary <- function(plp, params, num.cores) {
  # Setup number of cores
  if(is.null(num.cores)) num.cores <- parallel::detectCores() - 1
  if(is.na(num.cores)) num.cores <- 1
  num.cores <- max(1, num.cores)
  num.cores <- min(parallel::detectCores() - 1, num.cores)
  cl <- swfscMisc::setupClusters(num.cores)
  freqs <- tryCatch({
    if(is.null(cl)) {     
      freqs <- lapply(1:nrow(plp), .compilePileup, plp = plp, params = params)
    } else {
      parallel::clusterExport(cl = cl, varlist = "plp", envir = environment())
      parallel::parLapplyLB(
        cl, 1:nrow(plp), .compilePileup, plp = plp, params = params
      )
    }
  }, finally = if(!is.null(cl)) parallel::stopCluster(cl) else NULL)
  
  stats::setNames(freqs, plp$pos)
}


#' @noRd
.compileInsertions <- function(plp.smry) {
  ins.df <- do.call(rbind, lapply(plp.smry, function(x) {
    if(is.null(x$ins.freq)) return(NULL)
    i <- rep(1, length(x$ins.freq))
    plp.ins <- x$plp[i, c("chrom", "pos")]
    plp.ins$insertion <- names(x$ins.freq)
    plp.ins$freq <- x$ins.freq
    plp.ins
  }))
  if(!is.null(ins.df)) {
    rownames(ins.df) <- NULL
    colnames(ins.df)[2] <- "ref.pos"
  }
  ins.df
}


#' @noRd
.compileSummaries <- function(smry, chrom.names) {
  cons.seq <- lapply(smry, function(x) x$cons.seq)
  cons.seq <- stats::setNames(cons.seq, chrom.names)
  
  plp <- do.call(rbind, lapply(smry, function(x) x$plp)) 
  plp <- plp[order(plp$chrom, plp$cons.pos), ]
  columns <- c("chrom", "ref.pos", "pos.order", "cons.pos", "ref", "consensus",
               "cons.freq", "cov", "A", "C", "G", "T", "N", "-")
  plp <- plp[, columns]
  colnames(plp)[c(5, 8)] <- c("reference", "coverage")
  
  insertions <- do.call(rbind, lapply(smry, function(x) x$insertions))
  insertions <- insertions[order(insertions$ref.pos, insertions$freq), ]
  
  list(cons.seq = cons.seq, plp = plp, insertions = insertions)
}