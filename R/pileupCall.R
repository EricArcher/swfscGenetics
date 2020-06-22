#' @rdname pileupCall
#' @title Call bases from pileup file
#' @description Reads a pileup formatted file created by \code{samtools mpileup}
#'   and calls bases for each chromosome listed. Base calling is controlled
#'   by coverage and frequency parameters
#'
#' @param folder folder containing pileup files from a run
#' @param fname filename of pileup file
#' @param min.cov.call minimum coverage for base calling. Sites with coverage below
#'   this are assigned \code{N}'s.
#' @param min.cov.freq minimum coverage above which \code{min.freq} is applied. 
#'   Sites below this value and >= than \code{min.cov} will only be called 
#'   if all reads agree.
#' @param min.base.freq minimum frequency of either the reference or alternate base
#'   for calling. If both bases are below this frequency, an \code{N} 
#'   is assigned.
#' @param min.ins.freq minimum frequency of insertion.
#' @param min.prob.freq minimum frequency for binomial probability.
#' @param min.binom.prob minimum probability from binomial distribution.
#' @param pattern text pattern for pileup files.
#' @param label label for run output files.
#' @param num.cores number of cores to use during processing. If \code{NULL}, 
#'   will default to \code{parallel::detectCores() - 1}.
#'
#' @return list with the following elements:
#' \tabular{ll}{ 
#'   \code{plp} \tab data frame.\cr 
#'   \code{base.lo} \tab data frame of bases and log-odds.\cr 
#' }
#'
#' @note The input pileup file should be the result of a call to 
#'   \code{samtools mpileup} on a single BAM file.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#' 
pileupCallRun <- function(min.cov.call, min.cov.freq, 
                          min.base.freq, min.ins.freq, 
                          min.prob.freq, min.binom.prob,
                          folder = ".", pattern = "\\.pileup$",
                          label = NULL, num.cores = NULL) {
  files <- dir(folder, pattern = pattern, full.names = TRUE)
  message("Processing ", length(files), " files")
  for(f in files) {
    out.fname <- paste0("cons.plp_", basename(f), ".csv")
    out.fname <- file.path(folder, out.fname)
    if(file.exists(out.fname)) next
    plp <- pileupCallFile(f, min.cov.call, min.cov.freq, min.base.freq,
                          min.ins.freq, num.cores)
    data.table::fwrite(plp$plp, file = out.fname)
  }
  
  plp <- .binomialBaseCall(
    dir(folder, pattern = "cons.plp_", full.names = TRUE),
    min.prob.freq, min.binom.prob
  )
  
  cons.seq <- sapply(
    split(plp, list(fname = plp$fname, chrom = plp$chrom)), 
    function(x) {
      x <- x$cons.base
      x[!is.na(x)]
    },
    simplify = FALSE
  )
  
  if(is.null(label)) label <- "all_consensus"
  label <- paste0(label, format(Sys.time(), "_%Y%m%d_%H%M"))
  
  ape::write.dna(
    ape::as.DNAbin(cons.seq), 
    file = paste0(label, ".fasta"), 
    format = "fasta", 
    nbcol = -1, 
    colsep = "", 
    indent = 0, 
    blocksep = 0
  )
  
  data.table::fwrite(plp, paste0(label, ".csv"))
  
  list(cons.seq = cons.seq, plp = plp)
}

#' @rdname pileupCall
#' @export
#' 
pileupCallFile <- function(fname, min.cov.call, min.cov.freq, 
                           min.base.freq, min.ins.freq, 
                           num.cores = NULL) {
  # convert parameters and make sure the values are in proper range
  params <- list(
    min.cov.call = min.cov.call,
    min.cov.freq = min.cov.freq,
    min.base.freq = min(max(0.5, as.numeric(min.base.freq)), 1),
    min.ins.freq = min(max(0.5, as.numeric(min.ins.freq)), 1)
  )
  
  message("Reading: ", basename(fname))
  plp.file <- pileupRead(fname)
  plp.file <- split(plp.file, plp.file$chrom)
  
  smry <- lapply(1:length(plp.file), function(i) {
    message(
      "   Processing pileup for chromosome '", names(plp.file)[[i]], "' ",
      "(", i, "/", length(plp.file), ")..."
    ) 
    
    # call bases and summarize frequencies
    suppressPackageStartupMessages(
      .plpSummary(plp.file[[i]], params, num.cores)
    )
  })
  smry <- stats::setNames(smry, names(plp.file))
  
  list(
    cons.seq = sapply(smry, function(x) {
      x <- x$cons.base
      x[!is.na(x)]
    }, simplify = FALSE),
    plp = cbind(fname = basename(fname), dplyr::bind_rows(smry))
  )
}

#' @noRd
#' 
.plpSummary <- function(plp, params, num.cores) {
  # Setup number of cores
  if(is.null(num.cores)) num.cores <- parallel::detectCores() - 1
  if(is.na(num.cores)) num.cores <- 1
  num.cores <- max(1, num.cores)
  num.cores <- min(parallel::detectCores() - 1, num.cores)
  cl <- swfscMisc::setupClusters(num.cores)
  smry <- tryCatch({
    if(is.null(cl)) {     
      freqs <- lapply(1:nrow(plp), .compilePileup, plp = plp, params = params)
    } else {
      parallel::clusterExport(cl = cl, varlist = "plp", envir = environment())
      parallel::parLapplyLB(
        cl, 1:nrow(plp), .compilePileup, plp = plp, params = params
      )
    }
  }, finally = if(!is.null(cl)) parallel::stopCluster(cl) else NULL)
  dplyr::bind_rows(smry)
}

#' @noRd
#' 
.compilePileup <- function(i, plp, params) {
  # message(plp$chrom[i], " ", i)
  bases <- plp$bases[i]
  # delete first position markers
  bases <- gsub("[\\^][[:graph:]]{1}", "", bases)
  # delete last position markers
  bases <- gsub("[\\$]{1}", "", bases)
  # remove deletion text
  bases <- .removeIndels(bases, "-")
  # change deletion padding markers to "-"
  bases <- gsub("\\*", "-", bases)
  
  insertions <- .compileInsertions(bases, plp$cov[i], params)

  # remove insertion text 
  bases <- .removeIndels(bases, "+")
  # replace sites with reference base
  bases <- gsub("[.|,]", plp$ref.base[i], bases)
  # convert bases to single element factor
  bases <- unlist(strsplit(bases, ""))
  bases <- factor(bases, levels = c("A", "C", "G", "T", "N", "-", "+"))
  
  #if(all(is.na(bases))) stop()
  
  base.freq <- if(all(is.na(bases))) table(bases) else {
    # extract quality probabilities and convert to log-odds
    lo <- swfscMisc::logOdds(qual2prob(plp$base.quals[i]))
    # convert log-odds to weights with minimum = 1
    wt <- lo - min(lo) + 1
    # weighted frequencies of bases
    base.freq <- tapply(wt, bases, sum)
    base.freq[is.na(base.freq)] <- 0
    prop.table(base.freq)
  }
  
  # call bases
  base <- if(plp$cov[i] < params$min.cov.call) {
    "1" 
  } else if(plp$cov[i] < params$min.cov.freq) {
    x <- names(base.freq)[base.freq == 1]
    if(length(x) == 0) "2" else x
  } else {
    max.base <- which.max(base.freq)
    if(base.freq[max.base] < params$min.base.freq) "3" else {
      names(base.freq)[max.base]
    }
  }
  
  df <- dplyr::bind_rows(
    cbind(
      data.frame(
        cons.pos = 0,
        cons.base = if(base %in% 1:3) "N" else base,
        n.code = if(base %in% 1:3) base else NA
      ),
      rbind(base.freq)
    ),
    insertions
  )
  
  cbind(plp[rep(i, nrow(df)), c("chrom", "ref.pos", "ref.base", "cov")], df)
  
  # list(    
  #   df = cbind(
  #     plp[rep(i, nrow(df)), c("chrom", "ref.pos", "ref.base", "cov")],
  #     df
  #   ),
  #   base.lo = data.frame(
  #     chrom = plp$chrom[i], 
  #     ref.pos = plp$ref.pos[i],
  #     base = bases, 
  #     log.odds = lo
  #   )
  # )
}

#' @noRd
#' 
.removeIndels <- function(bases, flag) {
  for(x in .indelText(bases, flag)) {
    if(flag == "+") x <- paste0("\\+", x)
    bases <- gsub(x, "", bases)
  }
  bases
}

#' @noRd
#' 
.indelText <- function(bases, flag) {
  num.pattern <- if(flag == "+") "\\+[0-9]+" else "-[0-9]+"
  base.pattern <- if(flag == "+") "[ACGTNacgtn*#]+" else "[ACGTNacgtn]+"
  matches <- gregexpr(paste0(num.pattern, base.pattern), bases) 
  match.text <- unlist(regmatches(bases, matches))
  sapply(match.text, function(x) {
    x.num <- as.numeric(gsub(base.pattern, "", x))
    x.bases <- gsub(num.pattern, "", x)
    x.bases <- substr(x.bases, 1, abs(x.num))
    paste0(x.num, x.bases)
  })
}

#' @noRd
#' 
.compileInsertions <- function(bases, cov, params) {
  # extract insertion text
  ins.text <- .indelText(bases, "+")
  # extract inserted bases
  ins.bases <- unlist(
    regmatches(ins.text, gregexpr("[ACGTNacgtn*#]+", ins.text))
  )
  
  if(is.null(ins.bases)) NULL else {
    ins.freq <- table(ins.bases) / cov
    i <- which.max(ins.freq)
    if(ins.freq[i] < params$min.ins.freq) NULL else {
      data.frame(
        cons.pos = 1, 
        cons.base = names(ins.freq)[i], 
        '+' = ins.freq[i],
        check.names = FALSE
      )
    }
  }
       
  # # compute insertion frequencies relative to coverage  
  # ins.mat <- if(is.null(ins.bases)) NULL else {
  #   x <- strsplit(toupper(ins.bases), "")
  #   mat <- matrix(nrow = length(x), ncol = max(sapply(x, length)))
  #   for(row.i in 1:length(x)) mat[row.i, 1:length(x[[row.i]])] <- x[[row.i]]
  #   mat
  # }
  # ins.freq <- if(is.null(ins.mat)) NULL else {
  #   freq <- apply(
  #     mat, 2, 
  #     function(col.x) table(factor(col.x, levels = c("A", "C", "G", "T", "N")))
  #   )
  #   t(freq) / cov
  # }
  # 
  # # select insertions
  # insertion <- if(is.null(ins.freq)) NULL else {
  #   apply(ins.freq, 1, function(x) {
  #     max.i <- which.max(x)
  #     if(x[max.i] < params$min.ins.freq) NA else colnames(ins.freq)[max.i]
  #   })
  # }
  # 
  # if(is.null(insertion)) NULL else {
  #   cbind(
  #     data.frame(
  #       cons.pos = 1:length(insertion),
  #       cons.base = insertion
  #     ),
  #     ins.freq
  #   )
  # }
}

#' @noRd
#' 
.binomialBaseCall <- function(cons.plp.fnames, min.prob.freq, min.binom.prob) {
  run.plp <- do.call(rbind, lapply(cons.plp.fnames, data.table::fread))
  base.vec <- c("A", "C", "G", "T", "N", "-")
  i <- which(run.plp$n.code == 3)
  for(pos in unique(run.plp$ref.pos[i])) {
    pool.prop <- run.plp %>% 
      dplyr::filter(.data$ref.pos == pos) %>% 
      dplyr::select(dplyr::all_of(c("cov", base.vec))) %>% 
      tidyr::pivot_longer(-.data$cov, names_to = "base", values_to = "freq") %>% 
      dplyr::mutate(freq = round(.data$freq * .data$cov)) %>% 
      dplyr::group_by(.data$base) %>% 
      dplyr::summarize(prob = sum(.data$freq) / sum(.data$cov), .groups = "drop") %>% 
      tibble::deframe()
    pool.prop <- pool.prop[base.vec]
    
    for(j in which(run.plp$n.code == 3 & run.plp$ref.pos == pos)) {
      cov <- run.plp$cov[j]
      base.freq <- round(unlist(run.plp[j, base.vec]) * cov)
      read.prop <- run.plp[j, base.vec]
      
      binom.prob <- sapply(names(pool.prop), function(x) { 
        stats::pbinom(base.freq[x], cov, pool.prop[x])
      })
      cond.1 <- read.prop > pool.prop & pool.prop > 0.5
      cond.2 <- binom.prob >= min.binom.prob & 
        read.prop >= min.prob.freq &
        (0.875 * read.prop) - (pool.prop + 0.175) > 0 
      is.good <- cond.1 | cond.2
      
      if(any(is.good)) {
        good.prop <- read.prop[which(is.good)]
        max.prop <- which.max(good.prop)
        run.plp$cons.base[j] <- names(good.prop)[max.prop]
        run.plp$n.code[j] <- NA
      } else {
        run.plp$n.code[j] <- 4
      }
    }
  }
  run.plp
}