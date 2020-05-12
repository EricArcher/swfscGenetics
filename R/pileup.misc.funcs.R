#' @noRd
.indelFreq <- function(coverage, bases, flag) {
  pattern <- paste("[", flag, "][0-9]+[ACGTNacgtn]+", sep = "")
  pos <- gregexpr(pattern, bases)[[1]]
  if(pos[1] == -1) NULL else {
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
      c(indel, indel.bases)
    })
    rownames(indels) <- c("text", "bases")
    list(
      freq = table(indels["bases", ]) / coverage,
      text = unique(indels["text", ])
    )
  }
}


#' @noRd
.calcFreqs <- function(plp) {
  freqs <- lapply(1:nrow(plp), function(i) {
    bases <- plp$bases[i]
    coverage <- plp$cov[i]
    # remove insertions
    in.freq <- .indelFreq(coverage, bases, "+")
    if(!is.null(in.freq$text)) {
      for(x in in.freq$text) bases <- gsub(x, "", bases, fixed = T)
    }
    # remove deletions
    del.freq <- .indelFreq(coverage, bases, "-")
    if(!is.null(del.freq$text)) {
      for(x in in.freq$text) bases <- gsub(x, "", bases, fixed = T)
    }
    # remove read markers
    bases <- gsub("[ ^][[:alnum:]|[:punct:]]|[$]", "", bases)
    # insert reference base
    bases <- gsub("[.|,]", plp$ref[i], bases)
    # insert deletions
    bases <- gsub("[*]", "-", bases)
    # create factor
    bases <- strsplit(bases, "")[[1]]
    bases <- factor(bases, levels = c("A", "C", "G", "T", "N", "-"))

    list(
      plp = plp[i, c("chrom", "pos", "cov")],
      base.freq = table(bases) / coverage,
      in.freq = in.freq$freq
    )
  })
  stats::setNames(freqs, plp$pos)
}


#' @noRd
.makeSelections <- function(pos.freqs, params) {
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

  # collect called bases
  bases <- sapply(pos.freqs, function(pos) {
    .consensusBase(pos$plp$cov, pos$base.freq, params)
  })

  # collect insertions
  ins <- sapply(pos.freqs, function(pos) {
    ins.base <- .consensusBase(pos$plp$cov, pos$in.freq, params)
    if(is.null(ins.base)) return(NULL)
    if(ins.base == "N") return(NULL)
    strsplit(ins.base, "")[[1]]
  }, simplify = FALSE)
  ins <- ins[!sapply(ins, is.null)]
  if(length(ins) == 0) ins <- NULL

  list(bases = bases, ins = ins)
}


#' @noRd
.plpSummary <- function(plp, cons.seq, pos.freqs) {
  plp$bases <- NULL
  plp$consensus <- cons.seq
  # collect frequency of called base
  plp$cons.freq <- sapply(1:nrow(plp), function(i) {
    called.base <- plp$consensus[i]
    if(called.base == "N") return(NA)
    pos.freqs[[i]]$base.freq[called.base]
  })
  plp$cons.diff <- plp$consensus != plp$ref
  # create table of all base frequencies at position
  base.freq.list <- lapply(pos.freqs, function(pos) pos$base.freq)
  freq.table <- data.frame(do.call(rbind, base.freq.list))
  colnames(freq.table)[6] <- "-"
  freq.table$pos <- as.numeric(names(pos.freqs))
  # format final data frame
  plp <- merge(plp, freq.table, by = "pos", all.x = T)
  plp$cons.pos <- match(as.character(plp$pos), names(cons.seq))
  columns <- c("chrom", "pos", "cons.pos", "ref", "consensus", "cons.diff",
               "cons.freq", "cov", "A", "C", "G", "T", "N", "-")
  plp <- plp[, columns]
  colnames(plp)[c(2, 4, 8)] <- c("ref.pos", "reference", "coverage")
  plp
}


#' @noRd
.insertNs <- function(cons.seq) {
  # create a data.frame of current consensus sequence
  cons.df <- data.frame(
    pos = as.numeric(names(cons.seq)),
    base = cons.seq,
    stringsAsFactors = FALSE
  )
  # column 'i' is used to make sure bases stay in proper order after final sort
  cons.df$i <- 1:nrow(cons.df)
  # find missing sites (reference position numbers not present)
  missing.sites <- sort(setdiff(1:max(cons.df$pos), cons.df$pos))
  # associate missing sites with N's and bind to bottom of consensus data.frame
  if(length(missing.sites) > 0) {
    cons.df <- rbind(
      cons.df,
      data.frame(
        pos = missing.sites,
        base = "N",
        i = 1:length(missing.sites),
        stringsAsFactors = FALSE
      )
    )
  }
  # sort data.frame by position
  cons.df <- cons.df[order(cons.df$pos, cons.df$i), ]
  # return new vector of consensus sequence
  cons.seq <- as.vector(cons.df$base)
  names(cons.seq) <- cons.df$pos
  cons.seq
}