#' @title Convert MicroHaplot file to rubias format
#' @description Converts a reference or mixture file from the MicroHaplot 
#'   package to the proper format for the rubias package
#'   
#' @param df data frame from MicroHaplot
#' @param sample.type determines type of MicroHaplot sample
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} Converted from code by
#'   John Horne (\url{https://johnbhorne.wordpress.com/2019/03/22/snp-haplotypes-for-mixed-stock-fishery-analysis-microhaplot-to-rubias-conversion-in-r/})
#'
#' @export
#' 
microhaplot2rubias <- function(df, sample.type = c("reference", "mixture")) {
  df <- df[, 1:8] 
  df <- df[order(df$group, df$indiv.ID, df$locus), ]
  
  rubias <- do.call(rbind, lapply(split(df, df$indiv.ID), function(id.df) {
    dups <- duplicated(id.df$locus)
    dups[which(dups) - 1] <- TRUE
    # n.times <- as.numeric(dups)              # count number of duplicates per locus
    # n.times[which(n.times==0)] <- 2 # loci w/ 0 dups are homozygous, still need 2 alleles
    # df3 <- id[rep(seq_len(nrow(id)), n.times),]  # dataframe w/ right number of entries
    # df3$locus <- make.unique(as.character(df3$locus)) # make unique 2nd allele name
    # gt <- spread(df3[,3:5], locus, haplo) # put everything in the right shape
  }))
  colnames(rubias)[1] <- "indiv"
  
  sample.type <- match.arg(sample.type) 
  group <- if(sample.type == "reference") df$group[1] else NA
  cbind(
    data.frame(
      sample_type = rep(sample.type, nrow(rubias)), 
      repunit = as.character(rep(group, nrow(rubias))),           
      collection = as.character(rep(df$group[1], nrow(rubias))),
      stringsAsFactors = FALSE
    ),
    rubias
  )
}
