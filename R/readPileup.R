readPileup <- function(fname) {
  col.names <- c("chrom", "pos", "ref", "cov", "bases", "quals")
  plp <- data.table::fread(fname, sep = "\t")
  colnames(plp)[1:6] <- col.names
  plp %>%
    dplyr::mutate(
      ref = toupper(ref),
      bases = toupper(bases)
    ) %>%
    select(all_of(col.names))
}