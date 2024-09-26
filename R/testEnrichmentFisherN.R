testEnrichmentFisherN <- function(
    nD, nQ, nDQ, nU, alternative = "greater") {
  #N_mask,N_query,N_overlap,N_universe
  nDmQ <- nD - nDQ
  nQmD <- nQ - nDQ
  nUmDQ <- nU - nQ - nD + nDQ
  
  if (alternative == "two.sided") {
    pvg <- phyper(
      nDQ-1, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
      lower.tail = FALSE, log.p = TRUE) / log(10)
    pvl <- phyper(
      nDQ, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
      lower.tail = TRUE, log.p = TRUE) / log(10)
    log10.p.value <- pmin(pmin(pvg, pvl) + log(2), 0) / log(10)
    ## log10.p.value <- log10(fisher.test(matrix(c(
    ##     nDQ, nDmQ, nQmD, nUmDQ), nrow = 2))$p.value)
  } else if (alternative == "greater") {
    log10.p.value <- phyper(
      nDQ-1, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
      lower.tail = FALSE, log.p = TRUE) / log(10)
  } else if (alternative == "less") {
    log10.p.value <- phyper(
      nDQ, nDQ + nQmD, nUmDQ + nDmQ, nDmQ + nDQ,
      lower.tail = TRUE, log.p = TRUE) / log(10)
  } else {
    stop("alternative must be either greater, less or two-sided.")
  }
  
  odds_ratio <- nDQ / nQmD / nDmQ * nUmDQ # can be NaN if 0
  odds_ratio[odds_ratio == Inf] <- .Machine$double.xmax
  odds_ratio[odds_ratio == 0] <- .Machine$double.xmin
  data.frame(
    estimate = log2(odds_ratio),
    p.value = 10**(log10.p.value),
    log10.p.value = log10.p.value,
    test = "Log2(OR)",
    nQ = nQ, nD = nD, overlap = nDQ,
    cf_Jaccard = nDQ / (nD + nQmD),
    cf_overlap = nDQ / pmin(nD, nQ), # Szymkiewicz–Simpson
    cf_NPMI = (log2(nD)+log2(nQ)-2*log2(nU))/(log2(nDQ)-log2(nU))-1,
    cf_SorensenDice = 2 * nDQ/(nD + nQ))
}