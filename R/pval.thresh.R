pval.thresh <- function(pvals, p.thresholds, beta, bfile, 
                        keep=NULL, remove=NULL, extract=NULL, exclude=NULL, 
                        chr=NULL, cluster=NULL, trace=0) {
  #' Fast way to do p-value thresholding (without looping over the thresholds)
  #' @keywords internal
  beta <- as.vector(beta)
  stopifnot(is.vector(pvals) & is.vector(beta))
  stopifnot(length(pvals) == length(beta))
  Pvals <- sort(unique(c(0,p.thresholds,1)))
  cut <- cut(pvals, Pvals, include.lowest = TRUE)
  stopifnot(!any(is.na(cut)))
  nlevels <- nlevels(cut)
  
  parsed <- parseselect(bfile, extract, exclude, keep, remove, chr)
  pbin <- as.integer(cut) - 1
  attr(pbin, "nbin") <- nlevels

  result <- pgs(bfile, keep=parsed$keep, extract=parsed$extract,
             weights = beta, pbin=pbin, cluster=cluster)

  return(result)
  
}