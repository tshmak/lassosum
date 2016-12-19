#' @title Function to convert p-values to correlation via the t-statistic
#' @param n Sample size
#' @param p Vector of p-values
#' @param sign A vector giving the sign of the correlations (e.g. the log odds ratios)
#' @return A vector of correlations
#' @export
p2cor <- function(p, n, sign=rep(1, length(p))) {
  
  stopifnot(length(n)==1 || length(n) == length(p))
  stopifnot(length(p) == length(sign))
  
  t <- sign(sign) * qt(p/2, df=n-2, lower.tail=F)
  
  return(t / sqrt(n - 2 + t^2))
  
}
