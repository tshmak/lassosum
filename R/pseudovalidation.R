#' @title Performs `pseudovalidation' to select the best \eqn{\lambda} value in lassosum
#' 
#' @param bfile A plink bfile stem
#' @param beta The matrix of estimated \eqn{\beta}s
#' @param cor The vector of correlations (\eqn{r})
#' @param sd The standard deviation of the SNPs
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
#' @param cluster A \code{cluster} object from the \code{parallel} package for parallel computing
#' @details A function to calculate  
#' \deqn{f(\lambda)=\beta'r/\sqrt{\beta'X'X\beta}} 
#' where \eqn{X} is the standardized genotype matrix divided by \eqn{\sqrt n}, 
#' and \eqn{r} is a vector of (shrunken) correlations. 
#' @note \itemize{
#' \item Missing genotypes are interpreted as having the homozygous A2 alleles in the 
#' PLINK files (same as the \code{--fill-missing-a2} option in PLINK). 
#' \item The number of rows in \code{beta} and the length of \code{cor} should be the same as the number of
#' SNPs in the bfile after extract/exclude/chr.
#' }
#' @keywords internal
pseudovalidation <- function(bfile, beta, cor, sd=NULL, 
                             keep=NULL, extract=NULL, exclude=NULL, remove=NULL, 
                             chr=NULL, cluster=NULL, ...) {

  stopifnot(is.numeric(cor))
  stopifnot(!any(is.na(cor)))
  if(any(abs(cor) > 1)) warning("Some abs(cor) > 1")
  if(any(abs(cor) == 1)) warning("Some abs(cor) == 1")

  beta <- as.matrix(beta)
  stopifnot(!any(is.na(beta)))
  if(length(cor) != nrow(beta)) stop("Length of cor does not match number of rows in beta")
  
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr)
  if(length(cor) != parsed$p) stop("Length of cor does not match number of selected columns in bfile")
  
  
  if(is.null(sd)) sd <- sd.bfile(bfile = bfile, 
                                 keep=parsed$keep, extract=parsed$extract, ...)
  stopifnot(length(sd) == length(cor))

    weight <- 1/sd
    weight[!is.finite(weight)] <- 0
    scaled.beta <- as.matrix(Diagonal(x=weight) %*% beta)
    pred <- pgs(bfile, keep=parsed$keep, extract=parsed$extract, 
                weights=scaled.beta, cluster=cluster)
    pred2 <- scale(pred, scale=F)
    bXXb <- colSums(pred2^2) / parsed$n
    bXy <- cor %*% beta 
    
    result <- as.vector(bXy / sqrt(bXXb))
    attr(result, "bXXb") <- bXXb
    attr(result, "bXy") <- bXy
    
    return(result)
    #' @return the results of the pseudovalidation, i.e. \eqn{f(\lambda)}
    
}
