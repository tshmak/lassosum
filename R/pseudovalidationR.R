#' @title Performs `pseudovalidation' to select the best \eqn{\lambda} value in 
#' lassosum (without PLINK bfile)
#' 
#' @param genotype.mat A genotype matrix (coded 0/1/2)
#' @param beta The matrix of estimated \eqn{\beta}s
#' @param cor The vector of correlations (\eqn{r})
#' @details A function to calculate  
#' \deqn{f(\lambda)=\beta'r/\sqrt{\beta'X'X\beta}} 
#' where \eqn{X} is the standardized genotype matrix divided by \eqn{\sqrt n}, 
#' and \eqn{r} is a vector of (shrunken) correlations. 
#' @note \itemize{
#' \item The number of rows in \code{beta} and the length of \code{cor} should be the 
#' same as the number of columns in \code{genotype.mat}.
#' }
#' @keywords internal
pseudovalidationR <- function(genotype.mat, beta, cor) {

  stopifnot(!any(is.na(genotype.mat)))
  stopifnot(is.numeric(cor))
  stopifnot(!any(is.na(cor)))
  if(any(abs(cor) > 1)) warning("Some abs(cor) > 1")
  if(any(abs(cor) == 1)) warning("Some abs(cor) == 1")

  beta <- as.matrix(beta)
  stopifnot(!any(is.na(beta)))
  if(length(cor) != nrow(beta)) stop("Length of cor does not match number of rows in beta")
  if(ncol(genotype.mat) != nrow(beta)) stop("Number of columns in genotype.mat does not match number of rows in beta")
  
  X <- scale(genotype.mat) / sqrt(nrow(genotype.mat)-1)
  bXy <- as.vector(crossprod(beta, cor))
  Xb <- X %*% beta
  bXXb <- as.vector(colSums(Xb^2))

  result <- as.vector(bXy / sqrt(bXXb))
  attr(result, "bXXb") <- bXXb
  attr(result, "bXy") <- bXy
  
  return(result)
  #' @return the results of the pseudovalidation, i.e. \eqn{f(\lambda)}
    
}
