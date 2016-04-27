#' @title Independent LASSO using summary statistics (a.k.a. soft-thresholding)
#' 
#' @param coef vector of regression coefficients (\eqn{r})
#' @param lambda a vector of \eqn{\lambda}s 
#' 
#' @details A function to find the minimum of \eqn{\beta} in  
#' \deqn{f(\beta)=\beta'\beta - 2\beta'r + 2\lambda||\beta||_1}
#' where \eqn{r} is the vector of regression coefficients. The analytical solution
#' is given by
#' \deqn{\hat{\beta}=sign(r)(max(|r| - \lambda))}
#' @export
indeplasso <- function(coef, lambda=exp(seq(log(0.001), log(0.1), length.out=20))) {
  results <- outer(coef, rep(NA, length(lambda)))
  conv <- rep(0, length(lambda))
  for(i in 1:length(lambda)) {
    results[,i] <- sign(coef) * pmax((abs(coef) - lambda[i]),0)
    conv[i] <- 1
  }
  
  #' @return A list with the following
  #' \item{lambda}{Same as \code{lambda} in input}
  #' \item{beta}{A matrix of estimates of \eqn{\beta}}
  #' \item{conv}{A vector of convergence indicators. 1=converged. 0=not converged.}
  
  return(list(lambda=lambda, beta=results, conv=conv))
}
