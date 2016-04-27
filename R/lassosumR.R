#' @title Function to obtain LASSO estimates of a regression problem given summary statistics
#' and a reference panel (without PLINK bfile) 
#' @details A function to find the minimum of \eqn{\beta} in  
#' \deqn{f(\beta)=\beta'R\beta - 2\beta'r + 2\lambda||\beta||_1}
#' where 
#' \deqn{R=(1-s)X'X/n + sI}
#' is a shrunken correlation matrix, with \eqn{X} being standardized reference panel.
#' \eqn{s} should take values in (0,1]. \eqn{r} is a vector of correlations. 
#' @note \itemize{
#' \item Missing values in \code{refpanel} are filled with 0. 
#' \item Unlike lassosum, we do not provide the options keep/remove/extract/exclude.
#' It is thus up to the user to ensure the SNPs in the reference panel corresponds 
#' to those in the correlations.
#' }
#' @param cor A vector of correlations (\eqn{r})
#' @param refpanel reference panel as \code{data.frame} or \code{matrix}
#' @param lambda A vector of \eqn{\lambda}s (the tuning parameter)
#' @param shrink The shrinkage parameter \eqn{s} for the correlation matrix \eqn{R} 
#' @param thr convergence threshold for \eqn{\beta}
#' @param init Initial values for \eqn{\beta}
#' @param trace An integer controlling the amount of output generated. 
#' @param maxiter Maximum number of iterations
#' @param ridge Produce ridge regression results also (slow if nrow(refpanel) > 2000)
#' 
#' @export

lassosumR <- function(cor, refpanel, 
                     lambda=exp(seq(log(0.001), log(0.1), length.out=20)), 
                     shrink=0.9, ridge=F, 
                     thr=1e-4, init=NULL, trace=0, maxiter=10000) {
  
  stopifnot(is.matrix(refpanel) || is.data.frame(refpanel))
  cor <- as.vector(cor)
  stopifnot(!any(is.na(cor)))
  stopifnot(length(cor) == ncol(refpanel))
  refpanel[is.na(refpanel)] <- 0
	N <- nrow(refpanel)
	p <- ncol(refpanel)
	X <- scale(refpanel) / sqrt(N-1) * sqrt(1-shrink)
	X[is.nan(X)] <- 0
	el <- elnetR(lambda, shrink, X, cor)
	if(ridge) {
	  svd <- svd(X, nu=0)
	  invert <- 1/(svd$d^2 + shrink)
	  VG <- svd$v %*% Diagonal(x=invert - 1/shrink)
	  VTr <- t(svd$v) %*% cor
	  VGVTr <- VG %*% VTr
	  # yTy <- sum(cor * VGVTr) + 1/shrink * sum(cor^2)
	  # rss <- el$loss + shrink * colSums(el$beta^2) + yTy
	  Ridge <- as.vector(VGVTr + 1/shrink * cor)
	} else Ridge <- NULL

	nparams <- colSums(el$beta != 0) 

	return(list(lambda=lambda, 
	            beta=el$beta,
	            conv=el$conv,
	            pred=el$pred,
	            loss=el$loss,
	            fbeta=el$fbeta,
	            sd=attr(X, "scaled:scale"),
	            shrink=shrink,
	            nparams=nparams, 
	            ridge=Ridge))
	#' @return A list with the following
	#' \item{lambda}{same as the lambda input}
	#' \item{beta}{A matrix of estimated coefficients}
	#' \item{conv}{A vector of convergence indicators. 1 means converged. 0 not converged.}
	#' \item{pred}{\eqn{=(1-s)X\beta}}
	#' \item{loss}{\eqn{=(1-s)\beta'X'X\beta/n - 2\beta'r}}
	#' \item{fbeta}{\eqn{=\beta'R\beta - 2\beta'r + 2\lambda||\beta||_1}}
	#' \item{sd}{The standard deviation of the reference panel SNPs}
	#' \item{shrink}{same as input}
	#' \item{nparams}{Number of non-zero coefficients}
	#' \item{ridge}{ridge regression estimates}
	
}
