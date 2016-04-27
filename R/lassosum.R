#' @title Function to obtain LASSO estimates of a regression problem given summary statistics
#' and a reference panel. 
#' @details A function to find the minimum of \eqn{\beta} in  
#' \deqn{f(\beta)=\beta'R\beta - 2\beta'r + 2\lambda||\beta||_1}
#' where 
#' \deqn{R=(1-s)X'X/n + sI}
#' is a shrunken correlation matrix, with \eqn{X} being standardized reference panel.
#' \eqn{s} should take values in (0,1]. \eqn{r} is a vector of correlations. 
#' \code{keep}, \code{remove} could take one of three 
#' formats: (1) A logical vector indicating which indivduals to keep/remove, 
#' (2) A \code{data.frame} with two columns giving the FID and IID of the indivdiuals
#' to keep/remove (matching those in the .fam file), or (3) a character scalar giving the text file with the FID/IID. 
#' Likewise \code{extract}, \code{exclude} can also take one of the three formats,
#' except with the role of the FID/IID data.frame replaced with a character vector of 
#' SNP ids (matching those in the .bim file). 
#' @note Missing genotypes are interpreted as having the homozygous A2 alleles in the 
#' PLINK files (same as the \code{--fill-missing-a2} option in PLINK). 
#' @param cor A vector of correlations (\eqn{r})
#' @param bfile PLINK bfile (as character, without the .bed extension)
#' @param lambda A vector of \eqn{\lambda}s (the tuning parameter)
#' @param shrink The shrinkage parameter \eqn{s} for the correlation matrix \eqn{R} 
#' @param thr convergence threshold for \eqn{\beta}
#' @param init Initial values for \eqn{\beta}
#' @param trace An integer controlling the amount of output generated. 
#' @param maxiter Maximum number of iterations
#' @param extract samples to extract
#' @param exclude samples to exclude
#' @param keep SNPs to keep
#' @param remove SNPs to remove
#' @param chr a vector of chromosomes
#' 
#' @export

lassosum <- function(cor, bfile, 
                     lambda=exp(seq(log(0.001), log(0.1), length.out=20)), 
                     shrink=0.9, 
                     thr=1e-4, init=NULL, trace=0, maxiter=10000, 
                     keep=NULL, extract=NULL, exclude=NULL, remove=NULL, 
                     chr=NULL) {
  
  stopifnot(is.numeric(cor))
  stopifnot(!any(is.na(cor)))
  if(any(abs(cor) > 1)) warning("Some abs(cor) > 1")
  if(any(abs(cor) == 1)) warning("Some abs(cor) == 1")
  
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr)

  if(length(cor) != parsed$p) stop("Length of cor does not match number of selected columns in bfile")
  # stopifnot(length(cor) == parsed$p)
  
  if(is.null(parsed$extract)) {
    extract2 <- list(integer(0), integer(0))
  } else {
	  # print(parsed$extract)
    extract2 <- selectregion(!parsed$extract)
    extract2[[1]] <- extract2[[1]] - 1
  }
  
  if(is.null(parsed$keep)) {
    keepbytes <- integer(0)
    keepoffset <- integer(0)
  } else {
    pos <- which(parsed$keep) - 1
    keepbytes <- floor(pos/4)
    keepoffset <- pos %% 4 * 2
  }
  
  if(is.null(init)) init <- rep(0.0, parsed$p) else {
    stopifnot(is.numeric(init) && length(init) == parsed$p)
  }
 # print(extract2[[1]])
 # print(extract2[[2]])
 # print(4000-sum(extract2[[2]]))
  
  init <- init + 0.0 # force R to create a copy
  
  order <- order(lambda, decreasing = T)

  results <- runElnet(lambda[order], shrink, fileName=paste0(bfile,".bed"), 
                      r=cor, N=parsed$N, P=parsed$P, 
                      col_skip_pos=extract2[[1]], col_skip=extract2[[2]],
                      keepbytes=keepbytes, keepoffset=keepoffset, 
                      thr=1e-4, x=init, trace=trace, maxiter=maxiter)
  results <- within(results, {
    conv[order] <- conv
    beta[,order] <- beta
    pred[,order] <- pred
    loss[order] <- loss
    fbeta[order] <- fbeta
    lambda[order] <- lambda
  })
  results$shrink <- shrink
  results$nparams <- colSums(results$beta != 0)
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
  
  return(results)

}
