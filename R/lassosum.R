#' @title Function to obtain LASSO estimates of a regression problem given summary statistics
#' and a reference panel
#' 
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
#' 
#' @note Missing genotypes are interpreted as having the homozygous A2 alleles in the 
#' PLINK files (same as the \code{--fill-missing-a2} option in PLINK). 
#' 
#' @param cor A vector of correlations (\eqn{r})
#' @param bfile PLINK bfile (as character, without the .bed extension)
#' @param lambda A vector of \eqn{\lambda}s (the tuning parameter)
#' @param shrink The shrinkage parameter \eqn{s} for the correlation matrix \eqn{R} 
#' @param thr convergence threshold for \eqn{\beta}
#' @param init Initial values for \eqn{\beta} as a vector of the same length as \code{cor}
#' @param trace An integer controlling the amount of output generated. 
#' @param maxiter Maximum number of iterations
#' @param blocks A vector to split the genome by blocks (coded as c(1,1,..., 2, 2, ..., etc.))
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
#' @param mem.limit Memory limit for genotype matrix loaded. Note that other overheads are not included. 
#' @param chunks Splitting the genome into chunks for computation. Either an integer 
#' indicating the number of chunks or a vector (length equal to \code{cor}) giving the exact split. 
#' @param cluster A \code{cluster} object from the \code{parallel} package for parallel computing
#' 
#' @export

lassosum <- function(cor, bfile, 
                     lambda=exp(seq(log(0.001), log(0.1), length.out=20)), 
                     shrink=0.9, 
                     thr=1e-4, init=NULL, trace=0, maxiter=10000, 
                     blocks=NULL,
                     keep=NULL, remove=NULL, extract=NULL, exclude=NULL, 
                     chr=NULL, 
                     mem.limit=4*10^9, chunks=NULL, cluster=NULL) {
  
  stopifnot(is.numeric(cor))
  stopifnot(!any(is.na(cor)))
  if(any(abs(cor) > 1)) warning("Some abs(cor) > 1")
  if(any(abs(cor) == 1)) warning("Some abs(cor) == 1")
  if(length(shrink) > 1) stop("Only 1 shrink parameter at a time.")
  
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr)
  if(is.null(blocks)) {
    Blocks <- list(startvec=0, endvec=parsed$p - 1)
  } else {
    Blocks <- parseblocks(blocks)
    stopifnot(max(Blocks$endvec)==parsed$p - 1)
  }

  if(length(cor) != parsed$p) stop("Length of cor does not match number of selected columns in bfile")
  # stopifnot(length(cor) == parsed$p)
  
  #### Group blocks into chunks ####
  chunks <- group.blocks(Blocks, parsed, mem.limit, chunks, cluster)
  if(trace > 0) {
    if(trace - floor(trace) > 0) {
      cat("Doing lassosum on chunk", unique(chunks$chunks), "\n")
    } else {
      cat("Calculations carried out in ", max(chunks$chunks.blocks), " chunks\n")
    }
  }
  if(length(unique(chunks$chunks.blocks)) > 1) {
    if(is.null(cluster)) {
      results.list <- lapply(unique(chunks$chunks.blocks), function(i) {
        lassosum(cor=cor[chunks$chunks==i], bfile=bfile, lambda=lambda, shrink=shrink, 
                 thr=thr, init=init[chunks$chunks==i], trace=trace-0.5, maxiter=maxiter, 
                 blocks[chunks$chunks==i], keep=parsed$keep, extract=chunks$extracts[[i]], 
                 mem.limit=mem.limit, chunks=chunks$chunks[chunks$chunks==i])
      })
    } else {
      Cor <- cor; Bfile <- bfile; Lambda <- lambda; Shrink=shrink; Thr <- thr; 
      Maxiter=maxiter; Mem.limit <- mem.limit ; Trace <- trace; Init <- init; 
      Blocks <- blocks
      # Make sure these are defined within the function and so copied to 
      # the child processes
      results.list <- parallel::parLapplyLB(cluster, unique(chunks$chunks.blocks), function(i) {
        lassosum(cor=Cor[chunks$chunks==i], bfile=Bfile, lambda=Lambda, 
                 shrink=Shrink, thr=Thr, init=Init[chunks$chunks==i], 
                 trace=trace-0.5, maxiter=Maxiter, 
                 blocks=Blocks[chunks$chunks==i], 
                 keep=parsed$keep, extract=chunks$extracts[[i]], 
                 mem.limit=Mem.limit, chunks=chunks$chunks[chunks$chunks==i])
      })
    }
    return(do.call("merge.lassosum", results.list))
  }

  #### Group blocks into chunks 
  
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
                      thr=1e-4, x=init, trace=trace, maxiter=maxiter,
                      startvec=Blocks$startvec, endvec=Blocks$endvec)
  results$sd <- as.vector(results$sd)
  results <- within(results, {
    conv[order] <- conv
    beta[,order] <- beta
    pred[,order] <- pred
    loss[order] <- loss
    fbeta[order] <- fbeta
    lambda[order] <- lambda
  })
  results$shrink <- shrink
  
  if(length(lambda) > 0) results$nparams <- as.vector(colSums(results$beta != 0)) else 
    results$nparams <- numeric(0)
  results$conv <- as.vector(results$conv)
  results$loss <- as.vector(results$loss)
  results$fbeta <- as.vector(results$fbeta)
  results$lambda <- as.vector(results$lambda)
  
  class(results) <- "lassosum"
  return(results)
  #' @return A list with the following
  #' \item{lambda}{same as the lambda input}
  #' \item{beta}{A matrix of estimated coefficients}
  #' \item{conv}{A vector of convergence indicators. 1 means converged. 0 not converged.}
  #' \item{pred}{\eqn{=\sqrt(1-s)X\beta}}
  #' \item{loss}{\eqn{=(1-s)\beta'X'X\beta/n - 2\beta'r}}
  #' \item{fbeta}{\eqn{=\beta'R\beta - 2\beta'r + 2\lambda||\beta||_1}}
  #' \item{sd}{The standard deviation of the reference panel SNPs}
  #' \item{shrink}{same as input}
  #' \item{nparams}{Number of non-zero coefficients}
  
  
}
