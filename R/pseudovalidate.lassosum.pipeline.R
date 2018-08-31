pseudovalidate.lassosum.pipeline <- function(ls.pipeline, test.bfile=NULL, 
                                       keep=NULL, remove=NULL, 
                                       trace=1, 
                                       destandardize=F, plot=T, 
                                       exclude.ambiguous=T, 
                                       cluster=NULL, 
                                       rematch=!is.null(test.bfile), 
                                       ...) {
  #' @title Function to perform pseudovalidation from a lassosum.pipeline object
  #' @param ls.pipeline A lassosum.pipeline object
  #' @param test.bfile The (\href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK bfile} for the test dataset 
  #' @param keep Participants to keep (see \code{\link{lassosum}} for more details)
  #' @param remove Participants to remove
  #' @param trace Controls amount of output
  #' @param destandardize Should coefficients from \code{\link{lassosum}} be 
  #' destandardized using test dataset standard deviations before being returned?
  #' @param plot Should the validation plot be plotted? 
  #' @param exclude.ambiguous Should ambiguous SNPs (C/G, A/T) be excluded? 
  #' @param cluster A \code{cluster} object from the \code{parallel} package for parallel computing
  #' @param rematch Forces a rematching of the ls.pipline beta's with the new .bim file
  #' @param ... parameters to pass to \code{\link{pseudovalidation}}
  #' @details Pseudovalidation is explained in Mak et al (2016). It helps 
  #' choosing a value of \code{lambda} and \code{s} in the absence of a validation
  #' phenotype. 
  #' @rdname pseudovalidate
  #' @export
  installed <- installed.packages()[,1]
  if(!("fdrtool" %in% installed)) 
    stop("Pseudovalidation requires fdrtool. Please install from CRAN.")
  
  stopifnot(class(ls.pipeline) == "lassosum.pipeline")

  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  
  if(!is.null(keep) || !is.null(remove)) if(is.null(test.bfile))
    stop("Please specify test.bfile if you specify keep or remove")
  
  rematch <- rematch # Forces an evaluation at this point
  if(is.null(test.bfile)) {
    test.bfile <- ls.pipeline$test.bfile
    if(is.null(keep) && is.null(remove))
      keep <- ls.pipeline$keep.test
  }
  
  if(destandardize) {
    if(ls.pipeline$destandardized) stop("beta in ls.pipeline already destandardized.")
    sd <- sd.bfile(test.bfile, extract=ls.pipeline$test.extract, 
                   keep=keep, remove=remove, trace=trace, ...)
    sd[sd <= 0] <- Inf # Do not want infinite beta's!
    ls.pipeline$beta <- lapply(ls.pipeline$beta, 
                   function(x) as.matrix(Matrix::Diagonal(x=1/sd) %*% x))
    recal <- T
  }
  
  parsed.test <- parseselect(test.bfile, keep=keep, remove=remove)
  recal <- !identical(ls.pipeline$test.bfile, test.bfile) || 
    !identical(parsed.test$keep, ls.pipeline$keep.test)
  
  if(rematch) {
    if(trace) cat("Coordinating lassosum output with test data...\n")
    
    if(length(test.bfile) > 1) stop("Multiple 'test.bfile's not supported here.")
    bim <- fread(paste0(test.bfile, ".bim"))
    bim$V1 <- as.character(sub("^chr", "", bim$V1, ignore.case = T))
    
    m <- matchpos(ls.pipeline$sumstats, bim, auto.detect.ref = F, 
                  ref.chr = "V1", ref.snp="V2", ref.pos="V4", ref.alt="V5", ref.ref="V6", 
                  rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                  silent=T)
    
    beta <- lapply(ls.pipeline$beta, function(x) 
      as.matrix(Matrix::Diagonal(x=m$rev) %*% x[m$order, ]))
    if(trace) cat("Calculating PGS...\n")
    toextract <- m$ref.extract
    pgs <- lapply(beta, function(x) pgs(bfile=test.bfile, weights = x, 
                                        extract=toextract, 
                                        keep=parsed.test$keep, 
                                        cluster=cluster,
                                        trace=trace-1))
    names(pgs) <- as.character(ls.pipeline$s)
    results <- c(results, list(pgs=pgs))
    
  } else {
    toextract <- ls.pipeline$test.extract
    if(is.null(ls.pipeline$pgs) || recal) {
      if(trace) cat("Calculating PGS...\n")
      if(length(test.bfile) > 1) stop("Multiple 'test.bfile's not supported here.")
      pgs <- lapply(ls.pipeline$beta, function(x) pgs(bfile=test.bfile, 
                                                      weights = x, 
                                                      keep=parsed.test$keep, 
                                                      extract=ls.pipeline$test.extract, 
                                                      cluster=cluster, 
                                                      trace=trace-1))
      names(pgs) <- as.character(ls.pipeline$s)
      results <- c(results, list(pgs=pgs))
    # } else if(is.null(parsed.test$keep)) {
    } else  {
      results <- c(results, list(pgs=ls.pipeline$pgs))
    # } else {
    #   pgs <- ls.pipeline$pgs
    #   for(i in 1:length(pgs)) {
    #     pgs[[i]] <- pgs[[i]][parsed.test$keep, ]
    #   }
    #   results <- c(results, list(pgs=pgs))
    }
    beta <- ls.pipeline$beta
  } 
  
  ### Pseudovalidation ###
  lambdas <- rep(ls.pipeline$lambda, length(ls.pipeline$s))
  ss <- rep(ls.pipeline$s, rep(length(ls.pipeline$lambda), length(ls.pipeline$s)))
  PGS <- do.call("cbind", results$pgs)
  BETA <- do.call("cbind", ls.pipeline$beta)
  
  if(trace) cat("Estimating local fdr ...\n")
  fdr <- fdrtool::fdrtool(ls.pipeline$sumstats$cor, statistic="correlation", 
                          plot=F)
  cor.shrunk <- ls.pipeline$sumstats$cor * (1 - fdr$lfdr)
  if(trace) cat("Performing pseudovalidation ...\n")
  pv <- pseudovalidation(test.bfile, 
                         beta=BETA, 
                         cor=cor.shrunk, 
                         extract=toextract, 
                         keep=keep, remove=remove,
                         sd=ls.pipeline$sd, 
                         cluster=cluster, ...)

  pv[is.na(pv)] <- -Inf
  best <- which(pv == max(pv))[1]
  best.s <- ss[best]
  best.lambda <- lambdas[best]
  best.pgs <- PGS[,best]
  len.lambda <- length(ls.pipeline$lambda)
  best.beta.s <- ceiling(best / len.lambda)
  best.beta.lambda <- best %% len.lambda
  best.beta.lambda[best.beta.lambda == 0] <- len.lambda
  best.beta <- beta[[best.beta.s]][,best.beta.lambda]
  
  validation.table <- data.frame(lambda=lambdas, s=ss, value=pv)
  results <- c(results, list(best.s=best.s, 
                             best.lambda=best.lambda,
                             best.pgs=best.pgs, 
                             best.beta=best.beta,
                             validation.table=validation.table, 
                             validation.type="pseudovalidation"))
  class(results) <- "validate.lassosum"
  if(plot) plot(results)
  return(results)
  
}

