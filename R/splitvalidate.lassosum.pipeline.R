splitvalidate.lassosum.pipeline <- function(ls.pipeline, test.bfile=NULL, 
                                       keep=NULL, remove=NULL, 
                                       pheno=NULL, covar=NULL, 
                                       trace=1, split=NULL, 
                                       ...) {
  
  #' @title Function to perform split-validation using output from lassosum.pipeline with external phenotype
  #' @param ls.pipeline A lassosum.pipeline object
  #' @param test.bfile The (\href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK bfile} for the test dataset 
  #' @param keep Participants to keep (see \code{\link{lassosum}} for more details)
  #' @param remove Participants to remove
  #' @param pheno A vector of phenotype or a \code{data.frame} with 3 columns, the first 2 columns being headed "FID" and "IID"
  #' @param covar A matrix of covariates or a \code{data.frame} with 3 or more columns, the first 2 columns being headed "FID" and "IID"
  #' @param trace Controls amount of output
  #' @param ... parameters to pass to \code{\link{validate.lassosum.pipeline}}
  #' @details Performs split-validation. Randomly split the test data into half for validation 
  #' and half for prediction. Standardize the best cross-predicted pgs and stack together. 
  #' @rdname splitvalidate
  #' @export
  stopifnot(class(ls.pipeline) == "lassosum.pipeline")
  
  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  
  # if(!is.null(keep) || !is.null(remove)) if(is.null(test.bfile)) 
  #   stop("Please specify test.bfile if you specify keep or remove")
  
  redo <- T
  if(is.null(test.bfile)) {
    test.bfile <- ls.pipeline$test.bfile
    redo <- F
  }
  
  if(is.null(keep) && test.bfile == ls.pipeline$test.bfile) 
    keep <- ls.pipeline$keep.test
  
  ### Pheno & covar ### 
  parsed.test <- parseselect(test.bfile, keep=keep, remove=remove)
  phcovar <- parse.pheno.covar(pheno=pheno, covar=covar, parsed=parsed.test, 
                               trace=trace)
  parsed.test <- phcovar$parsed
  pheno <- phcovar$pheno
  covar <- phcovar$covar

  ### Split ###
  if(is.null(split)) {
    split <- sample(1:parsed.test$n %% 2 + 1)
  } else {
    stopifnot(length(split) == parsed.test$n)
    stopifnot(all(sort(unique(split)) == 1:max(split)))
  }

  ### Split-validation ###
  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  best.s <- best.lambda <- numeric(0)
  best.pgs <- pheno * NA
  best.beta <- numeric(0)
  validation.table <- data.frame()
  for(s in 1:max(split)) {
    if(is.null(parsed$keep)) {
      keep <- split == s
      pheno <- pheno[keep]
      if(!is.null(covar)) covar <- covar[keep,]
    }
    v <- validate(ls.pipeline, keep=keep, pheno=pheno, covar=covar, 
                  test.bfile=test.bfile, trace=trace, ...)
    best.s <- c(best.s, v$best.s)
    best.lambda <- c(best.lambda, v$best.lambda)
    best.pgs[keep] <- scale(v$best.pgs)
    best.beta <- cbind(best.beta, v$best.beta)
    validation.table <- rbind(validation.table, v$validation.table)
  }
  
  results <- c(results, list(split=split,
                             best.s=best.s, 
                             best.lambda=best.lambda,
                             best.pgs=best.pgs, 
                             best.beta=best.beta, 
                             validation.table=validation.table, 
                             validation.type=v$validation.type))
  class(results) <- "validate.lassosum"
  return(results)
  
}
