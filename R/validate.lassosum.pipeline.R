validate.lassosum.pipeline <- function(ls.pipeline, test.bfile=NULL, 
                              keep=NULL, remove=NULL, 
                              pheno=NULL, covar=NULL, 
                              validate.function=function(x, y) 
                                cor(x, y, use="complete"),
                              trace=1, 
                              destandardize=F, plot=T, 
                              exclude.ambiguous=T, 
                              cluster=NULL, ...) {
  
  #' @title Function to validate output from lassosum.pipeline with external phenotype
  #' @param ls.pipeline A lassosum.pipeline object
  #' @param test.bfile The (\href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK bfile} for the test dataset 
  #' @param keep Participants to keep (see \code{\link{lassosum}} for more details)
  #' @param remove Participants to remove
  #' @param pheno A vector of phenotype
  #' @param covar A matrix of covariates
  #' @param validate.function Function with which to perform validation
  #' @param trace Controls amount of output
  #' @param destandardize Should coefficients from \code{\link{lassosum}} be 
  #' destandardized using test dataset standard deviations before being returned?
  #' @param plot Should the validation plot be plotted? 
  #' @param exclude.ambiguous Should ambiguous SNPs (C/G, A/T) be excluded? 
  #' @param cluster A \code{cluster} object from the \code{parallel} package for parallel computing
  #' @param ... parameters to pass to \code{\link{sd.bfile}}
  #' @details Chooses the best \code{lambda} and \code{s} by validating 
  #' polygenic score against an external phenotype in the testing dataset. 
  #' If \code{pheno} is not specified, then the sixth column in the testing 
  #' dataset \href{https://www.cog-genomics.org/plink2/formats#fam}{.fam}\code{.fam} file is used. 
  #' @export
  stopifnot(class(ls.pipeline) == "lassosum.pipeline")
  
  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  
  if(!is.null(keep) || !is.null(remove)) if(is.null(test.bfile)) 
    stop("Please specify test.bfile if you specify keep or remove")
  
  redo <- T
  if(is.null(test.bfile)) {
    test.bfile <- ls.pipeline$test.bfile
    keep <- ls.pipeline$keep.test
    remove <- NULL
    redo <- F
  }
  
  if(destandardize) {
    if(ls.pipeline$destandardized) stop("beta in ls.pipeline already destandardized.")
    sd <- sd.bfile(test.bfile, extract=ls.pipeline$test.extract, 
                   keep=keep, remove=remove, ...)
    sd[sd <= 0] <- Inf # Do not want infinite beta's!
    ls.pipeline$beta <- lapply(ls.pipeline$beta, 
                               function(x) as.matrix(Matrix::Diagonal(x=1/sd) %*% x))
    redo <- T
  }
  
  if(redo) {
    ### Input Validation ### 
    extensions <- c(".bed", ".bim", ".fam")
    for(i in 1:length(extensions)) {
      if(!file.exists(paste0(test.bfile, extensions[i]))) {
        stop(paste0("File ", test.bfile, extensions[i], " not found."))
      }
    }
    ### Input Validation (end) ### 
    
    if(trace) cat("Coordinating lassosum output with test data...\n")
    
    bim <- fread(paste0(test.bfile, ".bim"))
    bim$V1 <- as.character(sub("^chr", "", bim$V1, ignore.case = T))
    
    m <- matchpos(ls.pipeline$sumstats, bim, auto.detect.ref = F, 
                       ref.chr = "V1", ref.pos="V4", ref.alt="V5", ref.ref="V6", 
                       rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                       silent=T)
    
    beta <- lapply(ls.pipeline$beta, function(x) 
      as.matrix(Matrix::Diagonal(x=m$rev) %*% x[m$order, ]))
    
    if(trace) cat("Calculating PGS...\n")
    
    pgs <- lapply(beta, function(x) pgs(bfile=test.bfile, weights = x, 
                                        extract=m$ref.extract, 
                                        keep=keep, remove=remove, 
                                        cluster=cluster))
    names(pgs) <- as.character(ls.pipeline$s)
    results <- c(results, list(pgs=pgs))

  } else {
    if(is.null(ls.pipeline$pgs)) {
      if(trace) cat("Calculating PGS...\n")
      pgs <- lapply(ls.pipeline$beta, function(x) pgs(bfile=test.bfile, 
                                          weights = x, 
                                          keep=keep, 
                                          cluster=cluster))
      names(pgs) <- as.character(ls.pipeline$s)
      results <- c(results, list(pgs=pgs))
    } else {
      results <- c(results, list(pgs=ls.pipeline$pgs))
    }
    beta <- ls.pipeline$beta
  } 

  ### Pheno ### 
  parsed.test <- parseselect(test.bfile, keep=keep, remove=remove)
  
  if(!is.null(pheno)) stopifnot(length(pheno) == parsed.test$n) else {
    fam <- fread(paste0(test.bfile, ".fam"))
    if(is.null(parsed.test$keep)) pheno <- fam$V6 else 
      pheno <- fam$V6[parsed.test$keep]
  }
  if(sd(pheno, na.rm = TRUE) == 0) stop("There's no variation in phenotype")
  
  ### Validate ###
  lambdas <- rep(ls.pipeline$lambda, length(ls.pipeline$s))
  ss <- rep(ls.pipeline$s, rep(length(ls.pipeline$lambda), length(ls.pipeline$s)))
  PGS <- do.call("cbind", results$pgs)

  ### covar ### 
  if(!is.null(covar)) {
    if(is.vector(covar)) covar <- matrix(covar, ncol=1)
    stopifnot(nrow(covar) == parsed.test$n) 
    for(i in ncol(PGS)) {
      PGS[,i] <- residuals(lm(PGS[,i] ~ ., data=covar))
    }
  }
  
  ### Validate (cont) ###
  suppressWarnings(cors <- as.vector(
    apply(PGS, MARGIN = 2, FUN=validate.function, pheno)))
  if(is.function(validate.function)) {
    funcname <- deparse(substitute(validate.function))
  } else if(is.character(validate.function)) {
    funcname <- validate.function
  } else {
    stop("What is validate.function? I can't figure out.")
  }

  cors[is.na(cors)] <- -Inf
  best <- which(cors == max(cors))[1]
  best.s <- ss[best]
  best.lambda <- lambdas[best]
  best.pgs <- PGS[,best]
  len.lambda <- length(ls.pipeline$lambda)
  best.beta.s <- ceiling(best / len.lambda)
  best.beta.lambda <- best %% len.lambda
  best.beta.lambda[best.beta.lambda == 0] <- len.lambda
  best.beta <- beta[[best.beta.s]][,best.beta.lambda]
  
  validation.table <- data.frame(lambda=lambdas, s=ss, value=cors)
  results <- c(results, list(best.s=best.s, 
                             best.lambda=best.lambda,
                             best.pgs=best.pgs, 
                             best.beta=best.beta, 
                             validation.table=validation.table, 
                             validation.type=funcname))
  class(results) <- "validate.lassosum"
  if(plot) plot(results)
  return(results)

}
