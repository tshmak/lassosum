xp.plink.linear <- function(bfile, nfolds=5, fold=NULL, 
                            pheno=NULL, covar=NULL, 
                            keep=NULL, remove=NULL, 
                            extract=NULL, exclude=NULL, 
                            chr=NULL, 
                            force=FALSE, 
                            fast=TRUE, 
                            ...) {
  #' @title Generate summary statistics for cross-prediction
  #' @param nfolds Number of folds
  #' @param fold A vector of fold ID
  #' @param pheno A vector of phenotype
  #' @param covar A vector/matrix of covariates
  #' @param keep,remove,extract,exclude,chr see \code{parseselect}
  #' @param force Force
  #' @param fast see details 
  #' @details If fast == TRUE, the summary statistics are calculated for 
  #'          each fold, and then combined together by averaging. For example,
  #'          assume there are 5 folds, then for fold 1, the summary stats
  #'          for fold 2, 3, 4, and 5 are averaged. If fast == FALSE, 
  #'          we actually run plink using the entire dataset minus the fold 1, 
  #'          and repeat this process for folds 2, 3, 4, and 5.
  #' @export
  
  parsed <- parseselect(bfile=bfile, keep=keep, remove=remove, 
                        extract=extract, exclude=exclude, chr=chr, 
                        export=TRUE)
  
  #### folds ####
  if(!is.null(fold)) {
    if(length(fold) != parsed$n) stop("Length of fold vector doesn't match number of samples")
    folds <- unique(fold)
    fold <- as.integer(factor(fold, levels=folds))
    nfolds <- length(folds)
  } else {
    fold <- sample((1:parsed$n) %% nfolds + 1)
  }
  
  t <- table(fold) 
  if(min(t) < 50 & !force) 
    stop(paste("The minimum fold has less than 50 observations.",
               "Perhaps reduce the number of folds,", 
               "or specify force=T to run anyway."))
  
  if(!is.null(parsed$keep)) {
    touse <- which(parsed$keep)  
  } else {
    touse <- 1:parsed$n
  }
  
  #### pheno ####
  if(is.null(pheno)) {
    if(is.null(parsed$fam)) parsed$fam <- read.table2(parsed$famfile)
    pheno <- parsed$fam[,6]
    if(!is.null(parsed$keep)) pheno <- pheno[parsed$keep]
  } else {
    if(length(pheno) != parsed$n) stop("Length of phenotype vector does not match number of samples")
  }

  u <- unique(pheno)
  u.nonNA <- u[!is.na(u)]
  if(all(u %in% c(-9, 0, 1, 2))) {
    binary <- TRUE
    pheno[pheno==-9] <- NA
    pheno <- pheno - 0.5 # To make it continuous in PLINK
  } else binary <- FALSE

  pheno.by.fold <- lapply(1:nfolds, function(i) pheno[fold == i])
  n.phenos <- sapply(pheno.by.fold, function(x) length(unique(x)))
  if(binary) {
    tt <- sapply(pheno.by.fold, function(x) min(table(x)))
    if(min(tt) < 20 & !force) stop("Some of the folds have less than 20 cases/controls. Perhaps reduce the number of folds.")
  }
  
  if(any(n.phenos < 2)) stop("Some of the folds have no variation in the phenotype!")
  
  #### Cross-prediction ####
  result <- list()
  for(i in 1:nfolds) {
    if(fast) training <- logical.vector(touse[fold == i], parsed$N) else 
      training <- logical.vector(touse[fold != i], parsed$N)
    result[[i]] <- plink.linear(bfile=bfile, 
                                pheno=pheno[training],
                                keep = training, 
                                extract= parsed$extract, ...)
  }
  
  cor <- lapply(result, function(x) x$BETA)
  chr <- lapply(result, function(x) x$CHR)
  pos <- lapply(result, function(x) x$BP)
  A1 <- lapply(result, function(x) x$A1)
  snp <- lapply(result, function(x) x$SNP)
  nmiss <- as.data.frame(lapply(result, function(x) x$NMISS))
  
  # Checks 
  stopifnot(all(sapply(chr, function(x) identical(x, chr[[1]]))))
  stopifnot(all(sapply(pos, function(x) identical(x, pos[[1]]))))
  stopifnot(all(sapply(A1, function(x) identical(x, A1[[1]]))))
  stopifnot(all(sapply(snp, function(x) identical(x, snp[[1]]))))
  
  chr <- chr[[1]]; pos <- pos[[1]]; A1 <- A1[[1]]; snp <- snp[[1]]
  
  if(fast) {
    Cor <- as.data.frame(cor)
    for(i in 1:nfolds) {
      cor[[i]] <- rowMeans(Cor[,-i, drop=FALSE])
    }
  }

  toreturn <- list(cor=cor, chr=chr, pos=pos, A1=A1, snp=snp, 
                fold=fold, 
                pheno.by.fold=pheno.by.fold, 
                keep=parsed$keep, extract=parsed$extract, 
                n=parsed$n, p=parsed$p, nonmiss=nmiss, 
                bfile=bfile)
  class(toreturn) <- "xp.plink.linear"
  return(toreturn)

}
