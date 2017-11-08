xp.lassosum <- function(xp.plink.linear, 
                        LDblocks=xp.plink.linear$chr, 
                        pseudovalidation=FALSE, 
                        ref.bfile=NULL, 
                        destandardize=TRUE, 
                        max.ref.bfile.n=5000, 
                        details=FALSE, 
                        keep.ref=NULL, 
                        exclude.ambiguous=FALSE,
                        validate.function="cor", 
                        plot=TRUE, 
                        trace=1, 
                        cluster=NULL, 
                        ...) {
  #' @title lassosum with cross-prediction
  #' @description We assume correlations are pre-calculated for the various 
  #' folds, using PLINK or otherwise (cos PLINK is a lot faster).
  #' @param xp.plink.linear An object from xp.plink.linear()
  #' @param LDblocks LD blocks. See \code{\link{lassosum.pipeline}}
  #' @param pseudovalidation Should pseudovalidation be performed on each fold?
  #' @param ref.bfile bfile of reference panel (if different from that used in \code{xp.plink.linear})
  #' @param destandardize Should coefficients be destandardized
  #' @param max.ref.bfile.n Maximum number of samples to use in ref.bfile (to improve speed)
  #' @param details Should a detailed output be given?
  #' @param keep.ref Participants to keep in the reference panel (see \code{parseselect} for expected input) 
  #' @param exclude.ambiguous Should ambiguous SNPs be excluded? 
  #' @param validate.function Function for validating polygenic score
  #' @param plot Whether a validation plot should be drawn
  #' @param ... Other parameters to pass to lassosum.pipeline()

  # ref.bfile <- NULL; max.ref.bfile.n <- 100; details <- TRUE; destandardize=TRUE
  ss <- xp.plink.linear # get a shorter name
  if(is.null(keep.ref)) {
    if(is.null(ref.bfile)) refsamplesize <- ss$n else {
      refsamplesize <- nrow.bfile(ref.bfile)
    }
    if(refsamplesize > max.ref.bfile.n) {
      touse <- sample(refsamplesize, max.ref.bfile.n)
      if(is.null(ss$keep)) {
        keep.ref <- logical.vector(touse, refsamplesize)
      } else {
        keep.ref <- ss$keep
        keep.ref[keep.ref] <- logical.vector(touse, refsamplesize)
      }
    } else {
      keep.ref <- ss$keep
    }
  } 
  nfolds <- length(ss$cor)
  l <- list()
  pv <- list()
  for(i in 1:nfolds) {
    # i <- 1
    if(trace > 0) cat("Processing fold", i, "of", nfolds, "\n")
    
    na <- is.na(ss$cor[[i]])
    # if(mean(na) > 0.1) {
    #   stop("There seems to be a lot of NA's in the correlations...")
    # }
    ss$cor[[i]][na] <- 0
    l[[i]] <- lassosum.pipeline(cor=ss$cor[[i]], 
                                chr=ss$chr, 
                                pos=ss$pos,
                                snp=ss$snp,
                                A1=ss$A1, 
                                LDblocks=LDblocks,
                                exclude.ambiguous=exclude.ambiguous, 
                                destandardize=destandardize, 
                                ref.bfile=ref.bfile, 
                                test.bfile=ss$bfile, 
                                keep.ref=keep.ref, 
                                cluster=cluster,
                                # Using the same reference 
                                # panel for all folds.
                                # This is neater and shouldn't
                                # affect performance much
                                # I think... 
                                keep.test=ss$fold == i, 
                                trace=trace-1, 
                                ...)  # no need to include extract
                                      # as they will not have been
                                      # included in ss$cor anyway. 
    if(pseudovalidation) {
      pv[[i]] <- pseudovalidate.lassosum.pipeline(l[[i]], 
               trace=trace - 1, plot=plot, 
               cluster=cluster)
    }
  }
  
  pgs <- l[[1]]$pgs
  for(s in 1:length(pgs)) {
    pgs[[s]] <- matrix(NA, nrow=ss$n, ncol=ncol(pgs[[1]]))
  }
  best.pgs.pv <- pheno <- fold <- rep(NA, ss$n)
  for(i in 1:nfolds) {
    fold[ss$fold == i] <- i
    pheno[ss$fold == i] <- ss$pheno.by.fold[[i]]
    for(s in 1:length(pgs)) {
      pgs[[s]][ss$fold == i, ] <- l[[i]]$pgs[[s]]
    }
    if(pseudovalidation) {
      best.pgs.pv[ss$fold == i] <- pv[[i]]$best.pgs
    }
  }
  
  # Validation 
  ll <- l[[1]]
  ll$pgs <- pgs
  ll$keep.test <- ss$keep
  v <- validate.lassosum.pipeline(ll, pheno=pheno, trace=trace-1, plot=plot,
                                  validate.function = validate.function, 
                                  cluster=cluster)

  tab <- data.frame(pheno=pheno, fold=fold, best.pgs=v$best.pgs)
  if(pseudovalidation) tab <- cbind(tab, best.pgs.pv)
  attr(tab, "keep.ref") <- keep.ref
  
  if(details) {
    attr(tab, "lassosum") <- l
    attr(tab, "validate") <- v
    if(pseudovalidation) attr(tab, "pseudovalidate") <- pv
  } 
  
  return(tab)
  
}
