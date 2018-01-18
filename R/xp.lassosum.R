xp.lassosum <- function(xp.plink.linear, 
                        LDblocks=xp.plink.linear$chr, 
                        pseudovalidation=FALSE,
                        Type2=FALSE, 
                        scale=TRUE, 
                        ref.bfile=NULL, 
                        destandardize=FALSE, 
                        max.ref.bfile.n=5000, 
                        details=FALSE, 
                        keep.ref=NULL, 
                        exclude.ambiguous=FALSE,
                        validate.function=function(x, y) 
                          cor(x,y, use="complete"),
                        plot=FALSE, 
                        trace=1, 
                        cluster=NULL, 
                        list.of.lpipe.output=FALSE, 
                        list.of.lpipe.input=NULL, 
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
  is.list <- check.class(ss, "xp.plink.linear", list.of.class=TRUE)
  
  if(is.null(list.of.lpipe.input)) {

    if(is.list) {
      call <- match.call()
      l <- do.call("xp.lassosum.list", as.list(call[-1]))
      ss <- ss[[1]]
    } else {
      if(is.null(keep.ref)) {
        if(is.null(ref.bfile)) {
          refsamplesize <- ss$n 
        } else {
          refsamplesize <- nrow.bfile(ref.bfile)
        }
        
        if(refsamplesize > max.ref.bfile.n) {
          touse <- sample(refsamplesize, max.ref.bfile.n)
        } else {
          touse <- 1:refsamplesize
        }
        
        if(is.null(ss$keep)) {
          keep.ref <- logical.vector(touse, refsamplesize)
        } else {
          keep.ref <- ss$keep
          keep.ref[keep.ref] <- logical.vector(touse, refsamplesize)
        }
      } else {
        stopifnot(is.logical(keep.ref))
        if(is.null(ref.bfile)) {
          stopifnot(length(keep.ref) == nrow.bfile(ss$bfile))
        } else {
          stopifnot(length(keep.ref) == nrow.bfile(ref.bfile))
        }
      }
      
      if(is.null(ss$keep)) {
        fold.test <- ss$fold
      } else {
        fold.test <- rep(0, ss$n)
        fold.test[ss$keep] <- ss$fold
      }
      
      nfolds <- length(ss$cor)
      l <- list()
      pv <- list()
      t2 <- list()
      best.pgs.t2 <- rep(NA, ss$n)
      test.bfile <- ss$bfile
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
                                    nomatch=TRUE, 
                                    LDblocks=LDblocks,
                                    exclude.ambiguous=exclude.ambiguous, 
                                    destandardize=destandardize, 
                                    ref.bfile=ref.bfile, 
                                    keep.ref=keep.ref, 
                                    # Using the same reference 
                                    # panel for all folds.
                                    # This is neater and shouldn't
                                    # affect performance much
                                    # I think... 
                                    test.bfile=test.bfile, 
                                    keep.test=fold.test == i, 
                                    cluster=cluster,
                                    trace=trace-1, 
                                    ...)  # no need to include extract
        # as they will not have been
        # included in ss$cor anyway. 
        if(i == 1) LDblocks <- l[[i]]$LDblocks
      }
      # For use by xp.lassosum.list() only 
      if(list.of.lpipe.output) {
        return(l)
      } 
    }
  } else {
    l <- list.of.lpipe.input # by fold
    check.class(l, "lassosum.pipeline")
  }
  
  # Is ss a list?
  if(is.list) ss <- ss[[1]]
  
  # Get results 
  result <- xp.lassosum.validate(l, ss, Type2=Type2, 
    pseudovalidation=pseudovalidation, scale=scale, plot=plot, 
    validate.function=validate.function, cluster=cluster, 
    trace=trace-1, details=details)
  result$split <- l[[1]]$split
  result$beta.split <- l[[1]]$beta.split
  result$ref.bfile <- l[[1]]$ref.bfile
  result$keep.ref <- l[[1]]$keep.ref
  result$test.bfile <- l[[1]]$test.bfile
  result$keep.test <- ss$keep  # l[[1]] is for fold 1 only
  result$fold <- ss$fold  # l[[1]] is for fold 1 only
  result$LDblocks <- l[[1]]$LDblocks
  result$destandardized <- l[[1]]$destandardized
  result$exclude.ambiguous <- l[[1]]$exclude.ambiguous
  result$time <- sum(unlist(lapply(l, function(x) x$time)))
  
  class(result) <- "xp.lassosum"

  return(result)

}
