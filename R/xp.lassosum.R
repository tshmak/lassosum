xp.lassosum <- function(xp.plink.linear, 
                        LDblocks=xp.plink.linear$chr, 
                        pseudovalidation=FALSE,
                        Type2=FALSE, 
                        ref.bfile=NULL, 
                        destandardize=TRUE, 
                        max.ref.bfile.n=5000, 
                        details=FALSE, 
                        keep.ref=NULL, 
                        exclude.ambiguous=FALSE,
                        validate.function=function(x, y) 
                          cor(x,y, use="complete"),
                        plot=FALSE, 
                        trace=1, 
                        cluster=NULL, 
                        list.of.lassosum.only=FALSE, 
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
  err.message <- "Input must be an xp.plink.linear object or a list of these objects."
  if(class(ss) != "xp.plink.linear") {
    if(is.list(ss)) {
      if(all(sapply(ss, class) == "xp.plink.linear")) {
        multiple.ss <- TRUE
        call <- match.call()
        l <- do.call("xp.lassosum.list", as.list(call[-1]))
        ss <- ss[[1]]
        # Need to replace test.bfile with a single bfile to trick validate.lassosum.pipeline
        # Otherwise will give an error
        real.test.bfiles <- l[[1]]$test.bfile
        for(i in 1:length(l)) {
          l[[i]]$test.bfile <- l[[i]]$test.bfile[1]
        }
      } else stop(err.message)
    } else stop(err.message)
  } else {
    multiple.ss <- FALSE
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
    t2 <- list()
    best.pgs.t2 <- rep(NA, ss$n)
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
      
    }
    # For use by xp.lassosum.list() only 
    if(list.of.lassosum.only) {
      return(l)
    } 
  }
  
  # Get results 
  
  result <- xp.lassosum.validate(l, ss, Type2=Type2, 
    pseudovalidation=pseudovalidation, plot=plot, 
    validate.function=validate.function, cluster=cluster, 
    trace=trace-1, details=details)

  ### Attributes see in xp.lassosum.validate() ###
  # if(details) {
  #   attr(tab, "lassosum") <- l
  #   attr(tab, "validate") <- v
  #   if(pseudovalidation) attr(tab, "pseudovalidate") <- pv
  #   if(Type2) attr(tab, "Type2") <- T2
  # } 
  
  attr(result, "keep.ref") <- keep.ref

  if(multiple.ss) {
    for(i in 1:length(l)) {
      l[[i]]$test.bfile <- real.test.bfiles
    }
  }
  
  return(result)
  
  
}
