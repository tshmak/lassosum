#' @title Perform cross-validation based on results from cross-prediction
#' @param xp.plink.linear An object returned from xp.plink.linear()
#' @param xp.lassosum An object returned from xp.lassosum()
#' @param return.lpipe For internal use
#' @export

xp.cv <- function(xp.plink.linear, xp.lassosum, 
                  return.lpipe=FALSE, ...) {
  
  ss <- xp.plink.linear
  l <- xp.lassosum 
  
  is.list <- check.class(ss, "xp.plink.linear", list.of.class=TRUE)

  if(is.list) {
    ls <- xp.cv.list(xp.plink.linear, xp.lassosum, ...)
  } else {
    if(class(xp.lassosum) != "xp.lassosum") 
      stop("Is this a xp.lassosum object?")
    
    #### Get ss ####
    chr <- ss$chr
    snp <- ss$snp
    pos <- ss$pos
    A1 <- ss$A1
    cor <- as.vector(rowMeans(as.data.frame(ss$cor)))
    cor[is.na(cor)] <- 0
    
    #### Cross-validation ####
    ls <- lassosum.pipeline(cor=cor, chr=chr, 
                            snp=snp, pos=pos, 
                            A1 = A1, nomatch=TRUE, 
                            LDblocks = l$LDblocks, 
                            lambda = l$best.lambda,
                            s=l$best.s, 
                            test.bfile = ss$bfile,
                            ref.bfile = l$ref.bfile, 
                            keep.ref=l$keep.ref, 
                            keep.test=l$keep.test, 
                            exclude.ambiguous = l$exclude.ambiguous, 
                            destandardize=l$destandardized, 
                            ...)
  }
  
  if(return.lpipe) return(ls)
  
  v <- validate.lassosum.pipeline(ls, pheno=l$pheno, plot=FALSE)
  
  return(v)

  #' @return A validate.lassosum.pipeline object
  
}
