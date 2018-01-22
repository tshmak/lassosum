#' @title Perform cross-validation based on results from cross-prediction
#' @param cp.plink.linear An object returned from cp.plink.linear()
#' @param cp.lassosum An object returned from cp.lassosum()
#' @param return.lpipe For internal use
#' @export

cp.cv <- function(cp.plink.linear, cp.lassosum, 
                  return.lpipe=FALSE, ...) {
  
  ss <- cp.plink.linear
  l <- cp.lassosum 
  
  is.list <- check.class(ss, "cp.plink.linear", list.of.class=TRUE)

  if(is.list) {
    ls <- cp.cv.list(cp.plink.linear, cp.lassosum, ...)
  } else {
    if(class(cp.lassosum) != "cp.lassosum") 
      stop("Is this a cp.lassosum object?")
    
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
