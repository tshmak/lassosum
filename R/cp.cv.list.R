#' @title cp.cv for a list of cp.plink.linear objects
#' @keywords internal

cp.cv.list <- function(list.of.cp.pl, # List of plink.linear
                       cp.lassosum, ...) {

  l <- cp.lassosum
  ldblocks <- l$split(l, l$LDblocks)  # split takes either 
                                      # lassosum.pipeline or cp.lassosum
                                      # as input  
  test.bfile <- l$test.bfile
  ref.bfile <- l$ref.bfile

  res <- list()
  for(i in 1:length(list.of.cp.pl)) {
    l$LDblocks <- ldblocks[[i]]
    
    if(length(test.bfile) == length(list.of.cp.pl)) {
      l$test.bfile <- test.bfile[i]
    } else {
      stopifnot(length(test.bfile) == 1)
    }
    
    if(length(ref.bfile) == length(list.of.cp.pl)) {
      l$ref.bfile <- ref.bfile[i]
    } else {
      stopifnot(length(ref.bfile) == 1)
    }
    
    res[[i]] <- cp.cv(list.of.cp.pl[[i]], l, return.lpipe=TRUE, ...)
  }  

  Res <- do.call("merge.lassosum.pipeline", res)
  return(Res)

  #' @return A lassosum.pipeline object
  
}
