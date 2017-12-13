#' @keywords internal

xp.cv.list <- function(list.of.xp.pl, # List of plink.linear
                       xp.lassosum, ...) {

  l <- xp.lassosum
  ldblocks <- l$split(l, l$LDblocks)  # split takes either 
                                      # lassosum.pipeline or xp.lassosum
                                      # as input  
  test.bfile <- l$test.bfile
  ref.bfile <- l$ref.bfile

  res <- list()
  for(i in 1:length(list.of.xp.pl)) {
    l$LDblocks <- ldblocks[[i]]
    
    if(length(test.bfile) == length(list.of.xp.pl)) {
      l$test.bfile <- test.bfile[i]
    } else {
      stopifnot(length(test.bfile) == 1)
    }
    
    if(length(ref.bfile) == length(list.of.xp.pl)) {
      l$ref.bfile <- ref.bfile[i]
    } else {
      stopifnot(length(ref.bfile) == 1)
    }
    
    res[[i]] <- xp.cv(list.of.xp.pl[[i]], l, return.lpipe=TRUE, ...)
  }  

  Res <- do.call("merge.lassosum.pipeline", res)
  return(Res)

  #' @return A lassosum.pipeline object
  
}
