cp.lassosum.list <- function(cp.plink.linear, 
                             trace=0,
                             keep.ref=NULL, ...) {
  #' For a list of cp.plink.linear objects
  #' @keywords internal

  ss.list <- cp.plink.linear 
  best.pgs.m2 <- rep(NA, ss.list[[1]]$n)
  
  # Checks 
  n <- ss.list[[1]]$n
  stopifnot(all(sapply(ss.list, function(x) 
    identical(x$n, n))))
  
  fold <- ss.list[[1]]$fold
  stopifnot(all(sapply(ss.list, function(x) 
    identical(x$fold, fold))))
  
  pheno.by.fold <- ss.list[[1]]$pheno.by.fold
  stopifnot(all(sapply(ss.list, function(x) 
    identical(x$pheno.by.fold, pheno.by.fold))))
  
  keep <- ss.list[[1]]$keep
  stopifnot(all(sapply(ss.list, function(x) 
    identical(x$keep, keep))))
  
  # Repeatedly call cp.lassosum
  l <- list()
  for(i in 1:length(ss.list)) {
    if(trace > 0) cat("Processing list item", i, "of", length(ss.list), "\n")
    if(is.null(keep.ref) & i > 1) {
      keep.ref <- l[[1]]$keep.ref # Keep the same keep.ref from before
    } 
    l[[i]] <- cp.lassosum(ss.list[[i]], trace=trace-1,
                          keep.ref=keep.ref, 
                          list.of.lpipe.output=TRUE, ...)
    gc()
  }
  
  # Organize by folds 

  ll <- organize.by.fold(l)
  
  return(ll)

}

  
  
