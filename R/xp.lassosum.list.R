xp.lassosum.list <- function(xp.plink.linear, 
                             trace=0,
                             keep.ref=NULL, ...) {
  #' For a list of xp.plink.linear objects
  #' @keywords internal

  ss.list <- xp.plink.linear 
  best.pgs.t2 <- rep(NA, ss.list[[1]]$n)
  
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
  
  # Repeatedly call xp.lassosum
  l <- list()
  for(i in 1:length(ss.list)) {
    if(trace > 0) cat("Processing list item", i, "of", length(ss.list), "\n")
    if(is.null(keep.ref) & i > 1) {
      keep.ref <- attr(l[[1]], "keep.ref")
    } 
    l[[i]] <- xp.lassosum(ss.list[[i]], trace=trace-1,
                          keep.ref=keep.ref, 
                          list.of.lassosum.only=TRUE, ...)
    gc()
  }
  
  # Organize by folds 
  xp.l <- l[[1]]
  nfolds <- length(xp.l)
  ll <- l.folds <- list()
  for(f in 1:nfolds) {
    l.folds[[f]] <- lapply(l, function(x) x[[f]])
    ll[[f]] <- do.call("merge.lassosum.pipeline", l.folds[[f]])
  }
  
  return(ll)

}

  
  
