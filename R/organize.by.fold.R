organize.by.fold <- function(list.of.list.of.lp) {
  #' lassosum.pipeline objects, orginally organized by fold then by chromosome, 
  #' are merged across chromosome
  #' @param list.of.list.of.lp A list of list of \code{lassosum.pipeline} objects
  #'                           from the output of \code{cp.lassosum}, when
  #'                           \code{list.of.lpipe.output=TRUE}
  #' @export
  
  l <- list.of.list.of.lp
  nfolds <- length(l[[1]])
  ll <- l.folds <- list()
  for(f in 1:nfolds) {
    # f <- 1
    l.folds <- lapply(l, function(x) x[[f]])
    ll[[f]] <- do.call("merge.lassosum.pipeline", l.folds)
    if(f == 1) ll <- rep(ll, nfolds) # memory pre-allocation 
  }
  return(ll)
}
