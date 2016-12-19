merge.lassosum <- function(...) {
  
  #' @title Merge lassosum results 
  #' @description e.g. when calculated over different blocks/chromosomes
  #' @method merge lassosum
  #' @export
  #' 
  ll <- list(...)
  stopifnot(all(sapply(ll, "class") == "lassosum"))
  stopifnot(all(sapply(ll, function(x) all(x$lambda == ll[[1]]$lambda))))
  shrink <- sapply(ll, function(x) x$shrink)
  stopifnot(all(shrink == shrink[1]))
  
  Cumsum <- function(...) {
    mat <- do.call("rbind", list(...))
    if(ncol(mat) > 0) return(as.vector(colSums(mat))) else
      return(numeric(0))
  }
  results <- list()
  results$lambda <- ll[[1]]$lambda
  results$beta <- do.call("rbind", lapply(ll, function(x) x$beta))
  results$conv <- do.call("pmin", lapply(ll, function(x) x$conv))
  pred <- do.call("Cumsum", lapply(ll, function(x) as.vector(x$pred)))
  results$pred <- matrix(pred, ncol=length(results$lambda), nrow=nrow(ll[[1]]$pred))
  results$loss <- do.call("Cumsum", lapply(ll, function(x) x$loss))
  results$fbeta <- do.call("Cumsum", lapply(ll, function(x) x$fbeta))
  results$sd <- do.call("c", lapply(ll, function(x) x$sd))
  results$shrink <- ll[[1]]$shrink
  results$nparams <- do.call("Cumsum", lapply(ll, function(x) x$nparams))
  class(results) <- "lassosum"
  return(results)
}