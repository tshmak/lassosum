#' @title Computes polygenic scores as a genotype matrix multiplied by a matrix of weights
#' @description Note that this is a low-level command. For applying estimated
#' betas to new bfiles, refer to https://github.com/tshmak/lassosum#apply-validated-betas-to-new-data 

pgs <- function(bfile, weights, ...) {
  #' @rdname pgs
  #' @export 
  return(.pgs(weights=weights, bfile=bfile, ...))
}
# This is to enable S3 parsing by the second argument
.pgs <- function(...) {
  #' @rdname pgs
  #' @export
  UseMethod(".pgs")
}
