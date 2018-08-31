pgs <- function(bfile, weights, ...) {
  #' @rdname pgs
  #' @export 
  return(.pgs(weights=weights, bfile=bfile, ...))
}
#' This is to enable S3 parsing by the second argument
.pgs <- function(...) {
  #' @rdname pgs
  #' @export
  UseMethod(".pgs")
}
