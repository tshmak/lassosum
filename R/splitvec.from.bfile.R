splitvec.from.bfile <- function(bfile) {
  #' Function to create split vectors from the info in vectors of bfiles
  #' @keywords internal
  pvec <- attr(bfile, "p")
  Pvec <- attr(bfile, "P")
  if(is.null(pvec) || is.null(Pvec)) {
    stop("Since v0.4.3 bfile should have the attributes 'p' and 'P' if it is a vector")
  } else {
    stopifnot(length(pvec) == length(bfile))
    stopifnot(length(Pvec) == length(bfile))
  }
  
  return(list(split_p=rep(1:length(bfile), attr(bfile, "p")), 
              split_P=rep(1:length(bfile), attr(bfile, "P"))))
  
}

