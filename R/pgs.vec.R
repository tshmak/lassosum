#' @title pgs for a list of bfiles
#' 
#' @keywords internal
pgs.vec <- function(bfile, weights, extract=NULL, exclude=NULL, 
                     ...) {

  pvec <- attr(bfile, "p")
  Pvec <- attr(bfile, "P")
  if(is.null(pvec) || is.null(Pvec)) {
    stop("Since v0.4.3 bfile should have the attributes 'p' and 'P' if it is a vector")
  } else {
    stopifnot(length(pvec) == length(bfile))
    stopifnot(length(Pvec) == length(bfile))
  }
  
  split1 <- rep(1:length(bfile), attr(bfile, "p"))
  split2 <- rep(1:length(bfile), attr(bfile, "P"))
  
  if(is.vector(weights)) weights <- matrix(weights, ncol=1)
  weights <- lapply(1:length(bfile)) function(i) weights[split1==i,]

  if(!is.null(extract)) {
    stopifnot(length(extract) == length(Pvec))
    extract <- lapply(1:length(bfile)) function(i) extract[split2==i]
  }
  if(!is.null(exclude)) {
    stopifnot(length(exclude) == length(Pvec))
    exclude <- lapply(1:length(bfile)) function(i) exclude[split2==i]
  }

  #### Start ####
  for(i in 1:length(bfile)) {
    PGS <- pgs(bfile[i], weights.list[[i]], 
               extract=extract[[i]], exclude=exclude[[i]], ...)
    if(i == 1) res <- PGS else res <- res + PGS
  }
  return(res)
  
}
  
