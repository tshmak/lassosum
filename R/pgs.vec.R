#' @title pgs for a list of bfiles
#' 
#' @keywords internal
pgs.vec <- function(bfile, weights, extract=NULL, exclude=NULL, 
                    trace=0, ...) {


  l <- splitvec.from.bfile(bfile)
  
  if(is.vector(weights)) weights <- matrix(weights, ncol=1)
  stopifnot(nrow(weights) == length(l$split_p))
  weights <- lapply(1:length(bfile), function(i) weights[l$split_p==i,])

  if(!is.null(extract)) {
    stopifnot(length(extract) == length(l$split_P))
    extract <- lapply(1:length(bfile), function(i) extract[l$split_P==i])
  }
  if(!is.null(exclude)) {
    stopifnot(length(exclude) == length(l$split_P))
    exclude <- lapply(1:length(bfile), function(i) exclude[l$split_P==i])
  }

  #### Start ####
  for(i in 1:length(bfile)) {
    if(trace > 0) cat("Processing ", bfile[i], "\n")
    PGS <- pgs(bfile[i], weights[[i]], 
               extract=extract[[i]], exclude=exclude[[i]], 
               trace=trace-1, ...)
    if(i == 1) res <- PGS else res <- res + PGS
  }
  return(res)
  
}
  
