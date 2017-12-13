#' @title pgs for a list of bfiles
#' 
#' @keywords internal
pgs.vec <- function(bfile.vec, weights.list, extract=NULL, exclude=NULL, 
                     ...) {
  
  #### Checks ####
  mess <- function(x) paste("If bfile is given as a list, so must", x)
  mess2 <- function(x) paste("Length of list of", x, 
                             "doesn't match length of list of bfiles")
  if(!is.list(weights.list)) stop(mess("weight"))
  if(length(weights.list) != length(bfile.list)) 
    stop(mess2("weights"))
  if(!is.null(extract)) {
    if(!is.list(extract)) stop(mess("extract"))
    if(length(extract) != length(bfile.list)) stop(mess2("extract"))
  }
  if(!is.null(exclude)) {
    if(!is.list(exclude)) stop(mess("exclude"))
    if(length(exclude) != length(bfile.list)) stop(mess2("exclude"))
  }
  
  #### Start ####
  for(i in 1:length(bfile.list)) {
    PGS <- pgs(bfile.list[i], weights.list[[i]], 
               extract=extract[[i]], exclude=exclude[[i]], ...)
    if(i == 1) res <- PGS else res <- res + PGS
  }
  return(res)
  
}
  
