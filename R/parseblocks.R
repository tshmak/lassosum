parseblocks <- function(vec) {

  #' @keywords internal
  if(is.factor(vec)) vec <- as.integer(vec)
  vec <- as.integer(factor(vec, levels=unique(vec)))
  blocks <- unique(vec)
  stopifnot(blocks == sort(blocks))
  stopifnot(min(blocks) == 1)
  stopifnot(max(blocks) == length(blocks))
  rle <- rle(vec)
  endvec <- cumsum(rle$lengths)
  startvec <- c(0, endvec[-length(endvec)])
  endvec <- endvec - 1 
  return(list("startvec"=startvec, "endvec"=endvec))
  
}