logical.vector <- function(positions, size) {
  
  #' @title Function to create a logical vector based on position and
  #' size, i.e. total length of vector
  #' @description It's basically the reverse of base::which()
  #'
  #' @param positions Positions to extract
  #' @param size Length of logical vector
  #' @export
  vec <- rep(F, size)
  vec[positions] <- T
  return(vec)
}