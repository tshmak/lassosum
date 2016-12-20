logical.vector <- function(positions, size) {
  ### Function to create a logical vector based on position and
  ### size, i.e. total length of vector
  ### It's basically the reverse of findpos
  vec <- rep(F, size)
  vec[positions] <- T
  return(vec)
}