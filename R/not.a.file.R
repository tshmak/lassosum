not.a.file <- function(scalar_character) {
  #' A function to protect a single SNP id from being interpreted as a file
  attr(scalar_character, "not.a.file") <- TRUE
  return(scalar_character)
}