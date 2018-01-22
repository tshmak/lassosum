check.class <- function(obj, class, list.of.class=TRUE) {
  #' Function to check whether obj belongs to class or
  #' is a list of class (if list.of.class==TRUE)
  #' @keywords internal
  
  err.message1 <- paste("Input must be a", class, "object")
  err.message2 <- paste(err.message1, "or a list of", class, "objects")

  is.list <- FALSE
  if(!is(obj, class)) {
    if(list.of.class) {
      if(is.list(obj)) {
        s <- sapply(obj, "class")
        if(!all(s==class)) {
          stop(err.message2)
        } else {
          is.list <- TRUE
          return(invisible(is.list))
        }
      } else {
        stop(err.message2)
      }
    } else {
      stop(err.message1)
    }
  }
  return(invisible(is.list))
}
