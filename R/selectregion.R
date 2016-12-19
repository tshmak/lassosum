#' @title Internal function to parse extract
#'
#' @param keep A boolean vector of which position to keep
#' @keywords internal
#' 
selectregion <- function(keep) {

  pos.keep <- which(keep)
  end <- pos.keep[-1]
  start <- pos.keep[-length(pos.keep)]
  diff.keep <- end - start
  skip <- diff.keep > 1
  starts <- c(start[1], end[skip])
  rle <- rle(keep)
  # print(rle)
  lengths <- rle$lengths[rle$values]
  return(list(starts, lengths))

}
