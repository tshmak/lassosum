#' @title Function to read a text file 
#' @keywords internal
read.table2 <- function(file, header=F, data.table=F, check.names=TRUE, ...) {
  return(data.table::fread(file, header=header, data.table=data.table, 
                           check.names=check.names, ...))
}
