#' @title Function to read a table depending on whether data.table is installed
#' @keywords internal
read.table2 <- function(file, header=F) {
  if("data.table" %in% rownames(installed.packages())) {
    return(data.table::fread(file, data.table=F))
  } else {
    return(read.table(file, header = header, stringsAsFactors = F))
  }
}
