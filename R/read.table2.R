#' @title Function to read a text file depending on whether data.table is installed
#' @keywords internal
read.table2 <- function(file, header=F) {
  if("data.table" %in% rownames(installed.packages())) {
    return(data.table::fread(file, data.table=F, header=header))
  } else {
    return(read.table(file, header = header, stringsAsFactors = F))
  }
}
