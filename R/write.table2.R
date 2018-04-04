write.table2 <- function(...)  {
  options <- list(...)
  if(is.null(options$quote)) options$quote <- F
  if(is.null(options$row.names)) options$row.names <- F
  if(is.null(options$col.names)) options$col.names <- F
  do.call("write.table", options)
}