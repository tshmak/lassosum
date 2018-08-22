parseargs <- function(required=character(0), ..., silent=F, include.others=TRUE, 
                      test=NULL) {

  #' @title Function to parse arguments from command line
  #' @keywords internal 
  
  opts <- list(...)
  if(!interactive()) {
    text <- commandArgs(trailingOnly = FALSE)
    # str(text)
    start <- which(text == "--args")
    if(length(start) == 0) {
      opts[[".SYSTEM"]] <- text
      text <- ""
    } else {
      opts[[".SYSTEM"]] <- text[1:start]
      text <- paste(text[-(1:start)], collapse = " ")
    }
  } else {
    if(!is.null(test)) Text <- test else 
      Text <- readline(prompt = "Arguments? ")
    text <- sub('^"(.*)"$', '\\1', Text)
    if(text == Text) {
      text <- sub("^'(.*)'$", '\\1', Text)
    }
  }
  sp <- strsplit(text, split="(^|[[:space:]])--")[[1]]
  opts[[".OTHER.before"]] <- sp[1]
# print(text)
# print(sp)
  sp <- sp[-1]
  if(length(sp) == 0) return(opts)
  
  sp2 <- strsplit(sp, split="[[:space:]]")
  
  Names <- sapply(sp2, function(x) x[1])
  names <- as.character(Names)
  if(!silent) message("Options: ")
  for(i in 1:length(sp2)) {
    if(length(sp2[[i]]) > 1) {
      val <- paste(sp2[[i]][-1], collapse=" ")
    } else {
      val <- TRUE
    }
    if(!silent) cat(paste0("--", Names[i], " ", val, "\n"))
    opts[[names[i]]] <- val
    
    if(i == length(sp2)) {  # Options not preceded by --
      if(length(sp2[[i]]) > 1) opts[[".OTHER.after"]] <- paste(sp2[[i]][-1], collapse=" ")
    }
  }
  if(!include.others) {
    exclude <- grepl("^\\.OTHER\\.", names(opts)) | grepl("^\\.SYSTEM", names(opts))
    opts <- opts[!exclude]
  }

  for(i in 1:length(opts)) {
    suppressWarnings(if(!is.logical(opts[[i]]) && !is.na(as.numeric(opts[[i]]))) 
      opts[[i]] <- as.numeric(opts[[i]]))
  }

  if(length(required) > 0) {
    for(v in required) {
      if(length(opts) == 0 || !(v %in% names(opts))) stop(paste0("--", v, " is required."))
    }
  }
  
  return(opts)
}
