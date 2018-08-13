parseargs <- function(..., silent=F, include.others=TRUE) {
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
    Text <- readline()
    text <- sub('^"(.*)"$', '\\1', Text)
    if(text == Text) {
      text <- sub("^'(.*)'$", '\\1', Text)
    }
  }
  # text <- "abc.sh --abc 234 --def-3 --uio 23/asdf/asdf"
  sp <- strsplit(text, split="[[:space:]]--")[[1]]
  opts[[".OTHER."]] <- sp[1]
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
      val <- sp2[[i]][2]
      if(!silent) cat(paste0("--", Names[i], " ", val, "\n"))
      if(grepl("\\.file[01]?[tsc]?$", Names[i])) {
        header <- !grepl("0[tsc]?$", Names[i])
        sep <- ifelse(grepl("t$", Names[i]), "\t", ifelse(grepl("c$", Names[i]), ",", ""))
        val <- read.table(val, header=header, sep=sep)
        if(ncol(df) == 1) {
          val <- val[,1]
        }
        names[i] <- sub("\\.file[01]?[tsc]?$", "", names[i])
      } else {
        if(!is.na(as.logical(val))) {
          val <- as.logical(val)
        } else {
          suppressWarnings(if(!is.na(as.numeric(val))) val <- as.numeric(val))
        }
      }
      opts[[names[i]]] <- val
      if(length(sp2[[i]]) > 2) {
        others <- sp2[[i]][-(1:2)]
        opts[[paste0(".OTHER.", names[i])]] <- others
        if(i == length(sp2)) {
          opts[[".OTHER.last"]] <- others
        }
      }
    } else {
      opts[[names[i]]] <- TRUE
    }
  }
  if(!include.others) {
    exclude <- grepl("^\\.OTHER\\.", names(opts)) | grepl("^\\.SYSTEM", names(opts))
    opts <- opts[!exclude]
  }
  
  return(opts)
}
