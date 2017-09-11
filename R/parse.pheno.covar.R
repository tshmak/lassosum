parse.pheno.covar <- function(input, parsed) {
  
  #' Parse the pheno and covariates input
  if(!is.null(input)) {
    if(is.character(input) && length(input) == 1) {
      ## A file is given
      toreturn <- input
    } else if(is.numeric(input) || is.data.frame(input)) {
      if(is.vector(input)) input <- matrix(input, ncol=1)
      if(is.matrix(input)) input <- as.data.frame(input)
      if(!all(c("FID", "IID") %in% colnames(input))) {
        if(is.null(parsed$fam)) parsed$fam <- read.table2(parsed$famfile) 
        input.data.frame <- parsed$fam[,1:2]
        stopifnot(nrow(input) == parsed$N || nrow(input) == parsed$n)
        if(nrow(input) == parsed$N) {
          input <- cbind(input.data.frame, input)
        } else {
          input <- cbind(input.data.frame[parsed$keep, ], input)
        }
      }
      write.table2(input,file=inputfile <- tempfile(pattern="input"))
      toreturn <- inputfile
    } else {
      stop("Don't know what to do yet...")
    }
    return(toreturn)
  }
  
} 

