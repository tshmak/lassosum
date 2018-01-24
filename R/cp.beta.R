#' @title Get the "overall" beta after cross-prediction
#' @param cp.plink.linear An object returned from cp.plink.linear()
#' @param cp.lassosum An object returned from cp.lassosum()
#' @param method Either the Method1 or Method2 Polygenic Score
#' @param save Save the pseudoinverse (in armadillo binary format)
#' @param load Load the pseudoinverse
#' @keywords internal

cp.beta <- function(cp.plink.linear, cp.lassosum, 
                    method=c("Method1", "Method2"), load="", 
                    save="", 
                    force=FALSE) {
  
  if(!is(cp.plink.linear, "cp.plink.linear")) 
    stop("Is this a cp.plink.linear object?")
  if(!is(cp.lassosum, "cp.lassosum")) 
    stop("Is this a cp.lassosum object?")
  
  ss <- cp.plink.linear
  l <- cp.lassosum 
  method <- match.arg(method)
  
  parsed <- parseselect(ss$bfile, extract=ss$extract, keep=ss$keep)
  
  if(parsed$n > parsed$p) {
    stop("cp.beta() doesn't work with n > p.")
  }
  
  if(load == "") {
    comp <- as.double(parsed$n)^2 * parsed$p
    mins <- round(comp / (2 * 1e10), 1)
    message(paste("This will take around", mins, "minutes (on a slow computer)."))
    if(mins > 60 & !force) {
      stop("This is expected to take very very long. Set force=TRUE if you really want to carry on.")
    } 
  }
  
  if(is.null(parsed$extract)) {
    extract2 <- list(integer(0), integer(0))
  } else {
    # print(parsed$extract)
    extract2 <- selectregion(!parsed$extract)
    extract2[[1]] <- extract2[[1]] - 1
  }
  
  if(is.null(parsed$keep)) {
    keepbytes <- integer(0)
    keepoffset <- integer(0)
  } else {
    pos <- which(parsed$keep) - 1
    keepbytes <- floor(pos/4)
    keepoffset <- pos %% 4 * 2
  }
  
  bedfile <- paste0(ss$bfile, ".bed")
  
  if(method == "Method1") {
    pred <- l$best.pgs
  } else if(method == "Method2") {
    pred <- l$best.pgs.m2
  }
  

  for(i in 1:length(ss$cor)) {
    na <- is.na(ss$cor[[i]])
    ss$cor[[i]][na] <- 0
  }
  
  meanbeta <- as.vector(rowMeans(as.data.frame(ss$cor)))
  obeta <- overallbeta(bedfile, parsed$N, parsed$P, 
              col_skip_pos=extract2[[1]], col_skip=extract2[[2]],
              keepbytes=keepbytes, keepoffset=keepoffset, 
              pred=pred, meanbeta=meanbeta, 
              save=save, load=load)
  return(as.vector(obeta))

  #' @return A vector of betas
  
}
