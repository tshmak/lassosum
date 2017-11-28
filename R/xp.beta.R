#' @title Get the "overall" beta after cross-prediction
#' @param xp.plink.linear An object returned from xp.plink.linear()
#' @param xp.lassosum An object returned from xp.lassosum()
#' @param type Either the Type1 or Type2 Polygenic Score
#' @param save Save the pseudoinverse (in armadillo binary format)
#' @param load Load the pseudoinverse
#' @export

xp.beta <- function(xp.plink.linear, xp.lassosum, 
                    type=c("Type1", "Type2"), load="", 
                    save="", 
                    force=FALSE) {
  
  if(!is.list(xp.plink.linear) != "xp.plink.linear") 
    stop("Is this a xp.plink.linear object?")
  if(!is.data.frame(xp.lassosum)) 
    stop("Is this a xp.lassosum object?")
  
  ss <- xp.plink.linear
  l <- xp.lassosum 
  type <- match.arg(type)
  
  parsed <- parseselect(ss$bfile, extract=ss$extract, keep=ss$keep)
  
  if(parsed$n > parsed$p) {
    stop("xp.beta() doesn't work with n > p.")
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
  
  if(type == "Type1") {
    pred <- l$best.pgs
  } else if(type == "Type2") {
    pred <- l$best.pgs.t2
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
