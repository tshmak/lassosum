#' @title Subset a lassosum.pipeline object by lambda and s
#'
#' @param lassosum.pipeline An object returned by lassosum.pipeline
#' @param s Value(s) of s to restrict to 
#' @param lambda Value(s) of lambda to restrict to 
#' @details This function is usually used to reapply a validated pgs to a new data.set. 
#' See example below. 
#' 
#' @rdname subset.lassosum.pipeline
#' @export
subset.lassosum.pipeline <- function(lassosum.pipeline, s=NULL, lambda=NULL) {
  
  err <- function(param, value) {
    stop(paste("There is no", param, "equalling", value, "in lassosum.pipeline"))
  }
  lp <- lassosum.pipeline
  
  if(!is.null(s)) {
    w <- which(lp$s %in% s)
    if(length(w) == 0) err("s", s)
    lp$s <- lp$s[w]
    lp$beta <- lp$beta[w]
    lp$pgs <- lp$pgs[w]
  }
  
  if(!is.null(lambda)) {
    w <- which(lp$lambda %in% lambda)
    if(length(w) == 0) err("lambda", lambda)
    lp$lambda <- lp$lambda[w]
    for(i in 1:length(lp$s)) {
      lp$beta[[i]] <- lp$beta[[i]][,w, drop=F]
      lp$pgs[[i]] <- lp$pgs[[i]][,w, drop=F]
    }
  }
  
  #' @return A lassosum.pipeline object
  class(lp) <- "lassosum.pipeline"
  return(lp)
  
  #' @examples 
  #' \dontrun{
  #'  ### Run lassosum using standard pipeline ### 
  #'  lp <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position, 
  #'                           A1=ss$A1, A2=ss$A2,
  #'                           ref.bfile=ref.bfile, test.bfile=test.bfile, 
  #'                           LDblocks = ld)
  #'  v <- validate(lp)
  #'  lp2 <- subset(lp, s=v$best.s, lambda=v$best.lambda)
  #'  v2 <- validate(lp2)
  #' }
}
