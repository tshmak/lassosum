plot.validate.lassosum <- function(obj, legend.x="bottomright", log="x") {
  #' @title Plot function for \code{validate.lassosum} objects
  #' @export
  t <- obj$validation.table
  if(!any(is.finite(t$value))) return(invisible(NULL)) # If nothing to plot
  plot(unlist(t$lambda), 
       unlist(t$value), 
       type="n", log=log,
       xlab="lambda", ylab=obj$validation.type)
  ss <- unique(t$s)
  for(i in 1:length(ss)) {
    subs <- subset(t, s == ss[i])
    points(subs$lambda, subs$value, col=i, type="o")
  }
  legend(x=legend.x, col=1:length(ss), pch=1, 
         legend=paste0("s=", ss))

}
