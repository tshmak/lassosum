plot.validate.lassosum <- function(obj) {
  #' @title Plot function for \code{validate.lassosum} objects
  #' @export
  t <- obj$validation.table
  plot(unlist(t$lambda), 
       unlist(t$value), 
       type="n", 
       xlab="lambda", ylab=obj$validation.type)
  ss <- unique(t$s)
  for(i in 1:length(ss)) {
    subs <- subset(t, s == ss[i])
    points(subs$lambda, subs$value, col=i, type="o")
  }
  legend(x="topright", col=1:length(ss), pch=1, 
         legend=paste0("s=", ss))

}
