#' @title Elastic net using summary statistics
#' @description Coordinate descent algorithm to solve: 
#' 0.5 x'X'Xx - x'b + lambda1 ||x||_1 + 0.5 lambda2 ||x||_2^2
#' Function to get elastic net solutions given X, a reference panel, and
#' b, regression coefficients
#' @keywords internal
elnetR <- function(lambda1, lambda2=0, X, b, thr=1e-4,
                     trace=0, maxiter=10000, 
                   blocks=NULL, 
                   x=NULL) {
  stopifnot(length(b) == ncol(X))
  diag <- colSums(X^2)

  if(length(lambda2) > 1) {
    nlambda2 <- length(lambda2)
    for(i in 1:nlambda2) {
      result <- elnetR(lambda1, lambda2[i], X, b, thr,
                         trace, maxiter, x)
      result <- list(fit=result, lambda2=lambda2[i])
      if(i == 1) Result <- rep(result, nlambda2) else
        Result[i] <- result

    }
    return(Result)
  }

  order <- order(lambda1, decreasing = T)
  lambda1a <- lambda1[order]
  conv <- lambda1a * NA
  len <- length(b)
  beta <- matrix(NA, len, length(lambda1))
  pred <- matrix(NA, nrow(X), length(lambda1))
  loss <- rep(NA, length(lambda1))
  fbeta <- loss

  if(is.null(x)) x <- b * 0.0 else {
    stopifnot(length(x) == len)
    x <- x + 0.0 # Making sure R creates a copy...
  }

  if(is.null(blocks)) {
    Blocks <- list(startvec=0, endvec=len - 1)
  } else {
    Blocks <- parseblocks(blocks)
    stopifnot(max(Blocks$endvec)==len - 1)
  }
  
  X <- as.matrix(X)
  yhat <- as.vector(X %*% x)

  for(i in 1:length(lambda1a)) {
    if(trace > 0) cat("lambda1: ", lambda1a[i], "\n")
    conv[i] <- repelnet(lambda1a[i], lambda2, diag, X, b,thr,x,yhat, trace-1,maxiter,
                        Blocks$startvec, Blocks$endvec)
    if(conv[i] != 1) stop("Not converging...")

    beta[,i] <- x
    pred[,i] <- yhat
    loss[i] <- sum(yhat^2) - 2* sum(b * x)
    fbeta[i] <- loss[i] + 2* sum(abs(x))*lambda1a[i] + sum(x^2)*lambda2
  }


  conv[order] <- conv
  beta[,order] <- beta
  pred[,order] <- pred
  loss[order] <- loss
  fbeta[order] <- fbeta

  return(list(lambda1=lambda1, lambda2=lambda2, beta=beta, conv=conv, pred=pred, loss=loss, fbeta=fbeta))

}
