xp.lassosum.validate <- function(list.of.lpipe, # by fold
                                 xp.plink.linear, 
                                 Method2, pseudovalidation, 
                                 scale=TRUE, 
                                 plot, 
                                 validate.function,
                                 cluster, 
                                 trace, 
                                 details, 
                                 list.of.bfiles) {
  #' For a list of xp.plink.linear objects
  #' @keywords internal

  # A function for scaling
  Scale <- function(x, scale) {
    if(scale) {
      x <- scale(x)
      x[is.nan(x)] <- 0
      return(x) 
    } else return(x)
  }
  
  # Get shorthands  
  l <- list.of.lpipe
  ss <- xp.plink.linear
  
  # Check
  if(!is(ss, "xp.plink.linear")) {
    if(is.list(ss) && all(sapply(ss, class) == "xp.plink.linear")) {
      ss <- ss[[1]]
    } else {
      stop("Unknown type for list.of.lpipe.")
    }
  }

  # Initialize stacked PGS matrix
  pgs <- l[[1]]$pgs
  for(s in 1:length(pgs)) {
    pgs[[s]] <- matrix(NA, nrow=ss$n, ncol=ncol(pgs[[1]]))
  }
  
  # Initialize 
  m2.fold <- best.pgs.pv <- best.pgs.m2 <- pheno <- fold <- rep(NA, ss$n)
  M2 <- pv <- list()
  M2.beta <- numeric(0)
  
  ### Loop over folds
  nfolds <- length(ss$pheno.by.fold)
  for(i in 1:nfolds) {
    
    # pseudovalidation
    if(pseudovalidation) {
      pv[[i]] <- pseudovalidate.lassosum.pipeline(l[[i]], 
                                                  trace=trace, plot=plot, 
                                                  cluster=cluster)
      best.pgs.pv[ss$fold == i] <- pv[[i]]$best.pgs
    }
    
    # Method 2
    if(Method2) {
      split <- sample(as.logical(round(seq(0,1, length=sum(ss$fold == i)))))
      split <- rep(list(split), 2)
      split[[2]] <- !split[[2]]
      pheno2 <- rep(list(ss$pheno.by.fold[[i]]), 2)
      pheno2[[1]][split[[1]]] <- NA
      pheno2[[2]][split[[2]]] <- NA
      
      m2 <- list()
      m2.pgs <- rep(NA, length(ss$pheno.by.fold[[i]]))
      for(ii in 1:2) {
        m2[[ii]] <- validate.lassosum.pipeline(l[[i]], pheno=pheno2[[ii]], 
                                               trace=trace, plot=plot, 
                                               validate.function=validate.function, 
                                               cluster=cluster)
        m2.pgs[split[[ii]]] <- Scale(m2[[ii]]$best.pgs[split[[ii]]], scale)
        M2.beta <- cbind(M2.beta, m2[[ii]]$best.beta)
      }
      best.pgs.m2[ss$fold == i] <- m2.pgs
      m2.fold[ss$fold == i] <- 
        as.numeric(paste0(i, ".", as.integer(split[[2]]) + 1))
      M2[[i]] <- m2
    }
    
    # Stacking pgs
    # fold[ss$fold == i] <- i # Seems to be useless... 
    pheno[ss$fold == i] <- ss$pheno.by.fold[[i]]
    for(s in 1:length(pgs)) {
      pgs[[s]][ss$fold == i, ] <- Scale(l[[i]]$pgs[[s]], scale)
    }
    
  }
  
  # Validation 
  ll <- l[[1]]
  ll$pgs <- pgs
  ll$keep.test <- ss$keep
  v <- validate.lassosum.pipeline(ll, pheno=pheno, trace=trace-1, plot=plot,
                                  validate.function = validate.function, 
                                  cluster=cluster)
  
  # betas
  lambda.pos <- which(l[[1]]$lambda == v$best.lambda)
  for(i in 1:nfolds) {
    beta <- l[[i]]$beta[[as.character(v$best.s)]][,lambda.pos]
    if(i == 1) Beta <- beta else Beta <- cbind(Beta, beta)
  }

  # Generate output table
  tab <- list(pheno=pheno, fold=ss$fold, 
              best.pgs=v$best.pgs, 
              best.beta=Beta,
              best.lambda=v$best.lambda, 
              best.s=v$best.s)
  if(pseudovalidation) tab$best.pgs.pv <- best.pgs.pv
  if(Method2) {
    tab$m2.fold <- m2.fold
    tab$best.pgs.m2 <- best.pgs.m2
    tab$best.beta.m2 <- M2.beta
  }

  if(details) {
    attr(tab, "validate") <- v
    attr(v, "remarks") <- 
      paste("best.beta represents the best.beta in the first fold, ", 
            "selected using the best lambda and s from validation.")

    attr(tab, "lassosum") <- l
    
    if(pseudovalidation) attr(tab, "pseudovalidate") <- pv
    if(Method2) attr(tab, "Method2") <- M2
  } 

  return(tab)

}