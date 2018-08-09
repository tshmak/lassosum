parse.pheno.covar <- function(pheno, covar, parsed, trace=0) {
  #' @keywords internal

  fam <- NULL
  keep <- parsed$keep
  
  update.keep <- function(old, new) {
    if(all(new)) {
      return(old)
    } else {
      if(is.null(old)) return(new) else {
        if(is.null(new)) return(old) else 
          return(old & new)
      }
    }
  }
  #### pheno ####
  if(is.data.frame(pheno)) {
    if(ncol(pheno) != 3) {
      stop(paste("A pheno data.frame must have 3 columns exactly",
                 "with the first 2 with headers 'FID' and 'IID'"))
    }
    colnames <- colnames(pheno) 
    if(!all(colnames[1:2] == c("FID", "IID"))) {
      stop(paste("The first two columns of the pheno", 
                 "data.frame must have headers 'FID' and 'IID'"))
    }
    if(is.null(fam)) fam <- fread(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    rownames(pheno) <- paste(pheno$FID, pheno$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(pheno))
    Pheno <- pheno[,3]
  } else {
    if(!is.null(pheno)) {
      stopifnot(length(pheno) == parsed.test$n)
    } else {
      fam <- fread(parsed$famfile)
      if(is.null(parsed$keep)) pheno <- fam$V6 else 
        pheno <- fam$V6[parsed$keep]
    }
  }
  
  #### covar ####
  if(is.data.frame(covar)) {
    colnames <- colnames(covar) 
    if(!all(colnames[1:2] == c("FID", "IID"))) {
      stop(paste("The first two columns of the covar", 
                 "data.frame must have headers 'FID' and 'IID'"))
    }
    if(is.null(fam)) fam <- fread(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    rownames(covar) <- paste(covar$FID, covar$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(covar))
    Covar <- covar[,-(1:2), drop=FALSE]
  } else {
    if(!is.null(covar)) {
      if(is.vector(covar)) covar <- matrix(covar, ncol=1)
      if(is.matrix(covar)) covar <- as.data.frame(covar)
    } 
  }
  
  #### updates ####
  parsed$keep <- update.keep(parsed$keep, keep)
  if(!is.null(parsed$keep)) parsed$n <- sum(parsed$keep)
  if(is.data.frame(pheno)) {
    if(!is.null(parsed$keep)) pheno <- Pheno[rownames(fam)[parsed$keep]] else 
      pheno <- Pheno[rownames(fam)]
    if(trace) message(length(pheno), " out of ", length(Pheno), " samples kept in pheno.")
  } 

  if(is.data.frame(covar)) {
    if(!is.null(parsed$keep)) covar <- Covar[rownames(fam)[parsed$keep],] else 
      covar <- Covar[rownames(fam),]
    if(trace) message(nrow(covar), " out of ", nrow(Covar), " samples kept in covar.")
  } 
  
  stopifnot(length(pheno) == parsed$n)
  if(sd(pheno, na.rm = TRUE) == 0) stop("There's no variation in phenotype")

  return(list(pheno=pheno, covar=covar, parsed=parsed))
  
}



