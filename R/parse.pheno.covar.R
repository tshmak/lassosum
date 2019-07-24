parse.pheno.covar <- function(pheno, covar, parsed, trace=0) {
  #' @keywords internal
  fam <- parsed[['fam']]
  keep <- parsed$keep
  # keep <- NULL
  pheno.df <- NULL
  
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
  if(!is.null(pheno) && is.character(pheno) && length(pheno) == 1) {
    if(file.exists(pheno)) pheno <- read.table2(pheno, header=T) else 
      stop(paste("Cannot find", pheno))
  }
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
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    pheno.df <- pheno
    colnames(pheno.df)[3] <- "pheno"
    rownames(pheno) <- paste(pheno$FID, pheno$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(pheno))
    Pheno <- as.data.frame(pheno)[,3]
    names(Pheno) <- rownames(pheno)
  } else {
    if(!is.null(pheno)) {
      stopifnot(length(pheno) == parsed$n)
    } else {
      fam <- read.table2(parsed$famfile)
      if(is.null(parsed$keep)) pheno <- fam$V6 else 
        pheno <- fam$V6[parsed$keep]
    }
  }
  
  #### covar ####
  user.covar <- FALSE
  if(!is.null(covar) && is.character(covar) && length(covar) == 1) {
    if(file.exists(covar)) covar <- read.table2(covar, header=T) else 
      stop(paste("Cannot find", covar))
  }
  if(is.data.frame(covar) & all(colnames(covar)[1:2] == c("FID", "IID"))) {
    user.covar <- TRUE
    covar <- as.data.frame(covar)
    colnames <- colnames(covar) 
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    rownames(covar) <- paste(covar$FID, covar$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(covar))
    Covar <- covar[,-(1:2), drop=FALSE]
  } else {
    if(!is.null(covar)) {
      if(is.vector(covar)) covar <- matrix(covar, ncol=1)
      if(is.matrix(covar)) covar <- as.data.frame(covar)
      Covar <- covar
    } 
  }
  
  #### updates ####
  parsed$keep <- update.keep(parsed$keep, keep)
  if(!is.null(parsed$keep)) parsed$n <- sum(parsed$keep)
  if(is.data.frame(pheno)) {
    if(!is.null(parsed$keep)) {
      names <- rownames(fam)[parsed$keep]
    } else {
      names <- rownames(fam)
    }
    pheno <- Pheno[names] 
    if(trace) {
      message(length(pheno), " out of ", length(Pheno), " samples kept in pheno.")
      # message(paste("Note that the order of best.pgs is the order given in the .fam file", 
      #               " rather than the pheno data.frame. Use v$best.pgs[v$order] to get", 
      #               " the pgs in the order of the phenotype."))
    }
    Order <- 1:length(pheno)
    names(Order) <- names
    pheno.df$order <- Order[names(Pheno)]
  } 

  if(user.covar) {
    if(!is.null(parsed$keep)) covar <- Covar[rownames(fam)[parsed$keep],,drop=F] else 
      covar <- Covar[rownames(fam),,drop=F]
    if(trace) message(nrow(covar), " out of ", nrow(Covar), " samples kept in covar.")
  } 
  
  if(length(pheno) == 0) {
    stop("No phenotype left. Perhaps the FID/IID do not match?")
  } else if(length(pheno) != parsed$n) {
    stop("The length of pheno does not match the number of samples.")
  }
  if(!is.null(covar) && nrow(covar) != parsed$n) {
    stop(paste("The dimension of the covar matrix does not match the number of samples used.", 
               "If your covariate is a data.frame with FID and IID, make sure they have headers."))
  }
  # if(sd(pheno, na.rm = TRUE) == 0) stop("There's no variation in phenotype")
  parsed$fam <- fam

  return(list(pheno=pheno, covar=covar, parsed=parsed, table=pheno.df))
  
}



