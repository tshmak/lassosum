#' @title Computes polygenic scores as a genotype matrix multiplied by a matrix of weights
#'
#' @param bfile A plink bfile stem
#' @param weights The weights for the SNPs (\eqn{\beta})
#' @param extract SNPs to extract (see \code{\link{parseselect}})
#' @param exclude SNPs to exclude (see \code{\link{parseselect}})
#' @param keep samples to keep (see \code{\link{parseselect}})
#' @param remove samples to remove (see \code{\link{parseselect}})
#' @param chr a vector of chromosomes
#' @param cluster A \code{cluster} object from the \code{parallel} package. 
#' For parallel processing. 
#' @param pbin p-value bins (for \code{\link{pval.thresh}}) 
#' @note \itemize{
#' \item Missing genotypes are interpreted as having the homozygous A2 alleles in the 
#' PLINK files (same as the \code{--fill-missing-a2} option in PLINK). 
#' \item The number of rows in \code{weights} should be the same as the number of
#' SNPs in the bfile after extract/exclude/chr.
#' }
#' @details A function to calculate \eqn{X\beta} where \eqn{X} is the genotype matrix
#' in the plink bfile. 
#' 

#' @export
pgs <- function(bfile, weights, keep=NULL, extract=NULL, exclude=NULL, remove=NULL, 
                   chr=NULL, cluster=NULL, pbin=NULL) {

  if(length(bfile) > 1) {
    call <- match.call()
    return(do.call("pgs.vec", as.list(call[-1])))
  }

  stopifnot(is.numeric(weights))
  stopifnot(!any(is.na(weights)))
  if(is.vector(weights)) weights <- matrix(weights, ncol=1)
  stopifnot(is.matrix(weights))
  
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr)
  if(nrow(weights) != parsed$p) stop("Number of rows in (or vector length of) weights does not match number of selected columns in bfile")
  # stopifnot(length(cor) == parsed$p)
  
  if(!is.null(cluster)) {
    nclusters <- length(cluster)
    if(nclusters > 1) {
      split <- ceiling(seq(1/parsed$p, nclusters, length=parsed$p))
      t <- table(split)
      compute.size <- as.double(min(t)) * parsed$N * ncol(weights)
      if(compute.size < 1e8) {
        # Too many clusters
        f <- 1e8 / compute.size
        recommended <- min(ceiling(nclusters / f), nclusters - 1)
        return(pgs(bfile, weights, keep=parsed$keep, extract=parsed$extract, 
                   cluster=cluster[1:recommended], pbin=pbin))
      }
      Bfile <- bfile # Define this within the function so that it is copied
                      # to the child processes
      l <- parallel::parLapply(cluster, 1:nclusters, function(i) {
        toextract <- if(!is.null(parsed$extract)) parsed$extract else 
          rep(TRUE, parsed$P)
        touse <- split == i
        toextract[toextract] <- touse
        
        pbin.touse <- pbin[touse]
        attr(pbin.touse, "nbin") <- attr(pbin, "nbin")
        return(pgs(Bfile, weights[touse, ], keep=parsed$keep, extract=toextract, 
                   pbin=pbin.touse))
      })
      result <- l[[1]]
      if(nclusters > 1) for(i in 2:nclusters) result <- result + l[[i]]
      return(result)
    }
  }
  
  if(is.null(parsed$extract)) {
    extract2 <- list(integer(0), integer(0))
  } else {
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
  
  bfile <- paste0(bfile, ".bed")
  
  if(is.null(pbin)) {
    return(multiBed3(bfile, parsed$N, parsed$P, weights,
                     extract2[[1]], extract2[[2]], 
                     keepbytes, keepoffset))
  } else {
    if(length(weights) != length(pbin)) {
      stop("Length of pbin doesn't match length of weights")
    }
    nbin <- attr(pbin, "nbin")
    if(is.null(nbin)) {
      stop("Perhaps call pval.thresh instead of directly calling pgs()?")
    }
    stopifnot(max(pbin) < nbin)
    return(multiBed4(bfile, parsed$N, parsed$P, 
                     weights, pbin, nbin, 
                     extract2[[1]], extract2[[2]], 
                     keepbytes, keepoffset))
    
  }
  #' @return A matrix of Polygenic Scores
  
}
