#' @title Computes polygenic scores as a genotype matrix multiplied by a matrix of weights
#'
#' @param weights The weights for the SNPs (\eqn{\beta})
#' @param bfile A plink bfile stem
#' @param extract SNPs to extract (see \code{\link{parseselect}})
#' @param exclude SNPs to exclude (see \code{\link{parseselect}})
#' @param keep samples to keep (see \code{\link{parseselect}})
#' @param remove samples to remove (see \code{\link{parseselect}})
#' @param chr a vector of chromosomes
#' @param cluster A \code{cluster} object from the \code{parallel} package. 
#' For parallel processing. 
#' @param trace Level of output
#' @param sparse Assumes sparse weights matrix
#' @note \itemize{
#' \item Missing genotypes are interpreted as having the homozygous A2 alleles in the 
#' PLINK files (same as the \code{--fill-missing-a2} option in PLINK). 
#' \item The number of rows in \code{weights} should be the same as the number of
#' SNPs in the bfile after extract/exclude/chr.
#' }
#' @details A function to calculate \eqn{X\beta} where \eqn{X} is the genotype matrix
#' in the plink bfile. 
#' 
#' @rdname pgs
#' @export
.pgs.default <- function(weights, bfile, keep=NULL, extract=NULL, exclude=NULL, remove=NULL, 
                   chr=NULL, cluster=NULL, trace=0, sparse=TRUE) {

  if(length(bfile) > 1) {
    return(pgs.vec(bfile=bfile, weights=weights, keep=keep, remove=remove,
                   extract=extract, exclude=exclude, chr=chr, 
                   cluster=cluster, trace=trace, sparse=sparse))
  }

  stopifnot(is.numeric(weights))
  stopifnot(!any(is.na(weights)))
  if(is.vector(weights)) weights <- matrix(weights, ncol=1)
  stopifnot(is.matrix(weights))
  
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr, order.important=TRUE)
  if(nrow(weights) != parsed$p) stop("Number of rows in (or vector length of) weights does not match number of selected columns in bfile")
  # stopifnot(length(cor) == parsed$p)
  
  if(!is.null(cluster)) {
    nclusters <- length(cluster)
    if(nclusters > 1) {
      CW <- cumsum(rowSums(weights != 0))
      split <- ceiling(CW / max(CW) * nclusters)
      split[split == 0] <- 1
      # split <- ceiling(seq(1/parsed$p, nclusters, length=parsed$p))
      t <- table(split)
      compute.size <- as.double(min(t)) * parsed$N * ncol(weights)
      if(compute.size < 1e8 || sum(t > 0) < nclusters) {
        # Too many clusters
        if(sum(t > 0) < nclusters) {
          return(pgs(bfile, weights, keep=parsed$keep, extract=parsed$extract, 
                     trace=trace, sparse=sparse))
        } else {
          f <- 1e8 / compute.size
          recommended <- min(ceiling(nclusters / f), nclusters - 1)
          return(pgs(bfile, weights, keep=parsed$keep, extract=parsed$extract, 
                     cluster=cluster[1:recommended], trace=trace, sparse=sparse))
        }
      }
      Bfile <- bfile # Define this within the function so that it is copied
                      # to the child processes
      l <- parallel::parLapply(cluster, 1:nclusters, function(i) {
        toextract <- if(!is.null(parsed$extract)) parsed$extract else 
          rep(TRUE, parsed$P)
        touse <- split == i
        toextract[toextract] <- touse
        
        return(pgs(Bfile, weights[touse, ], keep=parsed$keep, extract=toextract, 
                   trace=trace, sparse=sparse))
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

  if(!sparse) {
    return(multiBed3(bfile, parsed$N, parsed$P, weights,
                     extract2[[1]], extract2[[2]], 
                     keepbytes, keepoffset, trace=trace))
  } else {
    ss <- Matrix::summary(Matrix::Matrix(t(weights), sparse = TRUE))
    nonzeros <- as.integer(table(factor(ss$j, levels=1:nrow(weights))))
    colpos <- ss$i - 1
    
    return(multiBed3sp(bfile, parsed$N, parsed$P, 
                       beta=ss$x, nonzeros=nonzeros, colpos=colpos, ncol=ncol(weights), 
                       extract2[[1]], extract2[[2]], 
                       keepbytes, keepoffset, trace=trace))
  }
  #' @return A matrix of Polygenic Scores
  
}
