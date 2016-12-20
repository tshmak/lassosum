#' @title Computes polygenic scores as a genotype matrix multiplied by a matrix of weights
#'
#' @param bfile A plink bfile stem
#' @param weights The weights for the SNPs (\eqn{\beta})
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
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
                   chr=NULL) {

  stopifnot(is.numeric(weights))
  stopifnot(!any(is.na(weights)))
  if(is.vector(weights)) weights <- matrix(weights, ncol=1)
  stopifnot(is.matrix(weights))
  
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr)
  if(nrow(weights) != parsed$p) stop("Number of rows in (or vector length of) weights does not match number
				     of selected columns in bfile")
  # stopifnot(length(cor) == parsed$p)
  
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
  return(multiBed3(bfile, parsed$N, parsed$P, weights,
                   extract2[[1]], extract2[[2]], 
                   keepbytes, keepoffset))
  #' @return A matrix of Polygenic Scores
  
}
