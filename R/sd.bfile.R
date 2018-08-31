#' @title Obtain the SNP-wise standard deviations from the PLINK bfile
#' @details \code{keep}, \code{remove} could take one of three 
#' formats: (1) A logical vector indicating which indivduals to keep/remove, 
#' (2) A \code{data.frame} with two columns giving the FID and IID of the indivdiuals
#' to keep/remove (matching those in the .fam file), or (3) a character scalar giving the text file with the FID/IID. 
#' Likewise \code{extract}, \code{exclude} can also take one of the three formats,
#' except with the role of the FID/IID data.frame replaced with a character vector of 
#' SNP ids (matching those in the .bim file). 

#' @param bfile plink file stem
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
#' @note Missing genotypes are interpreted as having the homozygous A2 alleles in the 
#' PLINK files (same as the \code{--fill-missing-a2} option in PLINK). 
#' @keywords internal
sd.bfile <- function(bfile, keep=NULL, remove=NULL, extract=NULL, exclude=NULL, 
                        chr=NULL, trace=0, ...) {
  
  if(trace > 0) cat("Calculating SD...\n")
  if(length(bfile) > 1) {
    l <- splitvec.from.bfile(bfile)
    
    if(!is.null(extract)) {
      stopifnot(length(extract) == length(l$split_P))
      extract <- lapply(1:length(bfile), function(i) extract[l$split_P==i])
    }
    if(!is.null(exclude)) {
      stopifnot(length(exclude) == length(l$split_P))
      exclude <- lapply(1:length(bfile), function(i) exclude[l$split_P==i])
    }
    sd <- l$split_p * NA
    for(i in 1:length(bfile)) {
      sd[l$split_p == i] <- sd.bfile(bfile[i], keep=keep, remove=remove, 
                                     extract=extract[[i]], exclude=exclude[[i]], 
                                     chr=chr,trace=trace-1, ...)
    }
    return(sd)
  }
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr)
  return(lassosum(cor = rep(0.0, parsed$p), bfile = bfile, lambda=numeric(0), 
                  shrink=1, keep=parsed$keep, extract=parsed$extract, 
                  blocks=1:parsed$p, ...)$sd)
}
