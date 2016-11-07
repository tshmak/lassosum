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
#' @param samples to remove
#' @param chr a vector of chromosomes
#' @note Missing genotypes are interpreted as having the homozygous A2 alleles in the 
#' PLINK files (same as the \code{--fill-missing-a2} option in PLINK). 
#' @export
sd.bfile <- function(bfile, extract=NULL, exclude=NULL, 
                        keep=NULL, remove=NULL, chr=NULL) {
  
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr)
  return(lassosum(cor = rep(0.0, parsed$p), bfile = bfile, lambda=numeric(0), shrink=1, 
                                 keep=parsed$keep, extract=parsed$extract)$sd)
}
