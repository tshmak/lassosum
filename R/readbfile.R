#' @title read a PLINK bfile file into a matrix
#' @param bfile PLINK bfile (as character, without the .bed extension)
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
#' @param fillmissing Whether to fill missing values with 0
#' @keywords internal

readbfile <- function(bfile, keep=NULL, extract=NULL, exclude=NULL, remove=NULL, 
                    chr=NULL, fillmissing=F) {
  
  
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr)
  
  if(is.null(parsed$extract)) {
    extract2 <- list(integer(0), integer(0))
  } else {
    # print(parsed$extract)
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
  
  bedfile <- paste0(bfile, ".bed")
  return(genotypeMatrix(bedfile, parsed$N, parsed$P, 
                        col_skip_pos=extract2[[1]], col_skip=extract2[[2]],
                        keepbytes=keepbytes, keepoffset=keepoffset, 
                        fillmissing=as.integer(fillmissing)))
  #' @return A genotype matrix of 0,1, and 2 (possibly with NaNs)
  
}
