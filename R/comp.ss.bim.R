#' @title Compare summary statistics and PLINK .bim data.frames
#' 
#' @param ss.df A data.frame of summary statistics
#' @param bim.df a data.frame from A PLINK .bim file
#' @param exclude.ambiguous should ambiguous SNPs (i.e. A/T or G/C) be excluded
#' @details Both \code{ss.df} and \code{bim.df} must be \code{data.frame}s with exactly 
#' three columns in this order: (1) SNP ids, (2) reference allele, (3) alternative allele
#' 
#' @keywords internal
#' #@export
comp.ss.bim <- function(ss.df, bim.df, exclude.ambiguous=F) {
  
  stopifnot(is.data.frame(ss.df))
  stopifnot(is.data.frame(bim.df))
  
  stopifnot(ncol(ss.df) == 3)
  stopifnot(ncol(bim.df) == 3)
  
  colnames(ss.df) <- c("rsid", "ref.ss", "alt.ss")
  colnames(bim.df) <- c("rsid", "ref", "alt")
  
  ss.df <- cbind(ss.df, index.ss=1:nrow(ss.df))
  bim.df <- cbind(bim.df, index=1:nrow(bim.df))
  
  #merged <- merge(bim.df, ss.df, all.x=F, all.y = F)
  merged <- merge(bim.df, ss.df, all=F)
  merged <- merged[order(merged$index),]
  
  merged$ok <- merged$ref == merged$ref.ss & merged$alt == merged$alt.ss
  merged$antiok <- merged$ref == merged$alt.ss & merged$alt == merged$ref.ss
  merged$eitherok <- merged$ok | merged$antiok
  
  # exclude.ambiguous <- F
  if(exclude.ambiguous) {
    merged$ambiguous <- with(merged, 
                             (ref == "A" & alt == "T") | 
                               (ref == "T" & alt == "A") | 
                               (ref == "G" & alt == "C") | 
                               (ref == "C" & alt == "G") ) 
  } else merged$ambiguous <- F
  
  ss.order <- merged$index.ss[merged$eitherok & !merged$ambiguous]
  bim.extract <- rep(F, nrow(bim.df))
  bim.extract[merged$index[merged$eitherok & !merged$ambiguous]] <- T
  rev <- ifelse(merged$antiok[merged$eitherok & !merged$ambiguous], -1,1)
  
  return(list(ss.order=ss.order, bim.extract=bim.extract, rev=rev))
  #' @return A list with 3 vectors 
  #' \item{ss.order}{the order the SNPs in the summary statistics must be rearranged in}
  #' \item{bim.extract}{the SNPs to be extracted from the PLINK bfile, i.e. that 
  #' are in the summary statistics data.frame}
  #' \item{rev}{A vector of 1 and -1 to facilitate the flipping SNPs effects for SNPs whose
  #' reference and alternative alleles are reversely coded between the two data.frames}
  
}
