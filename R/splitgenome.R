splitgenome <- function(CHR, POS, ref.CHR, ref.breaks, details=T, right=TRUE) {
  #' Function to split a set of SNPs by their position using a reference
  #' Break points should be given for the reference by chromosome
  #' @keywords internal
  
  CHR <- as.character(CHR)
  ref.CHR <- as.character(ref.CHR)
  POS <- as.integer(POS)
  ref.breaks <- as.integer(ref.breaks)
  
  stopifnot(all(!is.na(POS)) && all(!is.na(ref.breaks)))
  stopifnot(all(POS >= 1) && all(ref.breaks >= 1))
  stopifnot(length(CHR) == length(POS))
  stopifnot(length(ref.CHR) == length(ref.breaks))
  
  
  chr <- (unique(CHR))
  chr.ref <- (unique(ref.CHR))
  included.chr <- chr %in% chr.ref
  # if(all(!included.chr)) stop("Cannot match the chromosomes. Make sure the notations are the same. e.g. 'chr1' vs 1 or chrX vs chr23.")
  if(!all(included.chr)) stop("Some chromosomes were not defined in the reference. Make sure the notations are the same. e.g. 'chr1' vs 1 or chrX vs chr23.")
  
  levels <- character(0)
  results <- character(length(POS))
  Details <- data.frame()
  for(C in chr.ref) {
    breaks <- sort(unique(ref.breaks[ref.CHR == C]))
    if(breaks[1] > 1) breaks <- c(1, breaks)
    if(breaks[length(breaks)] < Inf) breaks <- c(breaks, Inf)
    cut <- cut(POS[CHR == C],include.lowest = T,breaks = breaks, right=right)
    levels <- c(levels, paste0(C, "_", levels(cut)))
    cut <- paste0(C, "_", as.character(cut))
    results[CHR == C] <- cut
    if(details) {
      df <- data.frame(chr=C, start=breaks[-length(breaks)], end=breaks[-1])
      Details <- rbind(Details, df)
    }
  }
  results <- factor(results, levels = levels)
  
  if(details) {
    Details$counts <- as.integer(table(results))
    attr(results, "details") <- Details
  }
  
  return(results)
  
}