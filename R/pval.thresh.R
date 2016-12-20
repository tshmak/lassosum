pval.thresh <- function(pvals, p.thresholds, beta, bfile, 
                        keep=NULL, remove=NULL, extract=NULL, exclude=NULL, 
                        chr=NULL, 
                        mem.limit=4*10^9, chunks=NULL, cluster=NULL, trace=0) {
  #' Fast way to do p-value thresholding (without looping over the thresholds)
  #' 
  beta <- as.vector(beta)
  stopifnot(is.vector(pvals) & is.vector(beta))
  stopifnot(length(pvals) == length(beta))
  Pvals <- sort(unique(c(0,p.thresholds,1)))
  cut <- cut(pvals, Pvals, include.lowest = TRUE)
  stopifnot(!any(is.na(cut)))
  nlevels <- nlevels(cut)
  
  parsed <- parseselect(bfile, extract, exclude, keep, remove, chr)
  groups <- group.blocks(parseblocks(1:parsed$p), parsed, mem.limit, chunks, cluster)
  nchunks <- length(unique(groups$chunks.blocks))
  if(trace > 0) cat("Genotype divided into", nchunks, "chunks\n")
  
  dummy <- Diagonal(x=beta) %*% Matrix(model.matrix(~ cut - 1),sparse = T)
  
  if(is.null(parsed$extract)) select <- rep(TRUE, parsed$P) else 
    select <- parsed$extract
  
  if(is.null(cluster)) {
    l <- lapply(unique(groups$chunks), function(i) {
      if(trace > 0) cat("Processing chunk", i, "\n")
      select[select] <- groups$chunks == i
      mat <- readbfile(bfile, keep=parsed$keep, extract=select, fillmissing = T)
      return(mat %*% dummy[groups$chunks == i,])
    })
  } else {
    l <- parLapplyLB(cluster, unique(groups$chunks), function(i) {
      if(trace > 0) cat("Processing chunk", i, "\n")
      select[select] <- groups$chunks == i
      mat <- readbfile(bfile, keep=parsed$keep, extract=select, fillmissing = T)
      return(mat %*% dummy[groups$chunks == i,])
    })
  }

  sum <- function(...) {
    s <- colSums(do.call("rbind", lapply(list(...), as.vector)))
    return(matrix(s, nrow(list(...)[[1]])))
  }
  
  mat <- do.call(sum, l)
  result <- t(apply(t(mat), 2, cumsum))
  return(result)
  
}