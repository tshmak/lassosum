group.blocks <- function(Blocks, parseselect, mem.limit=2^32*8, chunks=NULL, 
                         cluster=NULL) {
  #' @title Group blocks into chunks so as not to exhaust memory or hit 
  #' matrix size ceiling
  #' 
  #' @param parseselect An object returned by parseselect()
  #' @keywords internal

  ncols <- Blocks$endvec - Blocks$startvec + 1
  max.size <- max(ncols) * parseselect$n * 8
  if(max.size > mem.limit) stop("Size of blocks too large. Suggest split into smaller blocks or increase mem.limit or use only a sample of individuals.")
  cum.size <- cumsum(ncols * parseselect$n * 8)
  required.memory <- sum(ncols) * parseselect$n * 8
  
  #### Determine number of chunks from cluster ####
  if(!is.null(cluster)) {
    stopifnot(inherits(cluster, "cluster"))
    if(is.null(chunks)) {
      nclusters <- length(cluster)
      mem.limit <- mem.limit / nclusters
      chunks <- max(c(nclusters, ceiling(required.memory / mem.limit)))
    }
  }
  
  #### Determine chunk split ####
  if(is.null(chunks)) {
    ok <- T; i <- 1
    groups <- rep(0, length(cum.size))
    while(any(groups == 0)) {
      groups[cum.size < mem.limit & groups==0] <- i
      cum.size <- cum.size - max(cum.size[groups == i])
      i <- i + 1
    }
    chunks <- rep(groups, ncols)
  } else {
    if(is.vector(chunks) && length(chunks) == parseselect$p) {
      table <- !(table(rep(1:length(ncols), ncols), as.vector(chunks)) == 0)
      if(any(rowSums(table) != 1)) {
        stop("The 'chunks' split doesn't agree with the 'blocks' split. Make sure 'blocks' are specified and that they are finer units than 'chunks'.")
      }
      groups <- as.vector(table %*% as.integer(colnames(table)))
    } else if(length(chunks) == 1 && is.vector(chunks)) {
      rough.split <- required.memory/chunks
      ok <- T; i <- 1
      groups <- rep(0, length(cum.size))
      while(any(groups == 0)) {
        newgroup <- cum.size <= rough.split & groups==0
        groups[newgroup] <- i
        current.index <- max(which(groups > 0), 0)
        add1 <- current.index + 1
        # add1 <- max(which(newgroup))+1
        if(add1 <= length(groups)) groups[add1] <- i
        cum.size <- cum.size - max(cum.size[groups == i])
        if(max(cum.size[groups == i]) * parseselect$n * 8 > mem.limit) {
          stop("Chunks too large. Either increase the number of chunks or mem.limit.")
        }
        i <- i + 1
      }
      chunks <- rep(groups, ncols)
    } else {
      stop("I don't know what to do with this 'chunks' input.")
    }
  }  
  
  if(is.null(parseselect$extract)) {
    extracts <- lapply(1:max(groups), function(i) chunks == i)
  } else {
    extracts <- lapply(1:max(groups), 
                       function(i) logical.vector(which(parseselect$extract)[chunks == i],
                                                  parseselect$P))
  }
  
  return(list(chunks.blocks=groups, chunks=chunks, extracts=extracts))
}