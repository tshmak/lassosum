merge.lassosum.pipeline <- function(...) {
  #' Merge a list of lassosum.pipeline objects
  #' @keywords internal
  
  l <- list(...)
  names(l) <- NULL # This will save a LOT of memory!!!
                   # Otherwise c() starts concantenating the names also
                   # and uses a LOT OF memory!!!

  lpipe <- l[[1]]
  if(length(l) > 1) {
    for(i in 2:length(l)) {
        # i <- 1; j <- 1
        
      # Checks
      checks <- c("lambda", "s", "keep.test", "keep.ref", 
                  "destandardized")
      for(m in 1:length(checks)) {
        if(!identical(lpipe[[checks[m]]],
                            l[[i]][[checks[m]]])) {
          stop(paste(m, "is not the same in the",
                     paste0(i,"th"), "element of the list."))
        }
      }

      # pgs
      for(k in 1:length(lpipe$beta)) {
        lpipe$pgs[[k]] <- lpipe$pgs[[k]] + l[[i]]$pgs[[k]]
      }
    }
    
    # vecs
    vecs <- c("test.extract", "also.in.refpanel", "sd", 
              "ref.bfile", "test.bfile", "LDblocks", "time")
    for(m in 1:length(vecs)) {
      lpipe[[vecs[m]]] <- do.call("c", lapply(l, function(x) x[[vecs[m]]]))
    }

    # beta
    for(k in 1:length(lpipe$beta)) {
      lpipe$beta[[k]] <- do.call("rbind", lapply(l, function(x) x$beta[[k]]))
    }
    
    # sumstats
    lpipe$sumstats <- do.call("rbind", lapply(l, function(x) l$sumstats))
    
    # Checks if same file name. 
    if(all(lpipe$ref.bfile == lpipe$ref.bfile[1])) {
      lpipe$ref.bfile <- lpipe$ref.bfile[1]
    }

    if(all(lpipe$test.bfile == lpipe$test.bfile[1])) {
      lpipe$test.bfile <- lpipe$test.bfile[1]
    }
  }
  
  split.vec <- sapply(l, function(x) nrow(x$beta[[1]]))
  lpipe$beta.split <- rep(1:length(l), split.vec)
  lpipe$split <- function(obj, vec) {
    # obj should be either lassosum.pipeline or xp.lassosum
    s <- obj$beta.split
    return(lapply(1:max(s), function(i) vec[s==i]))
  }
  class(lpipe) <- "lassosum.pipeline"
  return(lpipe)
}
