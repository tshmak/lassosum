merge.lassosum.pipeline <- function(...) {
  #' Merge a list of lassosum.pipeline objects
  #' @keywords internal
  
  l <- list(...)

  xp.l <- l[[1]]
  if(length(l) > 1) {
    for(i in 2:length(l)) {
        # i <- 1; j <- 1
        
      # Checks
      checks <- c("lambda", "s", "keep.test", "destandardized")
      for(m in 1:length(checks)) {
        if(!identical(xp.l[[checks[m]]],
                            l[[i]][[checks[m]]])) {
          stop(paste(m, "is not the same in the",
                     paste0(i,"th"), "element of the list."))
        }
      }

      # beta
      for(k in 1:length(xp.l$beta)) {
        xp.l$beta[[k]] <- rbind(xp.l$beta[[k]],
                                     l[[i]]$beta[[k]])
      }

      # pgs
      for(k in 1:length(xp.l$beta)) {
        xp.l$pgs[[k]] <- xp.l$pgs[[k]] + l[[i]]$pgs[[k]]
      }

      # vecs
      vecs <- c("test.extract", "also.in.refpanel", "sd", "test.bfile")
      for(m in 1:length(vecs)) {
        xp.l[[vecs[m]]] <- c(xp.l[[vecs[m]]],
                                    l[[i]][[vecs[m]]])
      }

    }
    if(all(xp.l$test.bfile == xp.l$test.bfile[1])) {
      xp.l$test.bfile <- xp.l$test.bfile[1]
    }
  }
  
  return(xp.l)
}
