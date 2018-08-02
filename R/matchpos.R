matchpos <- function(tomatch, ref.df, 
                     auto.detect.tomatch=T, auto.detect.ref=T, 
                     chr="", ref.chr="", 
                     pos="", ref.pos="",
                     snp="", ref.snp="",
                     ref="", ref.ref="",
                     alt="", ref.alt="",
                     exclude.ambiguous=F, 
                     silent=F, rm.duplicates=F) {

  #' @title Function to match a set of variants to a reference
  #' by chromosome and/or position and/or SNP/rsIDs 
  #' @description 
  #' If REF/ALT alleles available, will also check them for reverse coding. 
  #' Specify the exclude.ambiguous option for possible sense/antisense ambiguity. 
  #' @param tomatch A data.frame or data.table
  #' @param ref.df A data.frame or data.table
  #' @param rm.dumplicates Remove SNPs with more than one match
  #' @keywords internal 

  tomatch <- data.table::as.data.table(tomatch)
  ref.df <- data.table::as.data.table(ref.df)
    
  #### Column matching function ####
  match.col <- function(test, target, default.targets, auto.detect=T, silent=F) {
    mustfind <- target != ""
    if(mustfind) {
      test2 <- test
      targets2 <- target
    } else if(auto.detect) {
      test2 <- tolower(substr(test,1,3))
      targets2 <- tolower(substr(default.targets,1,3))
      target <- default.targets[1]
    } else return(integer(0))
    
    pos <- test2 == targets2[1]
    if(length(targets2) > 1) {
      for(i in 2:length(targets2)) {
        pos <- pos | test2 == targets2[i]
      }
    }
    pos <- which(pos)
    if(length(pos) > 1) {
      vars <- paste(test[pos], collapse = ", ")
      mess <- paste0("Multiple columns for ", target, ". Possible matches: ", vars)
      stop(mess)
    } else if(length(pos) == 1) {
      if(!silent) cat(paste(test[pos], "identified as", target, "\n")) 
    } else if(mustfind) {
      stop(paste("No column identified for ", target))
    } else {
      if(!silent) cat("No column identified for ", target, "\n")
    }
    
    return(pos)
  }
  
  #### Identify columns ####
  colnames.tomatch <- colnames(tomatch)
  if(!silent) cat("For 'tomatch' data.frame: \n")
  chr.col <- match.col(colnames.tomatch, chr, c("CHROMOSOME"), auto.detect.tomatch, silent=silent)
  pos.col <- match.col(colnames.tomatch, pos, c("POSITION", "BP"), auto.detect.tomatch, silent=silent)
  snp.col <- match.col(colnames.tomatch, snp, c("SNP", "rsid"), auto.detect.tomatch, silent=silent)
  alt.col <- match.col(colnames.tomatch, alt, c("alt", "a1"), auto.detect.tomatch, silent=silent)
  ref.col <- match.col(colnames.tomatch, ref, c("ref", "a2"), auto.detect.tomatch, silent=silent)
  
  colnames.ref <- colnames(ref.df)
  if(!silent) cat("\nFor 'reference' data.frame: \n")
  ref.chr.col <- match.col(colnames.ref, ref.chr, c("CHROMOSOME"), auto.detect.ref, silent=silent)
  ref.pos.col <- match.col(colnames.ref, ref.pos, c("POSITION", "BP"), auto.detect.ref, silent=silent)
  ref.snp.col <- match.col(colnames.ref, ref.snp, c("SNP", "rsid"), auto.detect.ref, silent=silent)
  ref.alt.col <- match.col(colnames.ref, ref.alt, c("alt", "a1"), auto.detect.ref, silent=silent)
  ref.ref.col <- match.col(colnames.ref, ref.ref, c("ref", "a2"), auto.detect.ref, silent=silent)
  
  
  #### Find the matching columns ####
  match.snp <- length(snp.col) > 0 && length(ref.snp.col) > 0
  match.chr.pos <- length(chr.col) > 0 && length(pos.col) > 0 && 
    length(ref.chr.col) > 0 && length(ref.pos.col) > 0 
  if(!match.snp && !match.chr.pos) {
    stop("We must be able to match either by variant ID (SNP) or by chromosome and position")
  } else if(match.chr.pos && match.snp) {
    match.cols <- c(snp.col, chr.col, pos.col) 
    ref.match.cols <- c(ref.snp.col, ref.chr.col, ref.pos.col)
  } else if(match.chr.pos) {
    match.cols <- c(chr.col, pos.col) 
    ref.match.cols <- c(ref.chr.col, ref.pos.col)
  } else {
    match.cols <- c(snp.col) 
    ref.match.cols <- c(ref.snp.col)
  }
    
  n.match.cols <- length(match.cols)
  match.cols.names <- colnames(tomatch)[match.cols]
  ref.match.cols.names <- colnames(ref.df)[ref.match.cols]
  
#   all.cols <- c(snp.col, chr.col, pos.col, alt.col, ref.col)
#   ref.all.cols <- c(ref.snp.col, ref.chr.col, ref.pos.col, ref.alt.col, ref.ref.col)
  
#   ncols.tomatch <- length(all.cols)
#   ncols.ref <- length(ref.all.cols)
  
  tomatch.ref <- length(ref.col) > 0
  tomatch.alt <- length(alt.col) > 0
  refer.ref <- length(ref.ref.col) > 0
  refer.alt <- length(ref.alt.col) > 0
  
  ref.alt.tomatch <- tomatch.ref && tomatch.alt
  ref.alt.ref <- refer.ref && refer.alt
  
  match.both.alleles <- ref.alt.tomatch && ref.alt.ref
  
  #### Matching by one allele #### 
  match1allele <- ""
  if(!match.both.alleles) {
    if(tomatch.ref && refer.ref) {
      match1allele <- "ref"
    } else if(tomatch.alt && refer.alt) {
      match1allele <- "alt"
    } else if(tomatch.ref && refer.alt) {
      match1allele <- "ref.alt"
    } else if(tomatch.alt && refer.ref) {
      match1allele <- "alt.ref"
    }
  }
  
  # if(match1allele != "") warning("Matching on 1 allele only!")
  
  tomatch$.index.tomatch <- 1:nrow(tomatch)
  ref.df$.index.ref <- 1:nrow(ref.df)

  #### MERGE ####
  merged <- merge(ref.df, tomatch, all=F, 
                  by.x=ref.match.cols.names, by.y=match.cols.names)
  setkey(merged, .index.ref)
  merged <- as.data.frame(merged)

  alt.col2 <- alt.col + ncol(ref.df) - sum(match.cols < alt.col)
  ref.col2 <- ref.col + ncol(ref.df) - sum(match.cols < ref.col)
  ref.alt.col2 <- sum(ref.match.cols > ref.alt.col) + ref.alt.col
  ref.ref.col2 <- sum(ref.match.cols > ref.ref.col) + ref.ref.col
  
  #### Change everything to capital letters ####
  if(tomatch.alt) merged[,alt.col2] <- toupper(merged[,alt.col2])
  if(tomatch.ref) merged[,ref.col2] <- toupper(merged[,ref.col2])
  if(refer.alt) merged[,ref.alt.col2] <- toupper(merged[,ref.alt.col2])
  if(refer.ref) merged[,ref.ref.col2] <- toupper(merged[,ref.ref.col2])

  #### Convert T -> A and G -> C if excluding ambiguous SNPs ####
  merged$ok123 <- T
  if(exclude.ambiguous) {
    if(match.both.alleles) {
      ### Ensure both ref and alt alleles switch at the same time ###
      ok1 <-  (merged[, ref.col2] == merged[, ref.ref.col2] & 
               merged[, alt.col2] == merged[, ref.alt.col2] ) 
      ok2 <- (merged[, ref.col2] == merged[, ref.alt.col2] & 
                merged[, alt.col2] == merged[, ref.ref.col2] )
      ok3 <- (merged[, ref.col2] != merged[, ref.ref.col2] & 
                merged[, ref.col2] != merged[, ref.alt.col2] & 
                merged[, alt.col2] != merged[, ref.ref.col2] & 
                merged[, alt.col2] != merged[, ref.alt.col2])
      merged$ok123 <- ok1 | ok2 | ok3
    }
    ### Convert T to A and G to C ###
    if(tomatch.alt) {
      merged[,alt.col2] <- gsub("T", "A", merged[,alt.col2])
      merged[,alt.col2] <- gsub("G", "C", merged[,alt.col2])
    }
    if(tomatch.ref) {
      merged[,ref.col2] <- gsub("T", "A", merged[,ref.col2])
      merged[,ref.col2] <- gsub("G", "C", merged[,ref.col2])
    }
    if(refer.alt) {
      merged[,ref.alt.col2] <- gsub("T", "A", merged[,ref.alt.col2])
      merged[,ref.alt.col2] <- gsub("G", "C", merged[,ref.alt.col2])
    }
    if(refer.ref) {
      merged[,ref.ref.col2] <- gsub("T", "A", merged[,ref.ref.col2])
      merged[,ref.ref.col2] <- gsub("G", "C", merged[,ref.ref.col2])
    }
  }
  
  #### Matching alleles ####
  if(match.both.alleles) {
    merged$ok <- merged[,ref.ref.col2] == merged[,ref.col2] & merged[,ref.alt.col2] == merged[,alt.col2]
    merged$antiok <- merged[,ref.ref.col2] == merged[,alt.col2] & merged[,ref.alt.col2] == merged[,ref.col2]
    merged$eitherok <- merged$ok | merged$antiok
  } else if(match1allele != "") {
    if(match1allele == "ref")  {
      merged$ok <- merged[, ref.ref.col2] == merged[, ref.col2]
      if(tomatch.alt) {
        merged$antiok <- merged[, ref.ref.col2] == merged[, alt.col2]
      } else if(refer.alt) {
        merged$antiok <- merged[, ref.alt.col2] == merged[, ref.col2]
      } else {
        merged$antiok <- !merged$ok
      }
    } else if(match1allele == "alt") {
      merged$ok <- merged[, ref.alt.col2] == merged[, alt.col2]
      if(tomatch.ref) {
        merged$antiok <- merged[, ref.alt.col2] == merged[, ref.col2]
      } else if(refer.ref) {
        merged$antiok <- merged[, ref.ref.col2] == merged[, alt.col2]
      } else {
        merged$antiok <- !merged$ok
      }
    } else if(match1allele == "ref.alt") {
      merged$antiok <- merged[, ref.col2] == merged[, ref.alt.col2]
      merged$ok <- !merged$antiok
    } else if(match1allele == "alt.ref") {
      merged$antiok <- merged[, alt.col2] == merged[, ref.ref.col2]
      merged$ok <- !merged$antiok
    } 
    merged$eitherok <- merged$ok | merged$antiok
  } else {
    merged$eitherok <- T
  }
  
  #### Exclude ambiguous SNPs ####
  # exclude.ambiguous <- F
  merged$ambiguous <- F
  if(exclude.ambiguous) {
    if(ref.alt.tomatch) {
      merged$ambiguous <- merged$ambiguous | merged[,ref.col2] == merged[,alt.col2] 
    } 
    if(ref.alt.ref) {
      merged$ambiguous <- merged$ambiguous | merged[,ref.ref.col2] == merged[,ref.alt.col2]
    }
  } 
  #### Find duplicated matches ####
  toinclude <- merged$eitherok & !merged$ambiguous & merged$ok123
  mismatches <- sort(unique(merged$.index.tomatch[!toinclude]))
  merged <- merged[toinclude,]
  dup.ref <- duplicated(merged$.index.ref)
  dup.tomatch <- duplicated(merged$.index.tomatch)
  
  if(any(dup.ref) || any(dup.tomatch)) {
    if(rm.duplicates) {
      which.index.ref <- unique(merged$.index.ref[dup.ref])
      duplicated <- merged$.index.ref %in% which.index.ref
      which.index.tomatch <- unique(merged$.index.tomatch[dup.tomatch])
      duplicated <- duplicated | merged$.index.tomatch %in% which.index.tomatch
      merged <- merged[!duplicated,]
      mismatches <- sort(unique(c(mismatches, merged$.index.tomatch[duplicated])))

    } else {
      print(merged[dup.ref | dup.tomatch, ])
      stop("Multiple matches found for some of the keys")
    }
  }

  #### Reverse allele coding ####
  rev <- NULL
  if(match.both.alleles || match1allele != "") {
    rev <- ifelse(merged$antiok, -1,1)
  } 
  
  #### order and ref.extract ####
  order <- merged$.index.tomatch
  ref.extract <- rep(F, nrow(ref.df))
  ref.extract[merged$.index.ref] <- T
  
  
  return(list(order=order, ref.extract=ref.extract, rev=rev, 
              mismatches=mismatches))

}