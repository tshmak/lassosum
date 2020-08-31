#!/usr/bin/env Rscript 
# Tim.load(lassosum)
library(lassosum)
args <- lassosum:::parseargs(include.others = FALSE)
opts <- args

#### Options specific to standalone version ####
standalone.opts <- c("nthreads", "data", "n", "pval", "beta", "OR", "pheno", "covar", 
                     "out", "lassosum.pipeline", "pseudovalidate", "validate", 
                     "splitvalidate", "validate.rds", "applyto", "seed", 
                     "debug")
for(i in standalone.opts) opts[[i]] <- NULL

#### seed ####
if(!is.null(args[['seed']])) {
  set.seed(args[['seed']])
} else {
  set.seed(1234)
}

#### applyto ####
if(!is.null(args[['applyto']])) {
  if(is.null(args[['lassosum.pipeline']])) 
    stop("--lassosum.pipeline must be specified with --applyto")
  if(is.null(args[['validate.rds']])) 
    stop("--validate.rds must be specified with --applyto")
}

#### out ####
if(is.null(args[['out']])) {
  out <- "lassosum"
} else {
  out <- args[['out']]
}

#### nthreads ####
if(!is.null(args[['nthreads']])) {
  opts[['cluster']] <- parallel::makeForkCluster(args[['nthreads']])
} 

parsefun <- function(input) {
  if(is.character(input)) {
    s <- strsplit(input, split=",")[[1]]
    suppressWarnings(ss <- as.numeric(s))
    if(!all(is.na(ss))) return(ss) else return(s)
  } else return(input)
}

#### test.bfile ####
if(!is.null(args[['test.bfile']])) {
  opts[['test.bfile']] <- parsefun(args[['test.bfile']])
} 

#### ref.bfile ####
if(!is.null(args[['ref.bfile']])) {
  opts[['ref.bfile']] <- parsefun(args[['ref.bfile']])
} 

if(is.null(args[['lassosum.pipeline']])) {

  #### lambda ####
  if(!is.null(args[['lambda']])) {
    opts[['lambda']] <- parsefun(args[['lambda']])
  } 
  
  #### s ####
  if(!is.null(args[['s']])) {
    opts[['s']] <- parsefun(args[['s']])
  } 
  
  #### data ####
  if(is.null(args[['data']])) stop("--data is required.")
  data <- lassosum:::read.table2(args[['data']], header=TRUE)  

  #### variables in data ####
  vars <- c("cor", "chr", "pos", "snp", "A1", "A2")
  for(v in vars) {
    if(!is.null(args[[v]])) {
      if(!is.null(data[[ args[[v]] ]])) opts[[v]] <- data[[ args[[v]] ]] else
        stop(paste("Cannot find", args[[v]], "in data"))
    } 
  }
  
  #### Others variables in data ####
  vars <- c("pval", "beta", "OR")
  for(v in vars) {
    if(!is.null(args[[v]])) {
      if(!is.null(data[[ args[[v]] ]])) assign(v, data[[ args[[v]] ]]) else 
        stop(paste("Cannot find", args[[v]], "in data"))
    } else assign(v, NULL)
  }
  
  #### n ####
  # n is special, cos it can take a number or be a variable in data
  if(!is.null(args[["n"]])) n <- if(is.character(args[['n']])) data[[ args[['n']] ]] else 
    args[['n']] else n <- NULL
  
  #### LDblocks ####
  if(!is.null(args[["LDblocks"]])) opts[['LDblocks']] <- if(!is.null(data[[ args[['LDblocks']] ]])) 
    data[[ args[['LDblocks']] ]] else opts[['LDblocks']] <- args[['LDblocks']]
  
  #### min.n ####
  if(!is.null(args[["min.n"]])) min.n <- args[["min.n"]] else min.n <- NULL
  
  #### p2cor ####
  if(is.null(beta) && !is.null(OR)) beta <- log(OR)
  if(is.null(opts[['cor']])) {
    if(!is.null(pval)) {
      if(!is.null(n) && !is.null(beta)) {
        if(!is.null(min.n)) {
          opts[['cor']] <- p2cor(pval, n=n, sign = beta, min.n = min.n) 
        } else {
          opts[['cor']] <- p2cor(pval, n=n, sign = beta) 
        }
      } else {
        stop("You need to specify n and beta/OR to convert pval to correlations.")
      }
    } else {
      stop("Either cor or pval must be provided.")
    }
  }
  
  #### debug ####
  if(!is.null(args[['debug']]) && interactive()) debug(lassosum.pipeline)
  
  #### Run lassosum.pipeline ####
  message("\nRunning lassosum.pipeline")
  # str(opts)
  lp <- do.call(lassosum.pipeline, opts)
  saveRDS(lp, file=paste0(out, ".lassosum.pipeline.rds"))
  
} else {
  lp <- readRDS(args[['lassosum.pipeline']])
}

#### Options for validation ####
opts2 <- list(ls.pipeline = lp)
opts2[['pheno']] <- args[['pheno']]
opts2[['covar']] <- args[['covar']]
opts2 <- c(opts2, opts)

opts2.tokeep <- list(pseudovalidate=names(as.list(args(lassosum:::pseudovalidate.lassosum.pipeline))), 
                     splitvalidate=names(as.list(args(lassosum:::splitvalidate.lassosum.pipeline))),
                     validate=names(as.list(args(lassosum:::validate.lassosum.pipeline))))

if(is.null(args[['validate.rds']])) {
  if(!is.null(args[['pheno']])) {
    args[['validate']] <- args[['splitvalidate']] <- TRUE
  }
  for(i in 1:length(opts2.tokeep)) {
    type <- names(opts2.tokeep)[i]
    tokeep <- names(opts2) %in% opts2.tokeep[[i]]
    
    if(!is.null(args[[type]])) {
      message("\nRunning ", type, ".lassosum.pipeline")
      if(!is.null(args[['debug']]) && interactive()) debug(paste0(type,".lassosum.pipeline"))
      v <- do.call(type, opts2[tokeep])
      saveRDS(v, file=paste(out, type, "rds", sep="."))
      lassosum:::write.table2(v[['results.table']], file=paste(out, type, "results.txt", sep="."), 
                              col.names=T)
    }
  }
} else {
  v <- readRDS(args[['validate.rds']])
}

#### Apply to new data ####
if(!is.null(args[['applyto']])) {
  if(is.character(args[['applyto']])) opts2[['test.bfile']] <- args[['applyto']]
  tokeep <- names(opts2) %in% opts2.tokeep[['validate']]
  opts2[['ls.pipeline']] <- subset(lp, s=v[['best.s']], lambda=v[['best.lambda']])
  v2 <- do.call("validate", opts2[tokeep])
  lassosum:::write.table2(v2[['results.table']], file=paste(out, "results.txt", sep="."), 
                          col.names=T)
}
