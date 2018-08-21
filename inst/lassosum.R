#!/usr/bin/env Rscript 
# Tim.load(lassosum)
library(lassosum)
args <- parseargs(required="data", include.others = FALSE)
opts <- args

#### Options specific to standalone version ####
standalone.opts <- c("nthreads", "data", "n", "pval", "beta", "OR", "pheno", "covar", 
                     "out", "lassosum.pipeline", "pseudovalidate", "validate", 
                     "splitvalidate", "debug")
for(i in standalone.opts) opts[[i]] <- NULL

#### out ####
if(is.null(args[['out']])) {
  out <- "lassosum"
} else {
  out <- args[['out']]
}

#### nthreads ####
if(!is.null(args$nthreads)) {
  opts$cl <- parallel::makeForkCluster(args$nthreads)
} 

#### lambda ####
if(!is.null(args[['lambda']]) && is.character(args[['lambda']])) {
  opts[['lambda']] <- as.numeric(strsplit(args[['lambda']], split=",")[[1]])
} 

#### s ####
if(!is.null(args[['s']]) && is.character(args[['s']])) {
  opts[['s']] <- as.numeric(strsplit(args[['s']], split=",")[[1]])
} 

#### data ####
data <- read.table2(args[['data']], header=TRUE)

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
if(is.null(args$lassosum.pipeline)) {
  message("Running lassosum.pipeline")
# str(opts)
  lp <- do.call(lassosum.pipeline, opts)
  saveRDS(lp, file=paste0(out, ".lassosum.pipeline.rds"))
} else {
  lp <- readRDS(args$lassosum.pipeline)
}

#### Options for validation ####
opts2 <- list(ls.pipeline = lp)
opts2$pheno <- args[['pheno']]
opts2$covar <- args[['covar']]
opts2 <- c(opts2, opts)

opts2.tokeep <- list(validate=names(as.list(args(validate.lassosum.pipeline))),
                     splitvalidate=names(as.list(args(splitvalidate.lassosum.pipeline))),
                     pseudovalidate=names(as.list(args(pseudovalidate.lassosum.pipeline))))

if(!is.null(args[['pheno']])) {
  args[['validate']] <- args[['splitvalidate']] <- TRUE
}
for(i in 1:length(opts2.tokeep)) {
  type <- names(opts2.tokeep)[i]
  tokeep <- names(opts2) %in% opts2.tokeep[[i]]
  
  if(!is.null(args[[type]])) {
    message("Running ", type, ".lassosum.pipeline")
    if(!is.null(args[['debug']]) && interactive()) debug(paste0(type,".lassosum.pipeline"))
    v <- do.call(type, opts2[tokeep])
    if(type %in% c("validate", "splitvalidate")) 
      saveRDS(v, file=paste(out, type, "rds", sep="."))
    write.table2(v$results.table, file=paste(out, type, "results.txt", sep="."), 
                 col.names=T)
  }
}

