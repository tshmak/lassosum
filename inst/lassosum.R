#!/usr/bin/env Rscript 
Tim.load(lassosum)
# library(lassosum)
args <- parseargs()

#### Options specific to standalone version ####
# nthreads, data, n, pval, beta, OR, pheno, covar, out, lassosum.pipeline, 
# pseudovalidate, validate, splitvalidate

#### out ####
if(is.null(args[['out']])) {
  out <- "lassosum"
} else {
  out <- args[['out']]
}

#### trace ####
if(!is.null(args$trace)) {
  trace <- args$trace
} else trace <- 1

#### nthreads ####
if(!is.null(args$nthreads)) {
  cl <- parallel::makeForkCluster(args$nthreads)
} else cl <- NULL

#### data ####
if(is.null(args$data)) {
  stop("--data must be specified")
}
data <- read.table2(args$data, header=TRUE)

#### variables in data ####
vars <- c("cor", "chr", "pos", "snp", "A1", "A2", "LDblocks", 
          "pval", "beta", "OR") # Other options specific to 
for(v in (vars)) if(!is.null(args[[v]])) assign(v, data[[ args[[v]] ]]) else 
  assign(v, NULL)

# n is special, cos it can take a number or be a variable in data
if(!is.null(args[["n"]])) n <- if(is.character(args[['n']])) data[[ args[['n']] ]] else 
  args[['n']] else n <- NULL
if(is.null(beta) && !is.null(OR)) beta <- log(OR)

#### p2cor ####
if(is.null(cor)) {
  if(!is.null(pval)) {
    if(!is.null(n) && !is.null(beta)) {
      if(!is.null(args[['min.n']])) {
        cor <- p2cor(pval, n=n, sign = beta, min.n = args[['min.n']]) 
      } else {
        cor <- p2cor(pval, n=n, sign = beta) 
      }
    } else {
      stop("You need to specify n and beta/OR to convert pval to correlations.")
    }
  } else {
    stop("Either cor or pval must be provided.")
  }
}

#### data ####
if(!is.null(args[['debug']]) && interactive()) debug(lassosum.pipeline)
if(is.null(args$lassosum.pipeline)) {
  message("Running lassosum.pipeline")
  opts <- list(cor=cor, chr=chr, pos=pos, snp=snp, A1 = A1, A2=A2,
               test.bfile = args$test.bfile, ref.bfile=args$ref.bfile, 
               LDblocks=args$LDblocks, lambda=args$lambda, 
               s=args[['s']], destandardize=args$destandardize, 
               exclude.ambiguous = args$exclude.ambiguous, 
               keep.ref=args$keep.ref, remove.ref=args$remove.ref, 
               keep.test=args$keep.test, remove.test=args$remove.test,
               sample=args$sample, cluster=cl, 
               max.ref.bfile.n = args$max.ref.bfile.n, 
               trace=trace)
  opts <- opts[!sapply(opts, is.null)]
# str(opts)
  lp <- do.call(lassosum.pipeline, opts)
  saveRDS(lp, file=paste0(out, ".lassosum.pipeline.rds"))
} else {
  lp <- readRDS(args$lassosum.pipeline)
}

#### pheno ####
if(!is.null(args[['pheno']])) {
  pheno <- read.table2(args[['pheno']], header=T)
} else pheno <- NULL

#### covar ####
if(!is.null(args[['covar']])) {
  covar <- read.table2(args[['covar']], header=T)
} else covar <- NULL

if(!is.null(pheno) || !is.null(args[['validate']])) {
  message("Running validate.lassosum.pipeline")
  v <- validate(lp, pheno=pheno, covar=covar, trace=trace)
  saveRDS(v, file=paste0(out, ".validate.rds"))
  write.table2(v$results.table, file=paste0(out, ".validate.results.txt"), 
               col.names=T)
}

if(!is.null(pheno) || !is.null(args[['splitvalidate']])) {
  message("Running splitvalidate.lassosum.pipeline")
  v <- splitvalidate(lp, pheno=pheno, covar=covar, trace=trace)
  saveRDS(v, file=paste0(out, ".splitvalidate.rds"))
  write.table2(v$results.table, file=paste0(out, ".splitvalidate.results.txt"), 
               col.names=T)
}

if(!is.null(args[['pseudovalidate']])) {
  message("Running pseudovalidate.lassosum.pipeline")
  v <- pseudovalidate(lp, trace=trace)
  saveRDS(v, file=paste0(out, ".pseudovalidate.rds"))
}

