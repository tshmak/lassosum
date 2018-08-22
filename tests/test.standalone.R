# Tmisc()
# Tim.load(lassosum)
# exec <- attachroot("~/WORK/myRpackages/lassosum/inst/lassosum.R")
# system(paste("chmod +x", exec))
# # system(exec)
# setwd(tempdir())
# f <- function(file="") paste0(system.file("data", package="lassosum"), "/", file)
# ss <- read.table2(f("summarystats.txt"), header=TRUE)
# options <- list(data=f("summarystats.txt"), chr="Chr", pos="Position",
#      A1="A1", A2="A2", pval="P_val", n=50000, OR="OR_A1",
#      test.bfile=f("testsample"), LDblocks="EUR.hg19",
#      validate=T, splitvalidate=T, pseudovalidate=T, nthreads=2)
# cmd <- parse.options(options = options, dot2dash=F)
# # opts <- paste(exec, cmd); args <- parseargs(test=opts, include.others=F); source.some(exec, start="opts <- args", run=T)
# # source(exec)
# system(paste(exec, cmd))
# head(read.table2("lassosum.validate.results.txt", header=T))
# 
# Tim.load(Rplink)
# fam <- read.fam(f("testsample"))
# fam2 <- fam[sample(nrow(fam), 150),]
# write.table2(fam2[,1:2], file=keep <- tempfile())
# options$keep.test <- keep
# system(paste(exec, parse.options(options = options, dot2dash=F)))
# lp <- readRDS("lassosum.lassosum.pipeline.rds")
# 
# fam3 <- fam[sample(nrow(fam), 170),]
# write.table2(fam3[,c(1:2,6)], col.names=TRUE, file=PHENO <- tempfile())
# options$lassosum.pipeline <- "lassosum.lassosum.pipeline.rds"
# options$pheno <- PHENO
# system(paste(exec, parse.options(options = options, dot2dash=F)))
# 
# fam4 <- fam[sample(nrow(fam), 170),]
# covar <- cbind(fam4[,1:2], rnorm(nrow(fam4)), rnorm(nrow(fam4)))
# write.table2(covar, col.names=TRUE, file=COVAR <- tempfile())
# options$covar <- COVAR
# system(paste(exec, parse.options(options = options, dot2dash=F)))
# 
# #### Apply to new data ####
# options2 <- list(lassosum.pipeline="lassosum.lassosum.pipeline.rds", 
#                  validate.rds="lassosum.pseudovalidate.rds", 
#                  applyto=T, keep=keep)
# cmd <- parse.options(options2, dot2dash = F)
# system(paste(exec, cmd))
# 
# options2 <- list(lassosum.pipeline="lassosum.lassosum.pipeline.rds", 
#                  validate.rds="lassosum.validate.rds", 
#                  applyto="/home/tshmak/WORK/myRpackages/lassosum/inst/data/refpanel")
# cmd <- parse.options(options2, dot2dash = F)
# system(paste(exec, cmd))
# # opts <- paste(exec, cmd); args <- parseargs(test=opts, include.others=F); source.some(exec, start="opts <- args", run=T)