# Tmisc()
# Tim.load(Rplink)
# setwd(attachroot("~/WORK/myRpackages/lassosum/tests/"))
# Tim.load(lassosum)
# 
# library(data.table)
# bfile <- "../inst/data/testsample"
# n <- nrow.bfile(bfile)
# ss <- read.table2("../inst/data/summarystats.txt", header=T)
# cor <- p2cor(ss$P_val, n=60000)
# lp <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position, A1 = ss$A1, A2=ss$A2,
#                         test.bfile = bfile, LDblocks="EUR.hg19")
# pheno <- rnorm(n)
# v <- validate(lp, pheno=pheno)
# fam <- read.fam(bfile)
# pheno2 <- cbind(fam[,1:2], pheno)
# v2 <- validate(lp, pheno=pheno2)
# sv <- splitvalidate(lp, pheno=pheno2)
# s <- sample(n, round(n*0.8))
# pheno3 <- pheno2[s, ]
# covar <- rnorm(n)
# covar2 <- cbind(fam[,1:2], covar, rnorm(n))
# s <- sample(n, round(n*0.8))
# covar3 <- covar2[s, ]
# 
# v <- validate(lp, pheno=pheno, covar=covar)
# v2 <- validate(lp, pheno=pheno2, covar=covar2)
# v3 <- validate(lp, pheno=pheno3, covar=covar3)
# print(rbind(v3$best.pgs[v3$results.table$order], v3$results.table$best.pgs))
# 
# sv <- splitvalidate(lp, pheno=pheno3, covar=covar3)
# 
# #### Testing subsetting ####
# lp2 <- subset(lp, lambda=v$best.lambda, s=v$best.s)
# vv <- validate(lp2)
