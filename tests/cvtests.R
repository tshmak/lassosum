# clear()
# Tim.load(Rplink)
# setwd(attachroot("/WORK/Projects/validation/lassosum/tests/"))
# load_all()
# 
# library(data.table)
# bfile("../tutorial/data/refpanel")
# n <- nrow.bfile(.bfile)
# set.seed(1000)
# pheno <- rnorm(n)
# nfolds <- 3
# ss <- xp.plink.linear(bfile=.bfile, nfolds=nfolds, pheno=pheno)
# 
# ld <- fread("../tutorial/data/Berisa.2015.EUR.bed")
# xp <- xp.lassosum(ss, LDblocks=ld, details=TRUE, Type2=TRUE)
# v <- attr(xp, "validate")
# beta <- xp.beta(ss, xp, save="cvtests.pinv")
# beta <- xp.beta(ss, xp, load="cvtests.pinv")
# pred <- pgs(ss$bfile, beta)
# diff <- xp$best.pgs - pred
# summary(diff)
# file.remove("cvtests.pinv")
