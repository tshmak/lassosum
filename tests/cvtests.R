# clear()
# Tim.load(Rplink)
# setwd(attachroot("~/WORK/myRpackages/lassosum/tests/"))
# load_all()
# 
# library(data.table)
# bfile("../inst/data/refpanel")
# n <- nrow.bfile(.bfile)
# set.seed(1000)
# pheno <- rnorm(n)
# nfolds <- 3
# ss <- cp.plink.linear(bfile=.bfile, nfolds=nfolds, pheno=pheno)
# 
# ld <- fread("../inst/data/Berisa.EUR.hg19.bed")
# cp <- cp.lassosum(ss, LDblocks=ld, details=TRUE, Method2=TRUE)
# v <- attr(cp, "validate")
# beta <- cp.beta(ss, cp, save="cvtests.pinv")
# beta <- cp.beta(ss, cp, load="cvtests.pinv")
# pred <- pgs(ss$bfile, beta)
# diff <- cp$best.pgs - pred
# summary(diff)
# file.remove("cvtests.pinv")
