clear()
load_all()
Tim.load(Rplink)
setwd(attachroot("/WORK/myRpackages/lassosum/tests/"))

library(data.table)
bfile("../tutorial/data/refpanel")
n <- nrow.bfile(.bfile)
set.seed(1000)
pheno <- rnorm(n)
nfolds <- 3
ss <- xp.plink.linear(bfile=.bfile, nfolds=nfolds, pheno=pheno)

ld <- fread("../tutorial/data/Berisa.2015.EUR.bed")
xp <- xp.lassosum(ss, LDblocks=ld, details=TRUE)
v <- attr(xp, "validate")