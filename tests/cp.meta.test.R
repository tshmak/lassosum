# clear()
# Tim.load(lassosum)
# Tim.load(Rplink)
# Tim.load(mysimtools)
# useWTCCC()
# pl <- cp.plink.linear(bfile=.bfile, chr=22, nfolds=2, 
#                       pheno=rnorm(nrow.bfile(.bfile)))
# 
# # Fake sumstats
# bim <- read.bim(.bfile)
# bim22 <- subset(bim, CHROM==22)
# n <- 1000
# samp <- sample(nrow(bim22), n)
# bim22 <- bim22[samp,]
# rev <- rbinom(n, 1, 0.5)
# bim22$A1 <- ifelse(rev, bim22$ALT, bim22$REF)
# bim22$A2 <- ifelse(rev, bim22$REF, bim22$ALT)
# sampsize <- 10000
# bim22$cor <- rnorm(n, sd=1/sampsize)
# 
# test <- with(bim22, cp.meta(pl, cor=cor, chr=CHROM, pos=POS, A1=A1, 
#                             n=sampsize))
# test2 <- cp.lassosum(test)