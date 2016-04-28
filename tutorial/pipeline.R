# setwd("./tutorial")
### Use data.table for speed ###
library(data.table)

### Read summary statistics file ###
ss <- fread("./data/summarystats.txt", data.table=F)

### Read .bim file ###
bim <- fread("./data/chr22a.bim", data.table=F)

### Select chromosome 22 only ###
ss.chr22 <- subset(ss, CHR==22) 	# Note that we recommend analyses chromosome by chromosome

### Compare ss and bim 
comp <- comp.ss.bim(ss.chr22[, c("SNP", "A2", "A1")], bim[, c("V2", "V6", "V5")]) 
### get vector of correlations from p-values ###
correlation <- with(ss.chr22, 
		    p2cor(p=P[comp$ss.order], 
			  n = NMISS[comp$ss.order], 
			  sign=log(OR[comp$ss.order]) * comp$rev))

### Choose lambda values ###
lambda <- exp(seq(log(0.001), log(0.1), length.out=20))

### Get beta estimates from lassosum ###
ls <- lassosum(cor=correlation, bfile="./data/chr22a", lambda=lambda, shrink=0.9, 
	       extract=comp$bim.extract)

### Get indeplasso estimates for SNPs that are not in reference 
correlation2 <- with(ss.chr22, p2cor(P, NMISS))
il <- indeplasso(correlation2, lambda = lambda)

### Get beta from both lasso and indeplasso 
beta <- il$beta
beta[comp$ss.order, ] <- ls$beta

######### pseudovalidation ###########
### Get shrunken estimates of cor
fdr <- fdrtool::fdrtool(correlation2, statistic="correlation", plot=F)
correlation2.shrunk <- correlation2 * (1 - fdr$lfdr)

### Download pseudovalidation .bim file ###
val.bim <- fread("./data/chr22b.bim", data.table=F) 

### Compare ss and val.bim 
comp2 <- comp.ss.bim(ss.chr22[, c("SNP", "A2", "A1")], val.bim[, c("V2", "V6", "V5")]) 

### pseudovalidation ###
pv <- pseudovalidation("./data/chr22b", 
                       beta=beta[comp2$ss.order, ] * outer(comp2$rev, rep(1,ncol(beta))), 
                       cor=correlation2.shrunk[comp2$ss.order], 
                       extract=comp2$bim.extract)
plot(lambda, pv, log="x")
best.lambda.pos <- which(pv == max(pv))
cat("Best lambda = ", lambda[best.lambda.pos], "\n")

### Obtain PGS for best lambda ###
PGS <- pgs(bfile = "./data/chr22b", weights = beta[comp2$ss.order] * comp2$rev, 
	   extract=comp2$bim.extract)

######### lassosum using a reference panel given as a matrix ###########
### load "./data/chr22a" as a matrix ###
chr22a <- readbfile("./data/chr22a", fillmissing=T)

### Get beta estimates from lassosum ###
lsR <- lassosumR(cor=correlation, refpanel=chr22a[,comp$bim.extract], 
                 lambda=lambda, shrink=0.9) 

### pseudovalidation using chr22b as target data ###
chr22b <- readbfile("./data/chr22b", fillmissing=T)
pvR <- pseudovalidationR(chr22b[, comp2$bim.extract], 
                        beta=beta[comp2$ss.order, ] * outer(comp2$rev, rep(1,ncol(beta))), 
                        cor=correlation2.shrunk[comp2$ss.order])

######### De-standardizing correlation coefficients to get regression coefficients ###########
### Obtain SNP-wise standard deviation of target dataset ###
sd <- sd.bfile(bfile = "./data/chr22b", extract=comp2$bim.extract)

### regression coefficients = correlation coefficients / sd(X) * sd(y) ###
reg.coef <- Matrix::Diagonal(x=1/sd) %*% beta[comp2$ss.order, ]

