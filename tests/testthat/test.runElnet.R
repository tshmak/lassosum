test_that("runElnet" , {
plinkfileStem <- "test"
plinkfile <- paste0(plinkfileStem, ".bed")
n <- nrow.bfile(plinkfileStem)
p <- ncol.bfile(plinkfileStem)
beta <- rnorm(p)
Xb <- pgs(plinkfileStem, beta)
X <- genotypeMatrix(plinkfile, n, p, 
                    integer(0), integer(0), integer(0), integer(0), 1)
Xb2 <- X %*% beta
stopifnot(all.equal(Xb2, Xb))
y <- Xb + rnorm(n)
corr <- cor(y, X)

lambda <- sort(exp(seq(log(0.001), log(0.1), length=20)), decreasing=T)
x <- rep(0.0, p)
L <-  runElnet(lambda, 0, plinkfile, 
         corr, n, p, 
         integer(0), integer(0), integer(0), integer(0), 
         1e-4, x, 0, 1e4, 0, p-1)
# L <-  runElnet(lambda, 0, plinkfile, 
#                corr, n, p, 
#                integer(0), integer(0), integer(0), integer(0), 
#                1e-4, x, 0, 1e4, c(0,5), c(4,p-1))
X2 <- scale(X) * sqrt(n/(n-1))
y <- scale(y) * sqrt(n/(n-1))
g <- glmnet::glmnet(X2, y = y,alpha=1, lambda = lambda, standardize = F, intercept = F, 
            thresh=1e-10)
g.beta <- as.matrix(g$beta)
max(abs(g.beta - L$beta))
})
