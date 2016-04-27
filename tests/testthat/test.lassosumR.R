test_that("runElnet" , {
plinkfileStem <- "test"
plinkfile <- paste0(plinkfileStem, ".bed")
n <- nrow.bfile(plinkfileStem)
p <- ncol.bfile(plinkfileStem)
beta <- rnorm(p)
Xb <- pgs(plinkfileStem, beta)
X <- genotypeMatrix(plinkfile, n, p, 
                    integer(0), integer(0), integer(0), integer(0))
Xb2 <- X %*% beta
stopifnot(all.equal(Xb2, Xb))
y <- Xb + rnorm(n)
corr <- cor(y, X)

lambda <- sort(exp(seq(log(0.001), log(0.1), length=20)), decreasing=T)
x <- rep(0.0, p)
shrink <- 0.4
L <-  runElnet(lambda, shrink, plinkfile, 
         corr, n, p, 
         integer(0), integer(0), integer(0), integer(0), 
         1e-4, x, 0, 1e4)
X2 <- scale(X) * sqrt(n/(n-1))
y <- scale(y) * sqrt(n/(n-1))

g <- lassosumR(corr, X, lambda, shrink)

max(abs(g$beta - L$beta))
expect_equal(g$beta, L$beta)
})
