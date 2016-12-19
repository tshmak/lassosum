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
shrink <- 0.4

test_that("runElnet" , {
L <-  runElnet(lambda, shrink, plinkfile, 
         corr, n, p, 
         integer(0), integer(0), integer(0), integer(0), 
         1e-4, x, 0, 1e4, 0, p-1)
X2 <- scale(X) * sqrt(n/(n-1))
y <- scale(y) * sqrt(n/(n-1))

g <- lassosumR(corr, X, lambda, shrink)

max(abs(g$beta - L$beta))
expect_equal(g$beta, L$beta)
})


test_that("Running by blocks" , {
  split <- floor(p/2)
  block1 <- rep(F, p); block1[1:p <= split] <- T
  block2 <- !block1
  blocks <- as.integer(block1); blocks[block2] <- 2

  g1 <- lassosumR(corr[block1], X[,block1], lambda, shrink)
  g2 <- lassosumR(corr[block2], X[,block2], lambda, shrink)
  g12.beta <- rbind(g1$beta, g2$beta)

  g <- lassosumR(corr, X, lambda, shrink,blocks = blocks)
  
  max(abs(g$beta - g12.beta))
  expect_equal(g$beta, g12.beta)
})
