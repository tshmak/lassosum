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

test_that("lassosum" , {
  ### Testing blocks... ###
  split <- floor(p/2)
  block1 <- rep(F, p); block1[1:p <= split] <- T
  block2 <- !block1
  blocks <- as.integer(block1); blocks[block2] <- 2
  
  g1 <- lassosum(corr[block1], plinkfileStem, lambda, shrink, extract=block1)
  g2 <- lassosum(corr[block2], plinkfileStem, lambda, shrink, extract=block2)
  g12 <- merge.lassosum(g1, g2)
  
  g <- lassosum(corr, plinkfileStem, lambda, shrink,blocks = blocks)
  
  max(abs(g$beta - g12$beta))
  expect_equal(g$beta, g12$beta)
})


