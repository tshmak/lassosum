#library(glmnet::glmnet)
library(testthat)

n <- 1000
p <- 20
# set.seed(1000)
X <- matrix(rnorm(n*p), n, p)
beta <- rnorm(p)
Xb <- X %*% beta
y <- Xb + rnorm(n) * p
y <- scale(y) * sqrt(n/(n-1))
X <- scale(X) * sqrt(n/(n-1))
r <- as.vector(t(X) %*% y) / n
summary(r)
summary(colSums(X^2))
summary(colSums(X))
summary(y)
sum(y^2)

test_that("LASSO comparison", {

lambda <- 0.01
alpha <- 1
lambda1 <- lambda * alpha
lambda2 <- lambda * (1-alpha)
gg <- glmnet::glmnet(X, y, standardize=F, intercept=F, lambda=lambda, alpha=alpha, thresh=1e-10)
as.vector(gg$beta)

x <- rep(0.0, p)
shrink <- 1.0
X.touse <- X /sqrt(n)* sqrt(shrink)
yhat <- X.touse %*% x
diag <- colSums(X.touse^2)
elnet(lambda1=lambda1, lambda2=lambda2, diag=diag, X=X.touse, r=r,
      1e-4,x,yhat, trace=1,100)

expect_equal(as.vector(x), as.vector(gg$beta), tol=1e-4, scale=1)
})

test_that(" elnet comparison" , {
lambda <- 0.01
alpha <- 0.7
lambda1 <- lambda * alpha
lambda2 <- lambda * (1-alpha)
gg <- glmnet::glmnet(X, y, standardize=F, intercept=F, lambda=lambda, alpha=alpha, thresh=1e-10)
as.vector(gg$beta)

x <- rep(0.0, p)
shrink <- 1.0
X.touse <- X /sqrt(n)* sqrt(shrink)
yhat <- X.touse %*% x
diag <- colSums(X.touse^2)
elnet(lambda1=lambda1, lambda2=lambda2, diag=diag, X=X.touse, r=r,
          1e-4,x,yhat, trace=1,100)
expect_equal(as.vector(x), as.vector(gg$beta), tolerance=1e-4, scale=1)
})

test_that("elnet test with shrinkage", {
lambda <- 0.01
alpha <- 0.7
lambda1 <- lambda * alpha
lambda2 <- lambda * (1-alpha)

x <- rep(0.0, p)
shrink <- 0.4
X.touse <- X /sqrt(n)* sqrt(shrink)
yhat <- X.touse %*% x
diag <- colSums(X.touse^2)
elnet(lambda1=lambda1, lambda2=lambda2, diag=diag, X=X.touse, r=r,
      1e-4,x,yhat, trace=1, 100)

gg <- glmnet::glmnet(X, y / shrink, standardize=F, intercept=F, lambda=lambda/shrink, 
             alpha=alpha, thresh=1e-10)
as.vector(gg$beta) - x
})

