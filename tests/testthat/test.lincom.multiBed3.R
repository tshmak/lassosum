test_that('pgs' , {
n <- nrow.bfile("test")
p <- ncol.bfile("test")
X <- genotypeMatrix("test.bed", n, p, 
                    integer(0), integer(0), integer(0), integer(0), 1)
expect_equal(sum(X), 27)
j <- rnorm(p)
sum(pgs("test",j))
beta.coef <- matrix(rnorm(p*3), p, 3)
system.time(a <- pgs("test", beta.coef))
})

