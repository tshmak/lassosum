lassosum
=======================

### Description

`lassosum` is a method for computing LASSO estimates of a linear regression problem given summary statistics from GWAS and Genome-wide meta-analyses, accounting for Linkage Disequilibrium (LD), via a reference panel.
The reference panel is assumed to be in PLINK [format](https://www.cog-genomics.org/plink2/).
Summary statistics are expected to be loaded into memory as a data.frame/data.table. 

### Installation

`lassosum` requires the following R packages: `RcppArmadillo`, `data.table`, `Matrix`. Install them by: 

```r
install.packages(c("RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
```
For Windows users, it would be easiest to download the following binary [lassosum_0.2.3.zip](https://github.com/tshmak/lassosum/releases/download/v0.2.3/lassosum_0.2.3.zip) and install using: 
```r
install.packages("/path/to/downloaded_binary_file.zip", repos=NULL)
```

For Mac and Linux users, we recommend downloading the source codes [lassosum_0.2.3.tar.gz](https://github.com/tshmak/lassosum/releases/download/v0.2.3/lassosum_0.2.3.tar.gz) and compiling on your computer. Mac users will need to install [Xcode](https://developer.apple.com/xcode/) to do this. After downloading, type:
```r
install.packages("/path/to/downloaded_source.tar.gz", repos=NULL, type="source")
```

If you have `devtools`, you can also type: 
```r
install_github("tshmak/lassosum@v0.2.3")
```
or
```r
install_github("tshmak/lassosum")
```
for the latest development version. Or you can clone the latest development version here and install yourself using `devtools`. 

### Warning!

Most functions in `lassosum` impute missing genotypes in PLINK bfiles with a homozygous A2 genotype, which is the same as using the `--fill-missing-a2` option in PLINK. It is the user's responsibility to filter out individuals and SNPs with too many missing genotypes beforehand. 

### Tutorial

In the following tutorial we make use of two dummy datasets, which can be downloaded [here](https://github.com/tshmak/lassosum/archive/v0.2.3.zip).

The data for this tutorial can be found in `tutorial/data` after unzipping. 
We will assume you have set your `R` working directory at `tutorial/` with 

```r
setwd("path/to/repository/tutorial")
```

First we read the summary statistics into R, and provide the `bfile` names of the refrence panel and the test data. If only the reference panel is provided then only the beta coefficients (no polygenic scores) are calculated. You can then apply these subsequently to a test dataset using `validate.lassosum.pipeline` or `pseudovalidate.lassosum.pipeline`. If no reference panel is provided, then the test data is taken as the reference panel. If no ld region file is provided, then `lassosum` is performed by chromosomes. We recommend you use the appropriate LD regions as defined in [Berisa and Pickrell (2015)](https://academic.oup.com/bioinformatics/article/32/2/283/1743626/Approximately-independent-linkage-disequilibrium) which are also included in our package. 

```r
library(data.table)

### Read summary statistics file ###
ss <- fread("./data/summarystats.txt")
head(ss)

### Specify the PLINK file stub of the reference panel ###
ref.bfile <- "./data/refpanel"

### Specify the PLINK file stub of the test data ###
test.bfile <- "./data/testsample"

### Read ld region file ###
ld <- fread("./data/Berisa.2015.EUR.bed")
```

To run `lassosum`, we need to input SNP-wise correlations. This can be converted from p-values via the `p2cor` function. 
```r
library(lassosum)
cor <- p2cor(p = ss$P_val, n = 60000, sign=log(ss$OR_A1))
```

Running lassosum using standard pipeline: 
```r
out <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position, 
                         A1=ss$A1, A2=ss$A2,
                         ref.bfile=ref.bfile, test.bfile=test.bfile, 
                         LDblocks = ld)

### Validation with phenotype ### 
v <- validate.lassosum.pipeline(out) # Use the 6th column in .fam file in test dataset for test phenotype

### pseudovalidation ###
# install.packages("fdrtool")
v <- pseudovalidate.lassosum.pipeline(out)

```

Parallel processing with the `parallel` package. Note that parallel processing is done by `LDblocks`. 
```r
library(parallel)
cl <- makeCluster(2, type="FORK")
out <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position, 
                         A1=ss$A1, A2=ss$A2,
                         ref.bfile=ref.bfile, test.bfile=test.bfile, 
                         LDblocks = ld, cluster=cl)
```
#### Including covariates in validation
It is possible to include covariates in validation (though not in pseudovalidation). To do so, you need to define an alternative `validate.function` to pass to `validate.lassosum.pipeline`. For example, suppose you have a matrix of covariates in `covar`. Define, e.g.: 
```r
FUN <- function(X, y) {
  res <- residuals(lm(X ~ covar))
	return(cor(res, y))
}
# covar <- rnorm(nrow.bfile(out$test.bfile)) # If you need a dummy for testing
```
Then run:
```r
v <- validate.lassosum.pipeline(out, validate.function=FUN)
```
Note that `X` and `y` are not variables in the R environment, but simply arguments in the function. Only `covar` must be defined. Moreover, it must have the same number of rows as the number of participants included in the analysis. 

##### A more technical note
`validate.lassosum.pipeline` performs `apply(PGS, MARGIN = 2, FUN=validate.function, pheno)))` to create a vector to assess the best PGS in predicting the phenotype. 


### Support
If there are any questions or problems with running or installing `lassosum`, please do email me at tshmak@hku.hk
