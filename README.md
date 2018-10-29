lassosum [![Build Status](https://travis-ci.org/tshmak/lassosum.svg?branch=master)](https://travis-ci.org/tshmak/lassosum)
=======================
*New!! A standalone version of `lassosum` is now available. Please see [here](https://github.com/tshmak/lassosum/blob/master/lassosum_standalone.md#lassosum-standalone-version-for-linux) for details.*

### Description

`lassosum` is a method for computing LASSO/Elastic Net estimates of a linear regression problem given summary statistics from GWAS and Genome-wide meta-analyses, accounting for Linkage Disequilibrium (LD), via a reference panel.
The reference panel is assumed to be in PLINK 1 [format](https://www.cog-genomics.org/plink/1.9/input#bed).
Summary statistics are expected to be loaded into memory as a data.frame/data.table. 

### Installation

`lassosum` requires the following R packages: `RcppArmadillo`, `data.table`, `Matrix`. Install them by: 

```r
install.packages(c("RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
```
For Windows users, it would be easiest to download the following binary [lassosum_0.4.3.zip](https://github.com/tshmak/lassosum/releases/download/v0.4.3/lassosum_0.4.3.zip) and install using: 
```r
install.packages("/path/to/downloaded_binary_file.zip", repos=NULL)
```

For Mac and Linux users, we recommend downloading the source codes [lassosum_0.4.3.tar.gz](https://github.com/tshmak/lassosum/releases/download/v0.4.3/lassosum_0.4.3.tar.gz) and compiling on your computer. Mac users will need to install [Xcode](https://developer.apple.com/xcode/) to do this. After downloading, type:
```r
install.packages("/path/to/downloaded_source.tar.gz", repos=NULL, type="source")
```

If you have `devtools`, you can also type: 
```r
install_github("tshmak/lassosum@v0.4.3")
```
or
```r
install_github("tshmak/lassosum")
```
for the latest development version. Or you can clone the latest development version here and install yourself using `devtools`. 

### Warning!

Most functions in `lassosum` impute missing genotypes in PLINK bfiles with a homozygous A2 genotype, which is the same as using the `--fill-missing-a2` option in PLINK. It is the user's responsibility to filter out individuals and SNPs with too many missing genotypes beforehand. 

### Tutorial

_If you are using lassosum for __cross prediction__, please refer to the manual [here](https://github.com/tshmak/crosspred)_

Otherwise, run the following: 
```r
library(lassosum)
setwd(system.file("data", package="lassosum")) # Directory where data and LD region files are stored
```

First we read the summary statistics into R, and provide the `bfile` names of the reference panel and the test data. If only the reference panel is provided then only the beta coefficients (no polygenic scores) are calculated. You can then apply these subsequently to a test dataset using `validate` or `pseudovalidate`. If no reference panel is provided, then the test data is taken as the reference panel. If no ld region file is provided, then `lassosum` is performed by chromosomes. We recommend you use the appropriate LD regions as defined in [Berisa and Pickrell (2015)](https://academic.oup.com/bioinformatics/article/32/2/283/1743626/Approximately-independent-linkage-disequilibrium) which are also included in our package. 

```r
library(data.table)

### Read summary statistics file ###
ss <- fread("summarystats.txt")
head(ss)

### Specify the PLINK file stub of the reference panel ###
ref.bfile <- "refpanel"

### Specify the PLINK file stub of the test data ###
test.bfile <- "testsample"

### Read ld region file ###
ld <- fread("Berisa.EUR.hg19.bed") # Replace EUR with ASN or AFR for Asian or African. Replace hg19 with hg38 for hg38 coordinates. 
```

To run `lassosum`, we need to input SNP-wise correlations. This can be converted from p-values via the `p2cor` function. 
```r
library(lassosum)
cor <- p2cor(p = ss$P_val, n = 60000, sign=log(ss$OR_A1))
# n is the sample size
```

Running lassosum using standard pipeline: 
```r
out <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position, 
                         A1=ss$A1, A2=ss$A2,
                         ref.bfile=ref.bfile, test.bfile=test.bfile, 
                         LDblocks = ld)

### Validation with phenotype ### 
v <- validate(out) # Use the 6th column in .fam file in test dataset for test phenotype
v <- validate(out, pheno=pheno) # Alternatively, specify the phenotype in the argument

# pheno <- rnorm(nrow.bfile(out$test.bfile)) # If you need a dummy for testing

### pseudovalidation ###
# install.packages("fdrtool")
v <- pseudovalidate(out)
```
Since v0.4.2, the `pheno` argument in `validate` can also take a `data.frame` with the first 2 columns headed by FID and IID, and the third column being the phenotype, or alternatively, a file name for such a data.frame. Moreover, a `v$results.table` object is also returned in `validate` and `pseudovalidate`, giving a table with the best PGS and the phenotype tabulated with the FID and IID (family and individual ID). 

A new feature since v0.4.2 is `split-validation`, where the test dataset (`test.bfile`) is split in half using one half for validation and the other half for calculating PGS. The PGS in the two halves are then standardised and stacked back together. This avoids overfitting due to the overlapping of the target and the validation dataset. See [this paper](https://www.biorxiv.org/content/early/2018/07/30/252270) for details. 
```r
### Split-validation ###
sv <- splitvalidate(out)
```

#### Parallel processing with the `parallel` package
Note that parallel processing is done by `LDblocks`. 
```r
library(parallel)
cl <- makeCluster(2, type="FORK")
out <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position, 
                         A1=ss$A1, A2=ss$A2,
                         ref.bfile=ref.bfile, test.bfile=test.bfile, 
                         LDblocks = ld, cluster=cl)
```
#### Including covariates in validation
It is possible to include covariates in validation or splitvalidation (though not in pseudovalidation). Simply pass the covariate matrix as an argument to `validate` or `splitvalidate`. 
```r 
v <- validate(out, covar=covar)
# covar <- rnorm(nrow.bfile(out$test.bfile)) # If you need a dummy for testing
```
Since v0.4.2, the `covar` argument in `validate` and `splitvalidate` can also take a `data.frame` with the first 2 columns headed by FID and IID, and the other columns being covariates (any headers). It can also be a file name for such a data.frame.

#### Apply validated betas to new data 
To apply the best lassosum predictor (indexed by `s` and `lambda`) to a new dataset, first subset the `lassosum.pipeline` object. Then `validate` again: 
```r 
out2 <- subset(out, s=v$best.s, lambda=v$best.lambda)
v2 <- validate(out2, covar=covar, test.bfile="Some_new_bfile")
```

### Support
If there are any questions or problems with running or installing `lassosum`, please do email me at <tshmak@hku.hk>. 
