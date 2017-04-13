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
For Windows users, it would be easiest to download the following binary [lassosum_0.2.zip](https://github.com/tshmak/lassosum/releases/download/v0.2.0/lassosum_0.2.zip) and install using: 
```r
install.packages("/path/to/downloaded_binary_file.zip", repos=NULL)
```

For Mac and Linux users, we recommend downloading the source codes [lassosum_0.2.tar.gz](https://github.com/tshmak/lassosum/releases/download/v0.2.0/lassosum_0.2.tar.gz) and compiling on your platform by:
```r
install.packages("/path/to/downloaded_source.tar.gz", repos=NULL, type="source")
```

Or you can clone the latest development version here and install yourself using `devtools`. 

### Warning!

Most functions in `lassosum` impute missing genotypes in PLINK bfiles with a homozygous A2 genotype, which is the same as using the `--fill-missing-a2` option in PLINK. It is the user's responsibility to filter out individuals and SNPs with too many missing genotypes beforehand. 

### Tutorial

In the following tutorial we make use of two dummy datasets, which can be downloaded [here](https://github.com/tshmak/lassosum/archive/v0.2.0.zip).

The data for this tutorial can be found in `tutorial/data` after unzipping. 
We will assume you have set your `R` working directory at `tutorial/` with 

```r
setwd("path/to/repository/tutorial")
```

First we read the summary statistics and genotyoe information of the refrence panel into R. (`read.table` is ok, but `fread` from the `data.table` package is much faster for large files.)



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

Parallel processing with the `parallel` package. Note that parallel processing is done by `LDblocks` so always define LDblocks when running parallel. 
```r
library(parallel)
cl <- makeCluster(2, type="FORK")
out <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position, 
                         A1=ss$A1, A2=ss$A2,
                         ref.bfile=ref.bfile, test.bfile=test.bfile, 
                         LDblocks = ld, cluster=cl)
```
### Support
If there are any questions or problems with running or installing `lassosum`, please do email me at tshmak@hku.hk
