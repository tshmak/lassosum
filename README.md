lassosum
=======================

### Description

`lassosum` is a method for computing LASSO estimates of a linear regression problem given summary statistics from GWAS and Genome-wide meta-analyses, accounting for Linkage Disequilibrium (LD), via a reference panel.
The reference panel is assumed to be in PLINK [format](https://www.cog-genomics.org/plink2/).
Summary statistics are expected to be loaded into memory as a data.frame/data.table. 

### Installation

`lassosum` requires the following R packages: `Rcpp`, `data.table`, `Matrix`. Install them by: 

```r
install.packages(c("Rcpp", "data.table", "Matrix"))
```
Download the following file [lassosum_0.2.tar.gz](https://github.com/tshmak/lassosum/releases/download/v0.2.0/lassosum_0.2.tar.gz). Install using 
```r
install.packages("/path/to/lassosum_0.2.tar.gz", repos=NULL, type="source")
```

Or you can install the latest development version using `devtools`
```r
devtools::install_github("tshmak/lassosum")
```
### Warning!

Most functions in `lassosum` impute missing genotypes in PLINK bfiles with a homozygous A2 genotype, which is the same as using the `--fill-missing-a2` option in PLINK. It is the user's responsibility to filter out individuals and SNPs with too many missing genotypes beforehand. 

### Tutorial

We advise the use of the packages `data.table` to import summary statistics text files.

In the following tutorial we make use of two dummy datasets, which can be downloaded from this repository.
You can download the repository via

```bash
git clone https://github.com/tshmak/lassosum
```

or just download the [ZIP file](https://github.com/tshmak/lassosum/archive/v0.2.0.zip).

The data for this tutorial is stored in `tutorial/data`. 
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
### Support
If there are any questions, please email Timothy Mak at tshmak@hku.hk
