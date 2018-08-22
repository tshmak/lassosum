lassosum (standalone version for Linux)
=======================

### Description
This page is for the standalone version of `lassosum` (written in R). For full details of `lassosum`, please refer to [this page](https://github.com/tshmak/lassosum). 

### Installation
Follow the instruction [here](https://github.com/tshmak/lassosum#installation) to install `lassosum` on R. Then add the `lassosum` path to the `$PATH` variable. The `lassosum` path can be obtained by typing the following in R:  
```r
> system.file(package="lassosum")
```
For example, on my computer, I would type 
```{sh}
$ PATH=/home/tshmak/WORK/Rpackages2/nonMRAN/lassosum/:$PATH
```

### A quick example
The following is a quick example to run `lassosum` from a Linux shell, assuming lassosum has been included in `$PATH`. 
```{sh}
$ lassosum --data summarystats.txt --chr Chr --pos Position \
        --A1 A1 --A2 A2 --pval P_val --n 50000 \
        --OR OR_A1 --test.bfile testsample \
        --LDblocks EUR.hg19 --pheno testsample.pheno.txt \
        --nthreads 2
```
To actually try out the above example, copy the relevant files from the directory given by 
```r
> system.file("data", package="lassosum")
```


### Support
If there are any questions or problems with running or installing `lassosum`, please do email me at <tshmak@hku.hk>. 
