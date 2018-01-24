Cross prediction using lassosum 
=======================
  
Cross prediction is a method for calculating polygenic scores in large cohorts without the use of summary statistics. For details, refer to this [paper](https://www.biorxiv.org/content/early/2018/01/23/252270). 

# Tutorial 

If **lassosum** is not yet installed, refer to the instruction [here](https://github.com/tshmak/lassosum#installation) for installation. If **lassosum** is not already loaded, load it via: 
```r
library(lassosum)
```
**lassosum** uses [Plink](https://www.cog-genomics.org/plink2/) to calculate summary statistics. Before the functions that call plink can be used, the link to the plink executable needs to be specified by: 
```r
options(lassosum.plink='/path/to/plink')
```

We assume that we have genotype in PLINK 1 [format](https://www.cog-genomics.org/plink/1.9/input#bed). For example, let's say our data files are: `mydata.bed`, `mydata.bim`, and `mydata.fam`. 

We first carry out summary statistics estimation for the different folds. This can be easily done by
```r
sumstats <- cp.plink.linear(bfile="mydata")
```
The default is 5-fold cross-prediction. The folds are randomly assigned to the samples. Use `nfolds` to specify the number of folds needed, or `fold` to allocate the fold yourself. Type `help(cp.plink.linear)` for more details. 

`sumstats` can then be fed into `cp.lassosum`. 
```r 
ld <- read.table(system.file("data/Berisa.EUR.hg19.bed", package="lassosum"), header=T)
cp <- cp.lassosum(sumstats, LDblocks=ld)
```
We recommend you use one of the LD blocks given in the package, if you do not have your own LD blocks defined. These LD blocks are given by the paper [Berisa and Pickrell (2015)](https://academic.oup.com/bioinformatics/article/32/2/283/1743626/Approximately-independent-linkage-disequilibrium), and are based on the 1000 Genome data. Replace EUR with ASN or AFR for Asian or African LD regions, respectively. hg38 coordinates are also available by [liftOver](https://genome.sph.umich.edu/wiki/LiftOver). Simply replace `hg19` with `hg38` above. 

By default, if the sample size is > 5000, `cp.lassosum` uses a random subset of 5000 as the reference panel. This is to reduce the computation burden. Increase or decrease this using the `max.ref.bfile.n` option. Alternatively, specify the exact sample to use using the `keep.ref` option. 

For further options, please see the manual
```r
help(cp.lassosum)
```
or email me at <tshmak@hku.hk>. 

### Running cross-prediction using separated .bed files
In large datasets, data is often stored across chromosomes in separate `.bed` files. To run cross prediction across different chromosomes, simply run `cp.plink.linear` separately for the different `.bed` files, then run `cp.lassosum` separately for each output from `cp.plink.linear`, specifying `list.of.lpipe.output=TRUE`. This will generate a list of `lassosum.pipeline` objects for each fold. However, **remember to set the random number seed to the same number before you run `cp.plink.linear` so that the folds are defined consistently, and also before `cp.lassosum` if your sample size is > `max.ref.bfile.n`, so that the same reference sample is chosen!** Then, use `organize.by.fold` to reorganize the all the outputs for the different chromosomes. Finally, run `cp.lassosum` with the `list.of.lpipe.input` option. Below is an example for two chromosomes. 
```r
set.seed(42)
sumstats1 <- cp.plink.linear("chr1")
set.seed(42)
sumstats2 <- cp.plink.linear("chr2")
set.seed(1000)
lp1 <- cp.lassosum(sumstats1, LDblocks=ld, list.of.lpipe.output=TRUE)
set.seed(1000)
lp2 <- cp.lassosum(sumstats2, LDblocks=ld, list.of.lpipe.output=TRUE)
lp <- organize.by.fold(list(lp1, lp2))
cp <- cp.lassosum(list.of.lpipe.input=lp)
```
### Running cross-validation using results from cross-prediction
Method 1 in cross-prediction essentially identifies the best `lambda` and `s` for use with **lassosum**. We can then apply this to the entire dataset. This is the cross-validation procedure in our paper, which can be achieved by running: 
```r
cv <- cp.cv(sumstats, cp)
```
`sumstats` can be a `list` of `cp.plink.linear` objects if they were obtained separately for different chromosomes. Note that the order of the list should match that given to `cp.lassosum` earlier. Note that `cp.cv` does not calculate summary statistics afresh using the entire data, but simply averages those obtained across the different folds in `cp.lassosum`. This makes it very fast. 

### Incorporating external summary statistics into cross-prediction 
Suppose external summary statistics are loaded into a `data.frame` called `dat` with the columns `chr`, `pos`, `alt`, `rsid`, `beta`, `p` for chromosomes, position, alternative allele, rsID, beta coefficients, and p-values respectively. We first need to convert the p-values to correlation coefficients using the `p2cor` command: 
```r
dat$cor <- p2cor(dat$p, n = 60000, sign=dat$beta)
```
where `n = 60000` is the sample size. Assuming `cp.pl` is an object from `cp.plink.linear`, we can then run 
```r
merged <- cp.meta(cp.pl, chr=chr, pos=pos, snp=rsid, A1=alt, cor=cor, n=60000)
```
to merge it with `cp.pl`. We can then proceed with `merged` as if it were a `cp.plink.linear` object in running, e.g., 
```r 
cp <- cp.lassosum(merged, LDblocks=ld)
```
Note that not all of the variables `chr`, `pos`, `snp`, `A1`, `A2` need to be specified. At least one of `A1` and `A2` must be specified. Otherwise specify as many of `chr`, `pos`, and `snp` as you need. 

### Using results from PRSice in cross-prediction 
I will write a function for this if there is demand. 

Please email me <tshmak@hku.hk> for any bug reports, comments, or suggestions. Thanks for using **lassosum**!



