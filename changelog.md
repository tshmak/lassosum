## V0.2.0 
* Initial release 

## v0.2.1
* Use a c++ function for counting lines instead of bash `wc`
* Corrects bug in `lassosum.pipeline()` when calling some functions in `data.table`

## v0.2.2
* includes some additional objects in the output of `validate.lassosum.pipeline()` and `pseudovalidate.lassosum.pipeline()`
* Allows the `validate.function` in `validate.lassosum.pipeline()` to be given as a character rather than a function. 

## v0.2.3
* `best.beta` is output in `validate.lassosum.pipeline()` and `pseudovalidate.lassosum.pipeline()`
* fixes bug when grouping blocks that are too small

## v0.2.4
* fixes bug when matching variable names in `matchpos()`
* Support for parallel processing in the `pgs()` function, (and therefore `validate.lassosum.pipeline()` and `pseudovalidate.lassosum.pipeline()`)
* Allows summary statistics to be matched by SNP id as well as chromosome/position 
* fixes bug when the `bfile` is too large (due to integer overflow in c++)
* adds a user-defined limit to the size of the bfile to discourage the use of large reference panels
* cpp functions are now exported (except `openPlinkBinaryFile()`)
* includes this changelog file
* fixes bug in plotting results from `validate.lassosum.pipeline()` and `pseudovalidate.lassosum.pipeline()`
* fixes bug in displaying "Running lassosum with s=1..." message. 
* renames hg19 LD region files
* adds hg38 LD region files (lifted over from the hg19 ones)

## v0.2.5
* Bug fix for `lassosum.pipeline()` in parsing LDblocks when given as a vector of factors
* Bug fix for `pgs()` when `extract` is specified 

## v0.3.0 
* Adds function for cross-prediction [link]()
* Allows multiple bfiles in many functions 

## v0.3.1
* Adds `cp.meta()` for meta-analysing external summary statistics with results from `cp.plink.linear()`
* `cp.plink.linear()` will weight correlation by sample size when averaging across folds
* A number of unrelated functions are now hidden. 
* Improved documentation for `cp` functions

## v0.3.2 
* fixes bug in cp.plink.linear covar option
* fixes bug when specifying keep
* adds `write.table2()` 

## v0.4.0 
* Removes support for cross-prediction. 
* Cross-prediction moved to a separate package
* recast `validate.lassosum.pipeline()`, `pseudovalidate.lassosum.pipeline()`, and `pgs` as S3 methods
* Adds option for including covariates in `validate.lassosum.pipelin()`
* bug fix in `matchpos()`
* Introduces c++ function for `pgs()` with sparse matrix. Used as default. 
* `p2cor()` will automatically code to NA those correlations with relatively small n.
* Always uses `data.table::fread()` for reading files

## v0.4.1
* Bug fix for `matchpos()`
* `pgs()` stops if keep or extract not in the same order as in .bim/.fam file

## v0.4.2
* Allows `covar` and `pheno` to take `data.frame`s or files as input. 
* You no longer have to re-enter `test.bfile` if `keep` or `remove` is used in `validate.lassosum.pipeline()` or `pseudovalidate.lassosum.pipeline()` 
* Adds `splitvalidate()` for split-validation. 
* Adds `subset.lassosum.pipeline()` for applying validated PGS to new data. 
* Allows `pheno` to be constant or NA if only one PGS is calculated in `validate.lassosum.pipeline()`. 
* `plot.validate.lassosum()` will abandon if there are no finite values. 
* **A standalone version of lassosum is available** 

## v0.4.3
* LDblocks does not need to be specified when s < 1
* Better tracing in SD calculation
* Bug fix in `validate.lassosum.pipeline()`
* Allows multiple `test.bfile` in `validate.lassosum.pipeline()`

## v0.4.4
* Bug fix
* Better documentation for Mac users

## v0.4.5
* [Bug fix] standalone version used to ignore --nthreads option when running `validate.lassosum.pipeline`
* [Bug fix] The `pheno` and `covar` option didn't use to like `data.table` and required `data.frame`
* [Bug fix] Pull request #20
* Adds seed option in standalone version 

