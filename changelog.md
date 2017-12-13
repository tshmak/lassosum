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

