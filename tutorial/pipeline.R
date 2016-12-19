# load_all(); setwd("./tutorial")
### Use data.table for speed ###
library(data.table)

### Locate reference panel .bfile (the .bed, .bim and .fam files without the extensions) ###
ref.bfile <-"./data/refpanel"

### Locate test sample .bfile ###
test.bfile <-"./data/testsample" 

### Read ld region file ###
ld <- fread("./data/fourier_ls-all.bed", data.table=F)

### Read summary statistics file ###
ss <- fread("./data/summarystats.txt", data.table=F)
head(ss)

### Run lassosum using standard pipeline ### 
out <- lassosum.pipeline(ref.bfile=ref.bfile, 
                         test.bfile=test.bfile, 
                         pvals=ss$P_val, 
                         sample.size = 60000, 
                         coef.sign = log(ss$OR_A1), 
                         chr = ss$Chr, 
                         pos=ss$Position, 
                         A1 = ss$A1, A2=ss$A2,
                         LDblocks = ld)
                         
### Validation with phenotype ### 
fake.pheno <- rnorm(nrow.bfile(test.bfile))
v <- validate.lassosum.pipeline(out, pheno=fake.pheno)

### pseudovalidation ###
# install.packages("fdrtool")
v <- pseudovalidate.lassosum.pipeline(out)

