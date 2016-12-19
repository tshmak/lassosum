# load_all(); setwd("./tutorial")
### Use data.table for speed ###
library(data.table)

### Locate reference panel .bfile (the .bed, .bim and .fam files without the extensions) ###
ref.bfile <-"./data/refpanel"

### Locate test sample .bfile ###
test.bfile <-"./data/testsample" 

### Read ld region file ###
ld <- fread("./data/Berisa.2015.EUR.bed")

### Read summary statistics file ###
ss <- fread("./data/summarystats.txt")
head(ss)

### Convert p-values to correlations, assuming a sample size of 60000 for the p-values ###
cor <- p2cor(p = ss$P_val, n = 60000, sign=log(ss$OR_A1))

### Run lassosum using standard pipeline ### 
out <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position, 
                         A1=ss$A1, A2=ss$A2,
                         ref.bfile=ref.bfile, test.bfile=test.bfile, 
                         LDblocks = ld)

### Validation with phenotype ### 
v <- validate.lassosum.pipeline(out) # Use the 6th column in .fam file in test dataset

### pseudovalidation ###
# install.packages("fdrtool")
v <- pseudovalidate.lassosum.pipeline(out)

