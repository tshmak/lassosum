clear()
setwd(attachroot("/WORK/2016-12/lassosum/tutorial/data"))
load_all(attachroot("/WORK/2016-12/lassosum/."))
Tim.load(Rplink)
Tim.load(SNPdata)
chr22 <- "/psychipc01/disk2/references/1000Genome/release/20130502_v5a/PLINK/chr22"
chr21 <- "/psychipc01/disk2/references/1000Genome/release/20130502_v5a/PLINK/chr21"
subsample(bfile=chr21, out="chr21", maf=0.05, fill.missing.a2=T, 
          keep=logical.vector(1:200,2504), p=1000, P=20000, silent=F)
subsample(bfile=chr22, out="chr22", maf=0.05, fill.missing.a2=T, 
          keep=logical.vector(1:200,2504), p=1000, P=20000, silent=F)
plink(bfile="chr21", cmd="--bmerge chr22", out="test", silent=F)

subsample(bfile="test", QC=F, p=1800, out="testsample")
subsample(bfile="test", QC=F, p=1800, out="refpanel")
test.fam <- read.table.tim("testsample.fam")
test.fam$V6 <- rbinom(nrow(test.fam), 1, 0.5)
write.table.tim(test.fam, file="testsample.fam", col.names=F)
library(data.table)
sumstats <- fread(attachroot("/DATA/summarystats/RA/RA_GWASmeta_European_v2.txt"))
bim <- fread("test.bim")
m <- matchpos(sumstats,bim, ref.chr = "V1", ref.pos="V4", rm.duplicates = T)
sumstats.tokeep <- sumstats[m$order]
colnames(sumstats.tokeep)[c(3,6,7,8,9)] <- c("Position", "OR_A1", "OR_95_CIlow", "OR_95_CIup", "P_val")
write.table.tim(sumstats.tokeep, file="summarystats.txt")
system("cp ~/DATA/ldetect/src/EUR/fourier_ls-all.bed .")
system("rm chr21.* chr22.* test.* *.log *.nosex")
