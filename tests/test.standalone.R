devel()
exec <- system.file("inst/lassosum.R", package="lassosum")
system(paste("chmod +x", exec))
# system(exec)
setwd(system.file("tests", package="lassosum"))
f <- function(file="") paste0(system.file("data", package="lassosum"), "/", file)
ss <- read.table2(f("summarystats.txt"), header=TRUE)
options <- list(data=f("summarystats.txt"), chr="Chr", pos="Position", 
     A1="A1", A2="A2", pval="P_val", n=50000, OR="OR_A1", 
     test.bfile=f("testsample"), LDblocks="EUR.hg19", 
     validate=T, splitvalidate=T, pseudovalidate=T)
cmd <- parse.options(options = options, dot2dash=F)
# print(paste(exec, cmd))f
# source(exec)
system(paste(exec, cmd))
head(read.table2("lassosum.validate.results.txt", header=T))

Tim.load(Rplink)
fam <- read.fam(f("testsample"))
fam2 <- fam[sample(nrow(fam), 50),]
write.table2(fam2[,1:2], file=keep <- tempfile())
options$keep.test <- keep
system(paste(exec, parse.options(options = options, dot2dash=F)))
lp <- readRDS("lassosum.lassosum.pipeline.rds")
