# Tmisc()
# setwd0("~/WORK/myRpackages/lassosum/inst/data/")
# for(pop in c("EUR", "ASN", "AFR")) {
#   # pop <- "EUR"
#   bed <- read.table.tim(paste("Berisa", pop, "hg19", "bed", sep="."))
#   l <- liftover(bed$chr, bed$start, bed$stop, from=19, to=38, bedformat = TRUE)
#   chrs <- unique(l$chr)
#   ll <- lapply(chrs, function(x) subset(l, chr==x))
#   lll <- list()
#   for(i in 1:length(ll)) {
#     # i <- 1
#     u <- unique(c(ll[[i]]$pos, ll[[i]]$pos2))
#     u <- sort(u[!is.na(u)])
#     lll[[i]] <- data.frame(chr=chrs[[i]], start=u[-length(u)], stop=u[-1])
#   }
#   new <- do.call(rbind, lll)
#   write.table.tim(new, file=paste("Berisa", pop, "hg38", "bed", sep="."))
#   
# }
