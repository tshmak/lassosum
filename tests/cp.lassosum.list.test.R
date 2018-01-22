# clear() # Takes time to load so don't restart often
# setwd(attachroot("/WORK/Projects/validation/simulations/BBsim/data/sim_1000/sumstats"))
# load("ss_0.1_0.5_2000_.1.RData")
# 
# load_all(attachroot("/WORK/Projects/validation/lassosum/."))
# ss.list <- result[[1]][[2]][1:2]
# 
# test <- cp.lassosum(ss.list, trace=3, s=1, Method2=TRUE)
# stop()
# 
# trace <- 3
# l <- list()
# for(i in 1:length(ss.list)) {
#   if(trace > 0) cat("Processing list item", i, "of", length(ss.list), "\n")
#   l[[i]] <- cp.lassosum(ss.list[[i]], trace=trace-1, list.of.lassosum.only=TRUE, s=1)
# }
# save(l, file="temp")
# 
# # Organize by folds 
# cp.l <- l[[1]]
# nfolds <- length(cp.l)
# ll <- l.folds <- list()
# for(f in 1:nfolds) {
#   l.folds[[f]] <- lapply(l, function(x) x[[f]])
#   ll[[f]] <- do.call("merge.lassosum.pipeline", l.folds[[f]])
# }
# 
# 
# 
