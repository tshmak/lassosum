#' @title Obtains the number of individuals in a PLINK bfile
#' 
#' @param bfile Plink file stem
#' @return an integer with the number of rows
#' #@keywords internal
#' @export
nrow.bfile <- function(bfile) {
	famfile <- paste0(bfile, ".fam")
	if(!file.exists(famfile)) 
		stop(paste0("Cannot find ", famfile)) 
	return(countlines(famfile))
	# if(Sys.info()["sysname"] == "Windows") {
	# 	wc.output <- shell(paste("wc -l", famfile), intern=T)
	# } else {
	# 	wc.output <- system(paste("wc -l", famfile), intern=T)
	# }
	# return(as.numeric(strsplit(wc.output, split = "\\s+")[[1]][1]))
}
