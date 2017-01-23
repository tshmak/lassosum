#' @title Obtains the number of column (SNPs) in a PLINK bfile
#' 
#' @param bfile Plink file stem
#' @return an integer with the number of columns
#' #@keywords internal
#' @export
ncol.bfile <- function(bfile) {
	bimfile <- paste0(bfile, ".bim")
	if(!file.exists(bimfile)) 
		stop(paste0("Cannot find ", bimfile)) 
	
	return(countlines(bimfile))
	# if(Sys.info()["sysname"] == "Windows") {
	# 	wc.output <- shell(paste("wc -l", bimfile), intern=T)
	# } else {
	# 	wc.output <- system(paste("wc -l", bimfile), intern=T)
	# }
	# return(as.numeric(strsplit(wc.output, split = "\\s+")[[1]][1]))
}
