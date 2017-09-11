plink.linear <- function(bfile, pheno, out=tempfile("lassosum.out"), 
                         keep=NULL, remove=NULL, 
                         extract=NULL, exclude=NULL, 
                         chr=NULL, 
                         covar=NULL, 
                         plink.cmd="--linear standard-beta", 
                         ext=".assoc.linear",
                         show.output.on.console=FALSE, 
                         ...) {
  
  #' @title Obtain standardized coefficients from linear regression in PLINK
  #' @param keep,remove,extract,exclude,chr see parseselect() 
  
  #### checks ####
  plink <- getOption("lassosum.plink")
  if(is.null(plink)) {
    stop(paste("plink executive for lassosum not yet specified.",
               "Please specify by typing", 
               "options(lassosum.plink='/path/to/plink')"))
  } 
  
  #### plink options ####
  options <- list(...)
  options$keep.allele.order <- ""
  options$allow.no.sex <- ""
  options$bfile <- bfile
  
  #### parse ####
  parsed <- parseselect(bfile=bfile, extract=extract, exclude=exclude, 
                      keep=keep, remove=remove, chr=chr, export=TRUE)
  
  #### keep ####
  if(!is.null(parsed$keep)) {
    if(is.null(parsed$fam)) parsed$fam <- read.table2(parsed$famfile) 
    tokeep <- tempfile("lassosum")
    write.table(parsed$fam[parsed$keep,], file=tokeep, 
                quote=FALSE, row.names = FALSE, col.names=FALSE)
    options$keep <- tokeep
  }
  
  #### extract ####
  if(!is.null(parsed$extract)) {
    if(is.null(parsed$bim)) parsed$bim <- read.table2(parsed$bimfile) 
    toextract <- tempfile("lassosum")
    write.table(parsed$bim[parsed$extract,], file=toextract, 
                quote=FALSE, row.names = FALSE, col.names=FALSE)
    options$extract <- toextract
  }
  
  #### pheno ####
  if(!is.null(pheno)) {
    options$pheno <- parse.pheno.covar(pheno, parsed)
    options$no.pheno <- ""
  }

  #### covar ####
  if(!is.null(covar)) options$covar <- parse.pheno.covar(covar, parsed)
  
  #### out ####
  options$out <- out
  
  #### Form plink command ####
  parse.plink.options <- function(options) {
    l <- length(options)
    names <- names(options)
    Names <- gsub("\\.", "-", names)
    cmd <- ""
    for(i in 1:l) {
      if(is.null(options[[i]])) options[[i]] <- ""
      cmd <- paste(cmd, paste0("--", Names[i]), options[[i]])
    }
    return(cmd)
  }
  cmd <- parse.plink.options(options)
    
  #### run plink ####
  cmd <- paste(plink, cmd, plink.cmd)
  system(cmd, show.output.on.console=show.output.on.console)
  
  tab <- read.table2(paste0(out, ext), header=TRUE)
  attr(tab, "out") <- options$out

  return(tab)

}


#### For me only ####
if(exists("attachroot")) {
  if(Sys.info()["sysname"] == "Windows") {
    options(lassosum.plink="D:/PLINK/plink.exe")
  } else {
    options(lassosum.plink="/home/tshmak/software/plink/v1.90b3.44/plink")
  }
}
