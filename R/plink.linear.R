plink.linear <- function(bfile, pheno, out=tempfile("lassosum.out"), 
                         keep=NULL, remove=NULL, 
                         extract=NULL, exclude=NULL, 
                         chr=NULL, 
                         covar=NULL, 
                         plink.cmd="--linear standard-beta", 
                         ext=".assoc.linear",
                         ignore.stdout=TRUE, 
                         ...) {
  
  #' @title Obtain standardized coefficients from linear regression in PLINK
  #' @param pheno Phenotype. see details.
  #' @param out The \code{--out} option in plink
  #' @param keep,remove,extract,exclude,chr see parseselect() 
  #' @param covar Covariates. see details
  #' @param plink.cmd The command passed to plink
  #' @param ext The extension for the output file generated from the plink command
  #' @param ignore.stdout Ignore stdout. see \code{\link[base]{base::system}}
  #' @param ... Other parameters to pass to plink. See details
  #' @export
  #' @details \code{pheno} and \code{covar} can take one of three formats:
  #'          \itemize{
  #'          \item a vector/matrix with length/number of rows equal the 
  #'          number of included samples, or
  #'          \item a data.frame with the 
  #'          first two columns labelled "FID" and "IID" and the other columns
  #'          giving the phenotype/covariates, or
  #'          \item the name of a file in a format 
  #'          understood by plink's \code{--pheno} or \code{--covar} options
  #'          } 
  #'          Only one column can be given for \code{pheno}, although covar can take 
  #'          more than one columns. 
  #'          
  #'          \code{...}: Other options to plink can be given by, e.g., 
  #'          \code{keep.allele.order=T} for \code{--keep-allele-order}, 
  #'          \code{maf=0.01} for \code{--maf 0.01}, etc.
  #'           
  
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
  system(cmd, ignore.stdout=ignore.stdout)
  
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
