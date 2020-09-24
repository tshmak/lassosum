#' @title Parse the keep/remove/extract/exclude/chr options
#' @details \code{keep} and \code{remove} can take one of three
#' formats:
#' 1. A logical vector indicating which indivduals to keep/remove,
#'
#' 2. A \code{data.frame} with two columns giving the FID and IID of the indivdiuals to keep/remove (matching those in the .fam file), or
#' 3. A character scalar giving the text file with the FID/IID.
#'
#' Note that these files should not contain headers.
#'
#' \code{extract} and \code{exclude} can be of format type __(1)__ or __(3)__ describe above, or a character vector of SNP ids (matching those in the .bim file).

#' @param bfile plink file stem
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
#' @param export Logical. Whether to include the \code{bim} and/or \code{fam} data in the returned object.
#' @param order.important Logical. Whether to strictly enforce that the order of \code{extract} and/or \code{keep} must match with \code{bim} or \code{fam}, respectively.
#'
#' @md
#' @export
parseselect <- function(bfile, extract=NULL, exclude=NULL,
                        keep=NULL, remove=NULL, chr=NULL,
                        export=FALSE, order.important=FALSE) {

  #### Introduced for parsing multiple bfiles ####
  if(length(bfile) > 1) {
    if(!is.null(chr) || !is.null(extract) || !is.null(exclude)) {
      stop("bfile cannot be a vector if we are parsing extract/exclude/chr. yet.")
    }
    bfile <- bfile[1]
  }

  #### Checks ####
  stopifnot(is.character(bfile) && length(bfile) == 1)
  bedfile <- paste0(bfile, ".bed")
  bimfile <- paste0(bfile, ".bim")
  famfile <- paste0(bfile, ".fam")
  stopifnot(file.exists(bedfile))
  stopifnot(file.exists(bimfile))
  stopifnot(file.exists(famfile))

  if(grepl("^~", bfile)) {
    stop("Don't use '~' as a shortcut for the home directory.")
  }

  p <- P <- ncol.bfile(bfile)
  n <- N <- nrow.bfile(bfile)
  bim <- NULL
  fam <- NULL

  #### extract ####
  if(!is.null(extract)) {
    if(is.logical(extract)) {
      stopifnot(length(extract) == P)
    } else {
      if(is.null(attr(extract, "not.a.file", exact=T))) {
        if(is.character(extract) && length(extract) == 1 && file.exists(extract)) {
          ### I'm interpreting this as a filename
          SNPs <- read.table2(extract)
          stopifnot(ncol(SNPs)==1)
          extract <- SNPs[[1]]
        }
      } else {
        attr(extract, "not.a.file") <- NULL
      }
      if(is.vector(extract)) {
        Extract <- as.character(extract)
        bim <- read.table2(bimfile)
        extract <- bim$V2 %in% Extract
        if(order.important) {
          if(!all(bim$V2[extract] == Extract)) {
            stop("Order of extract SNPs does not match that in .bim file.")
          }
        }
      } else {
        stop("I don't know what to do with this type of input for extract")
      }
    }

    p <- sum(extract)
  }

  #### exclude ####
  if(!is.null(exclude)) {
    if(is.logical(exclude)) {
      stopifnot(length(exclude) == P)
    } else {
      if(is.null(attr(exclude, "not.a.file", exact=T))) {
        if(is.character(exclude) && length(exclude) == 1 && file.exists(exclude)) {
          ### I'm interpreting this as a filename
          SNPs <- read.table2(exclude)
          stopifnot(ncol(SNPs)==1)
          exclude <- SNPs[[1]]
        }
      } else {
        attr(exclude, "not.a.file") <- NULL
      }
      if(is.vector(exclude)) {
        exclude <- as.character(exclude)
        if(is.null(bim)) bim <- read.table2(bimfile)
        exclude <- bim$V2 %in% exclude
      } else {
        stop("I don't know what to do with this type of input for exclude")
      }

    }

    if(is.null(extract)) extract <- !exclude else
      extract <- extract & !exclude

    p <- sum(extract)
  }

  #### chr ####
  if(!is.null(chr)) {

    stopifnot(is.vector(chr))
    chr <- as.character(chr)

    if(is.null(bim)) bim <- read.table2(bimfile)
    bimchr <- bim$V1
    bimchr[bimchr==""]
    extract.chr <- bim$V1 %in% chr

    if(is.null(extract)) extract <- extract.chr else
      extract <- extract & extract.chr

    p <- sum(extract)

  }

  #### keep ####
  if(!is.null(keep)) {
    if(is.logical(keep)) {
      stopifnot(length(keep) == N)
    } else {
      if(is.character(keep) && length(keep) == 1 && file.exists(keep)) {
        ### I'm interpreting this as a filename
        keep <- read.table2(keep)
      }
      if(inherits(keep, "data.frame")) {
        keep <- as.data.frame(keep)
        stopifnot(ncol(keep)==2)
        fam <- read.table2(famfile)
        famID <- paste(fam[,1], fam[,2], sep=".")
        keepID <- paste(keep[,1], keep[,2], sep=".")
        keep <- famID %in% keepID
        if(order.important) {
          if(!all(famID[keep] == keepID)) {
            stop("Order of keep doesn't match that in .fam file")
          }
        }
      } else {
        stop("I don't know what to do with this type of input for keep")
      }

    }
    n <- sum(keep)
  }

  #### remove ####
  if(!is.null(remove)) {
    if(is.logical(remove)) {
      stopifnot(length(remove) == N)
    } else {
      if(is.character(remove) && length(remove) == 1 && file.exists(remove)) {
        ### I'm interpreting this as a filename
        remove <- read.table2(remove)
      }
      if(inherits(remove, "data.frame")) {
        remove <- as.data.frame(remove)
        stopifnot(ncol(remove)==2)
        if(is.null(fam)) fam <- read.table2(famfile)
        famID <- paste(fam[,1], fam[,2], sep=".")
        removeID <- paste(remove[,1], remove[,2], sep=".")
        remove <- famID %in% removeID
      } else {
        stop("I don't know what to do with this type of input for remove")
      }
    }

    if(is.null(keep)) keep <- !remove else
      keep <- keep & !remove

    n <- sum(keep)
  }

  if(n==0) stop("No individuals left after keep/remove! Make sure the FID/IID are correct.")
  if(p==0) stop("No SNPs left after extract/exclude/chr! Make sure the SNP ids are correct.")

  if(!export) {
    return(list(keep=keep, extract=extract,
                N=N, P=P, n=n, p=p, bfile=bfile,
                bimfile=bimfile, famfile=famfile,
                bim=NULL, fam=NULL))
  } else {
    return(list(keep=keep, extract=extract,
                N=N, P=P, n=n, p=p, bfile=bfile,
                bimfile=bimfile, famfile=famfile,
                bim=bim, fam=fam))
  }
  #' @return A list containing:
  #' \item{keep}{Either NULL or a logical vector of which individuals to keep}
  #' \item{extract}{Either NULL or a logical vector of which SNPs to extract}
  #' \item{N}{Number of rows in the PLINK bfile}
  #' \item{P}{Number of columns in the PLINK bfile}
  #' \item{n}{Number of rows in the PLINK bfile after keep}
  #' \item{p}{Number of columns in the PLINK bfile after extract}
  #' \item{bfile}{File path to bfile stub}
  #' \item{bimfile}{File path to bim file}
  #' \item{famfile}{File path to fam file}
  #' \item{bim}{Either \code{NULL} or a data frame of bim data}
  #' \item{fam}{Either \code{NULL} or a data frame of fam data}

}
