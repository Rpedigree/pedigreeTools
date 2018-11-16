#### "pedigree" class methods

#' Constructor for pedigree objects
#'
#' A simple constructor for a pedigree object.  The main point for the
#' constructor is to use coercions to make the calls easier.
#'
#' @param sire integer vector or factor representation of the sires
#' @param dam integer vector or factor representation of the dams
#' @param label character vector of labels
#' @return an pedigree object of class \linkS4class{pedigree}
#' @note \code{sire}, \code{dam} and \code{label} must all have the
#'   same length and all labels in \code{sire} and \code{dam} must occur
#'   in \code{label}
#' @export
#' @examples
#' ped <- pedigree(sire=c(NA,NA,1,1,4,5), dam=c(NA,NA,2,NA,3,2), label=1:6)
pedigree <- function(sire, dam, label) {
    n <- length(sire)
    labelex <- c(label, NA, 0)
    stopifnot(n == length(dam),
              n == length(label),
              all(sire %in% labelex),
              all(dam %in% labelex))
    sire <- as.integer(factor(sire, levels = label))
    dam <- as.integer(factor(dam, levels = label))
    sire[sire < 1 | sire > n] <- NA
    dam[dam < 1 | dam > n] <- NA
    new("pedigree", sire = sire, dam = dam,
        label = as.character(label))
}

setAs("pedigree", "sparseMatrix", # representation as T^{-1}
      function(from) {
	  sire <- from@sire
	  n <- length(sire)
	  animal <- seq_along(sire)
	  j <- c(sire, from@dam)
	  ind <- !is.na(j)
	  as(new("dtTMatrix", i = rep.int(animal, 2)[ind] - 1L,
		 j = j[ind] - 1L, x = rep.int(-0.5, sum(ind)),
		 Dim = c(n,n), Dimnames = list(from@label, NULL),
		 uplo = "L", diag = "U"), "dtCMatrix")
      })

## these data frames are now storage efficient but print less nicely
setAs("pedigree", "data.frame",
      function(from)
      data.frame(sire = from@sire, dam = from@dam,
		 row.names = from@label))

#' Convert a pedigree to a data frame
#'
#' Express a pedigree as a data frame with \code{sire} and
#' \code{dam} stored as factors.  If the pedigree is an object of
#' class pedinbred then the inbreeding coefficients are
#' appended as the variable \code{F}
#'
#' @param x a pedigree object of class \linkS4class{pedigree}
#' @return a data frame
ped2DF <- function(x) {
    stopifnot(is(x, "pedigree"))
    lab <- x@label
    lev <- seq_along(lab)
    ans <- data.frame(sire = factor(x@sire, levels = lev, labels = lab),
                      dam  = factor(x@dam,  levels = lev, labels = lab),
                      row.names = lab)
    if (is(x, "pedinbred")) ans <- cbind(ans, F = x@F)
    ans
}

setMethod("show", signature(object = "pedigree"),
	  function(object) print(ped2DF(object)))

setMethod("head", "pedigree", function(x, ...)
	  do.call("head", list(x = ped2DF(x), ...)))

setMethod("tail", "pedigree", function(x, ...)
	  do.call("tail", list(x = ped2DF(x), ...)))

#' @useDynLib pedigreeTools pedigree_chol
setMethod("chol", "pedigree",
          function(x, pivot, LINPACK) {
              ttrans <- Matrix::solve(Matrix::t(as(x, "dtCMatrix")))
              .Call(pedigree_chol, x,
                    as(.Call("Csparse_diagU2N", Matrix::t(ttrans), PACKAGE = "Matrix"),
                       "dtCMatrix"))
          })

#' Inbreeding coefficients from a pedigree
#'
#' Create the inbreeding coefficients according to the algorithm given
#' in "Comparison of four direct algorithms for computing inbreeding
#' coefficients" by Mehdi Sargolzaei and Hiroaki Iwaisaki, Animal
#' Science Journal (2005) 76, 401--406.
#'
#' @param ped an object that inherits from class \linkS4class{pedigree}
#' @return the inbreeding coefficients as a numeric vector
#' @export
#' @useDynLib pedigreeTools pedigree_inbreeding
#' @examples
#' ped <- pedigree(sire=c(NA,NA,1,1,4,5), dam=c(NA,NA,2,NA,3,2), label=1:6)
#' inbreeding(ped)
inbreeding <- function(ped) {
    stopifnot(is(ped, "pedigree"))
    .Call(pedigree_inbreeding, ped)
}

#' Diagonal of D in the A = TDT' factorization.
#'
#' Determine the diagonal factor in the decomposition of the
#' relationship matrix A as TDT' where T is unit lower triangular.
#'
#' @param ped an object that inherits from class \linkS4class{pedigree}
#' @return a numeric vector
#' @export
#' @examples
#' ped <- pedigree(sire=c(NA,NA,1,1,4,5), dam=c(NA,NA,2,NA,3,2), label=1:6)
#' Dmat(ped)
Dmat <- function(ped)
{
    F <- inbreeding(ped)
    sire <- ped@sire
    dam <- ped@dam
    Fsire <- ifelse(is.na(sire), -1, F[sire])
    Fdam <-  ifelse(is.na(dam), -1, F[dam])
    ans <- 1 - 0.25 * (2 + Fsire + Fdam)
    names(ans) <- ped@label
    ans
}

#' Relationship factor from a pedigree
#'
#' Determine the right Cholesky factor of the relationship matrix for
#' the pedigree \code{ped}, possibly restricted to the specific labels
#' that occur in \code{labs}.
#'
#' @param ped a pedigree that includes the individuals who occur in svec
#' @param labs a character vector or a factor giving the labels to
#'   which to restrict the relationship matrix. If \code{labs} is a
#'   factor then the levels of the factor are used as the labels.
#'   Default is the complete set of labels in the pedigree.
#' @return an object that inherits from \linkS4class{CHMfactor}
#' @export
#' @examples
#' ped <- pedigree(sire=c(NA,NA,1,1,4,5), dam=c(NA,NA,2,NA,3,2), label=1:6)
#' relfactor(ped)
relfactor <- function(ped, labs)
{
    stopifnot(is(ped, "pedigree"))
    if (missing(labs))                  # square case
        return(Matrix::Diagonal(x = sqrt(Dmat(ped))) %*%
               Matrix::solve(Matrix::t(as(ped, "sparseMatrix"))))
    labs <- factor(labs) # drop unused levels from a factor
    stopifnot(all(labs %in% ped@label))
    rect <- Matrix::Diagonal(x = sqrt(Dmat(ped))) %*%
        Matrix::solve(Matrix::t(as(ped, "sparseMatrix")), # rectangular factor
              as(factor(ped@label, levels = ped@label),"sparseMatrix"))
    tmpA<-Matrix::crossprod(rect)
    tmp<- ped@label %in% labs
    tmpA<-tmpA[tmp,tmp]

    orlab <- order(as.numeric(factor(labped<-ped@label[tmp], levels=labs, ordered=T)))
    labped<- as.character(labped[orlab])
    tmpA  <- tmpA[orlab, orlab]
    stopifnot(all.equal(as.character(labped), as.character(labs)))
    relf<-Matrix::chol(tmpA)
    dimnames(relf)[[1]]<- dimnames(relf)[[2]]<-labs
    relf
}

#' Inverse of the Relationship Matrix
#'
#' @param ped a pedigree that includes the individuals who occur in svec
#'   which to restrict the relationship matrix. If \code{labs} is a
#'   factor then the levels of the factor are used as the labels.
#'   Default is the complete set of labels in the pedigree.
#' @return an object that inherits from \linkS4class{CHMfactor}
#' @export
#' @examples
#' ped <- pedigree(sire=c(NA,NA,1,1,4,5), dam=c(NA,NA,2,NA,3,2), label=1:6)
#' getAInv(ped)
getAInv <- function(ped)
{
    stopifnot(is(ped, "pedigree"))
    T_Inv <- as(ped, "sparseMatrix")
    D_Inv <- Matrix::diag(1/Dmat(ped))
    aiMx<-Matrix::t(T_Inv) %*% D_Inv %*% T_Inv
    dimnames(aiMx)[[1]]<-dimnames(aiMx)[[2]] <-ped@label
    aiMx
}

#' Additive Relationship Matrix
#'
#' Returns the additive relationship matrix for the pedigree \code{ped}.
#'
#' @param ped a pedigree that includes the individuals who occur in svec
#'   which to restrict the relationship matrix. If \code{labs} is a
#'   factor then the levels of the factor are used as the labels.
#'   Default is the complete set of labels in the pedigree.
#' @return an object that inherits from \linkS4class{CHMfactor}
#' @keywords array
#' @export
#' @examples
#' ped <- pedigree(sire=c(NA,NA,1,1,4,5), dam=c(NA,NA,2,NA,3,2), label=1:6)
#' getA(ped)
getA <- function(ped)
{
    stopifnot(is(ped, "pedigree"))
    aMx<-Matrix::crossprod(relfactor(ped))
    dimnames(aMx)[[1]]<-dimnames(aMx)[[2]] <-ped@label
    aMx
}

#' Counts number of generations of ancestors for one subject. Use recursion.
#'
#' @param pede data frame with a pedigree and a column for the number of
#'   generations of each subject.
#' @param id subject for which we want the number of generations.
#' @param ngen number of generation
#' @return a data frame object with the pedigree and generation of
#'   ancestors for subject id.
getGenAncestors <- function(pede, id, ngen=NULL){
    j <- which(pede$id==id)
    parents <- c(pede$sire[j], pede$dam[j])
    parents <- parents[!is.na(parents)]
    np <- length(parents)
    if(np==0)
    {
        pede$gene[j] <-0
        return(pede)
    }
    ## get the number of generations in parent number one
    tmpgenP1 <- pede$gene[pede$id==parents[1]]
    if( is.na(tmpgenP1))
    {
        #if ngen is not null, and not cero, ngen<- ngen-1
        #if ngen is cero, do not call recurrsively anymore
        pede <- getGenAncestors(pede, parents[1])
        genP1  <- 1 + pede$gene[pede$id==parents[1]]
    }  else {
        genP1 <- 1 + tmpgenP1
    }
    ## find out if there is a parent number two
    if (np==2){
        tmpgenP2 <- pede$gene[pede$id==parents[2]]
        if( is.na(tmpgenP2))
        {
            pede <- getGenAncestors(pede, parents[2])
            genP2  <- 1 + pede$gene[pede$id==parents[2]]
        }  else {
            genP2 <- 1 + tmpgenP2
        }
        genP1 <- max(genP1, genP2)
    }
    pede$gene[j] <- genP1
    ## print(paste('id:', id, ', gen:', genP1, ', row:', j))
    pede
}

#' Edits a disordered or incomplete pedigree.
#'
#' 1_ add labels for the sires and dams not listed as labels before.
#' 2_ order pedigree based on recursive calls to getGenAncestors.
#'
#' @param sire integer vector or factor representation of the sires
#' @param dam integer vector or factor representation of the dams
#' @param label character vector of labels
#' @param verbose logical to print the row of the pedigree that the
#'   function is ordering. Default is FALSE.
#' @return a data frame with the pedigree ordered.
#' @export
#' @examples
#' pede <- data.frame(sire=as.character(c(NA,NA,NA,NA,NA,1,3,5,6,4,8,1,10,8)),
#'                    dam=as.character(c(NA,NA,NA,NA,NA,2,2,NA,7,7,NA,9,9,13)),
#'                    label=as.character(1:14))
#' pede <- pede[sample(replace=FALSE, 1:14),]
#' pede <- editPed(sire=pede$sire, dam=pede$dam, label=pede$label)
#' ped <- with(pede, pedigree(label=label, sire=sire, dam=dam))
editPed <- function(sire, dam, label, verbose = FALSE)
{
    nped <- length(sire)
    if (nped != length(dam))  stop("sire and dam have to be of the same length")
    if (nped != length(label)) stop("label has to be of the same length than sire and dam")
    tmp <- unique(sort(c(as.character(sire), as.character(dam))))

    missingP <-NULL
    if(any(completeId <- ! tmp %in% as.character(label))) missingP <- tmp[completeId]
    labelOl <- c(as.character(missingP),as.character(label))
    sireOl <- c(rep(NA, times=length(missingP)),as.character(sire))
    damOl  <- c(rep(NA, times=length(missingP)),as.character(dam))
    sire <- as.integer(factor(sireOl, levels = labelOl))
    dam  <- as.integer(factor(damOl, levels = labelOl))
    nped <-length(labelOl)
    label <-1:nped
    sire[!is.na(sire) & (sire<1 | sire>nped)] <- NA
    dam[!is.na(dam) & (dam < 1 | dam > nped)] <- NA
    pede <- data.frame(id= label, sire= sire, dam= dam, gene=rep(NA, times=nped))
    noParents <- (is.na(pede$sire) & is.na(pede$dam))
    pede$gene[noParents] <- 0
    for(i in 1:nped){
        if(verbose) print(i)
        if(is.na(pede$gene[i])){
            id <-pede$id[i]
            pede <-getGenAncestors(pede, id)
        }}
    ord<- order(pede$gene)
    ans<-data.frame(label=labelOl, sire=sireOl, dam=damOl, gene=pede$gene,
                    stringsAsFactors =F)
    ans[ord,]
}

#' Subsets a pedigree for a specified vector of individuals upto a 
#' specified number of previous generations using Recursion.
  
#' @param ped Data Frame pedigree to be subset
#' @param selectVector Vector of individuals to select from pedigree
#' @param ngen Number of previous generations of parents to select starting from selectVector. 
  
#' @return Returns Subsetted pedigree as a DataFrame. 
#' @export


prunePed <- function(ped,selectVector,ngen=2){
  
  ped <- as.matrix(ped)

  returnPed <- matrix(c(NA,NA,NA),nrow=1,ncol=3)
                          
  findBase <- ped[,"label"] %in% selectVector 
  basePed <- ped[findBase,]
  findSire <- ped[,"label"] %in% basePed[,"sire"]
  findDam <- ped[,"label"] %in% basePed[,"dam"]
  
  newSelVec <- ped[findSire|findDam,"label"]
  newSelVec <- newSelVec[!(newSelVec %in% selectVector)]
  
  if(ngen!=-1){
  returnPed <- basePed
  returnPed <- unique(rbind(returnPed,prunePed(ped,newSelVec,ngen-1)))
  returnPed <- returnPed[rowSums(is.na(returnPed))!=3,]
  }else{return(returnPed)}

return(as.data.frame(returnPed))

}


