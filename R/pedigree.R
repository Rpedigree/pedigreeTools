#### "pedigree" class methods

#' @title Constructor for pedigree objects
#'
#' @description A simple constructor for a pedigree object. The main point for
#'   the constructor is to use coercions to make the calls easier.
#'
#' @param sire integer vector or factor representation of the sires
#' @param dam integer vector or factor representation of the dams
#' @param label character vector of individual labels
#' @return an pedigree object of class \linkS4class{pedigree}
#' @note \code{sire}, \code{dam} and \code{label} must all have the
#'   same length and all labels in \code{sire} and \code{dam} must occur
#'   in \code{label}
#' @export
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
#' ped
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

#' @title Convert a pedigree to a data frame
#'
#' @description Express a pedigree as a data frame with \code{sire} and
#'   \code{dam} stored as factors. If the pedigree is an object of
#'   class \code{\link{pedinbred}} then the inbreeding coefficients are
#'   appended as the variable \code{F}
#'
#' @param x \code{\link{pedigree}}
#' @return a data frame
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
#' ped2DF(ped)
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

#' @title Inbreeding coefficients from a pedigree
#'
#' @description Create the inbreeding coefficients according to the algorithm
#'   given in "Comparison of four direct algorithms for computing inbreeding
#'   coefficients" by Mehdi Sargolzaei and Hiroaki Iwaisaki, Animal Science
#'   Journal (2005) 76, 401--406.
#'
#' @param ped \code{\link{pedigree}}
#' @return the inbreeding coefficients as a numeric vector
#' @export
#' @useDynLib pedigreeTools pedigree_inbreeding
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
#' inbreeding(ped)
inbreeding <- function(ped) {
    stopifnot(is(ped, "pedigree"))
    .Call(pedigree_inbreeding, ped)
}

#' @title Diagonal of D in the A = TDT' factorization.
#'
#' @description Determine the diagonal factor in the decomposition of the
#'   relationship matrix A as TDT' where T is unit lower triangular.
#'
#' @param ped \code{\link{pedigree}}
#' @return a numeric vector
#' @export
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
#' (D <- Dmat(ped))
#'
#' # Test for correctness
#' DExp <- c(1.00, 1.00, 0.50, 0.75, 0.50, 0.46875)
#' stopifnot(!any(abs(D - DExp) > .Machine$double.eps))
Dmat <- function(ped) {
    F <- inbreeding(ped)
    sire <- ped@sire
    dam <- ped@dam
    Fsire <- ifelse(is.na(sire), -1, F[sire])
    Fdam <- ifelse(is.na(dam), -1, F[dam])
    ans <- 1 - 0.25 * (2 + Fsire + Fdam)
    names(ans) <- ped@label
    ans
}

#' @title Relationship factor from a pedigree
#'
#' @description Determine the right Cholesky factor of the relationship matrix
#'   for the pedigree \code{ped}, possibly restricted to the specific labels
#'   that occur in \code{labs}.
#'
#' @param ped \code{\link{pedigree}}
#' @param labs a character vector or a factor giving individual labels to
#'   which to restrict the relationship matrix. If \code{labs} is a
#'   factor then the levels of the factor are used as the labels.
#'   Default is the complete set of individuals in the pedigree.
#' @return matrix (\linkS4class{dtCMatrix} - upper triangular sparse)
#' @export
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
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
              as(factor(ped@label, levels = ped@label), "sparseMatrix"))
    tmpA <- Matrix::crossprod(rect)
    tmp <- ped@label %in% labs
    tmpA <- tmpA[tmp, tmp]

    labped <- ped@label[tmp]
    orlab <- order(as.numeric(factor(labped, levels=labs, ordered = TRUE)))
    labped <- as.character(labped[orlab])
    tmpA <- tmpA[orlab, orlab]
    stopifnot(all.equal(as.character(labped), as.character(labs)))
    relf <- Matrix::chol(tmpA)
    dimnames(relf) <- list(labs, labs)
    relf
}

# TODO: Add a test for L matrix

#' @title Inverse relationship factor from a pedigree
#'
#' @description Get inverse of the right Cholesky factor of the relationship
#'   matrix for the pedigree \code{ped}.
#'
#' @param ped \code{\link{pedigree}}
#' @return matrix (\linkS4class{dtCMatrix} - triangular sparse)
#' @export
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
#' getRelFactorInv(ped)
getRelFactorInv <- function(ped) {
    stopifnot(is(ped, "pedigree"))
    T_Inv <- as(ped, "sparseMatrix") # dtCMatrix (lower triangular sparse)
    DSq_Inv <- Matrix::Diagonal(x = 1 / sqrt(Dmat(ped))) # ddiMatrix (diagonal sparse)
    L_Inv <- T_Inv %*% DSq_Inv # dtCMatrix (triangular sparse)
    dimnames(L_Inv) <- list(ped@label, ped@label)
    L_Inv
}
# return(Matrix::Diagonal(x = sqrt(Dmat(ped))) %*%
# Matrix::solve(Matrix::t(as(ped, "sparseMatrix"))))

# TODO: Add a test for L_Inv matrix

#' @title Inverse of the Additive Relationship Matrix
#'
#' @description Returns the inverse of additive relationship matrix for the
#'   pedigree.
#'
#' @param ped \code{\link{pedigree}}
#' @return matrix (\linkS4class{dsCMatrix} - symmetric sparse)
#' @export
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
#' getAInv(ped)
#'
#' # Test for correctness
#' AInvExp <- matrix(data = c(1.0000, 0.0000, 0.5000, 0.5000, 0.5000, 0.2500,
#'                            0.0000, 1.0000, 0.5000, 0.0000, 0.2500, 0.6250,
#'                            0.5000, 0.5000, 1.0000, 0.2500, 0.6250, 0.5625,
#'                            0.5000, 0.0000, 0.2500, 1.0000, 0.6250, 0.3125,
#'                            0.5000, 0.2500, 0.6250, 0.6250, 1.1250, 0.6875,
#'                            0.2500, 0.6250, 0.5625, 0.3125, 0.6875, 1.1250),
#'                   byrow = TRUE, nrow = 6)
#' stopifnot(!any(abs(AInv - AInvExp) > .Machine$double.eps))
getAInv <- function(ped)
{
    stopifnot(is(ped, "pedigree"))
    A_Inv <- Matrix::crossprod(getRelFactorInv(ped)) # dsCMatrix (symmetric sparse)
    dimnames(A_Inv) <- list(ped@label, ped@label)
    A_Inv
}

# TODO: Add a test for AInv matrix
# https://github.com/Rpedigree/pedigreeTools/issues/3

#' @title Additive Relationship Matrix
#'
#' @description Returns the additive relationship matrix for the pedigree.
#'
#' @param ped \code{\link{pedigree}}
#' @return matrix (\linkS4class{dsCMatrix} - symmetric sparse)
#' @export
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
#' (A <- getA(ped))
#'
#' # Test for correctness
#' AExp <- matrix(data = c(1.0000, 0.0000, 0.5000, 0.5000, 0.5000, 0.2500,
#'                         0.0000, 1.0000, 0.5000, 0.0000, 0.2500, 0.6250,
#'                         0.5000, 0.5000, 1.0000, 0.2500, 0.6250, 0.5625,
#'                         0.5000, 0.0000, 0.2500, 1.0000, 0.6250, 0.3125,
#'                         0.5000, 0.2500, 0.6250, 0.6250, 1.1250, 0.6875,
#'                         0.2500, 0.6250, 0.5625, 0.3125, 0.6875, 1.1250),
#'                byrow = TRUE, nrow = 6)
#' stopifnot(!any(abs(A - AExp) > .Machine$double.eps))
getA <- function(ped) {
    stopifnot(is(ped, "pedigree"))
    aMx <- Matrix::crossprod(relfactor(ped))
    dimnames(aMx) <- list(ped@label, ped@label)
    aMx
}

# TODO: Add a test for A matrix

#' @title Counts number of generations of ancestors for one subject. Use recursion.
#'
#' @param pede data frame with a pedigree and a column for the number of
#'   generations of each subject.
#' @param id subject for which we want the number of generations.
#' @param ngen number of generation
#' @return a data frame object with the pedigree and generation of
#'   ancestors for subject id.
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
#' getGenAncestors(ped)
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

#' @title Edits a disordered or incomplete pedigree
#'
#' 1_ add labels for the sires and dams not listed as labels before.
#' 2_ order pedigree based on recursive calls to \code{\link{getGenAncestors}}.
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

#' @title Subsets a pedigree for a specified vector of individuals up to a
#' specified number of previous generations using recursion.

#' @param ped Data Frame pedigree to be subset
#' @param selectVector Vector of individuals to select from pedigree
#' @param ngen Number of previous generations of parents to select starting from selectVector.

#' @return Returns Subsetted pedigree as a DataFrame.
#' @export
#' @examples
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
#'                 dam =  c(NA, NA, 2, NA, 3, 2),
#'                 label = 1:6)
prunePed <- function(ped,selectVector,ngen=2){

  ped <- as.matrix(ped)

  returnPed <- matrix(c(NA,NA,NA),nrow=1,ncol=3)

  findBase <- ped[,"label"] %in% selectVector
  basePed <- ped[findBase,]
  findSire <- ped[,"label"] %in% basePed[,"sire"]
  findDam <- ped[,"label"] %in% basePed[,"dam"]

  newSelVec <- ped[findSire|findDam,"label"]
  newSelVec <- newSelVec[!(newSelVec %in% selectVector)]

  if (ngen != -1) {
    returnPed <- basePed
    returnPed <- unique(rbind(returnPed,prunePed(ped,newSelVec,ngen-1)))
    returnPed <- returnPed[rowSums(is.na(returnPed))!=3,]
  } else {
    return(returnPed)
  }

  return(as.data.frame(returnPed))
}
