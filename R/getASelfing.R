#' Extends the pedigree according to number of selfing cycles
#' and also optionally computes the Additive Relationship Matrix for that pedigree.
#' @param ID is a vector of individual IDs
#' @param Par1 vector of IDs of one of the parents
#' @param Par2 vector of IDs of the other parent
#' @param nCycles vector that indicates number of selfing cycles for each individual.
#' @param sepChar character, used for expanded pedigree IDs
#' @param verbose logical, print progress
#' @param fileNewPed Output csv file (comma separated value) with columns 'label', 'sire', 'dam', with the full pull pedigree expanded taking into account the selfing cycles
#' @param computeA Indicates if the A matrix is to be computed
#' @return Returns A matrix computed for the extended pedigree if computeA=TRUE
#' @export
getASelfing = function(ID, Par1, Par2, nCycles, verbose = FALSE, sepChar = "-F", fileNewPed = NULL, computeA = TRUE, ...) {
                       
    ## Create initial pedigree object
    ped <- pedigree(sire = Par1, dam = Par2, label = ID, selfing_generation = nCycles)
    
    ## Expand pedigree to account for selfing generations
    expandedPed <- expandPedigreeSelfing(ped, sepChar = sepChar, verbose = verbose)
    
    ## Write expanded pedigree to file if requested
    if (!is.null(fileNewPed)) {
        write.csv(ped2DF(expandedPed), file = fileNewPed, row.names = FALSE)
    }
    
    if (computeA) {
        ## Compute A matrix
        A = getA(expandedPed, labs = expandedPed@label[expandedPed@expanded == FALSE])
        return(A)
    }
    
    ## If computeA is FALSE, return the expanded pedigree
    return(expandedPed)
}
