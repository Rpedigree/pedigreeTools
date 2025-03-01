#' Pedigree class
#'
#' @export
setClass("pedigree", 
    representation = list(
        sire = "integer", 
        dam = "integer", 
        label = "character", 
        generation = "integer",
        selfing_generation = "integer",
        expanded = "logical"  # New slot
    ),
    prototype = list(
        generation = integer(0),
        selfing_generation = integer(0),
        expanded = logical(0)  # Default empty logical vector
    ),
    validity = function(object) {
        n <- length(sire <- object@sire)
        if (length(dam <- object@dam) != n)
            return("sire and dam slots must be the same length")
        if (length(object@label) != n)
            return("'label' slot must have the same length as 'sire' and 'dam'")
        if (length(object@generation) > 0 && length(object@generation) != n)
            return("If provided, 'generation' slot must have the same length as 'sire' and 'dam'")
        if (length(object@selfing_generation) > 0 && length(object@selfing_generation) != n)
            return("If provided, 'selfing_generation' slot must have the same length as 'sire' and 'dam'")
        if (length(object@expanded) > 0 && length(object@expanded) != n)
            return("If provided, 'expanded' slot must have the same length as 'sire' and 'dam'")
        if(n == 0) return(TRUE)
        animal <- 1:n
        snmiss <- !is.na(sire)
        dnmiss <- !is.na(dam)
        if (any(sire[snmiss] >= animal[snmiss]) ||
            any(dam[dnmiss]  >= animal[dnmiss]))
            return("the sire and dam must precede the offspring")
        if (any(sire[snmiss] < 1 | sire[snmiss] > n) |
            any(dam[dnmiss] < 1 | dam[dnmiss] > n))
            return(paste("Non-missing sire or dam must be in [1,",
                         n, "]", sep = ''))
        TRUE
    })
