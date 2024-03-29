% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pedigree.R
\name{prunePed}
\alias{prunePed}
\title{Subsets a pedigree for a specified vector of individuals up to a
specified number of previous generations using recursion.}
\usage{
prunePed(ped, selectVector, ngen = 2)
}
\arguments{
\item{ped}{Data Frame pedigree to be subset}

\item{selectVector}{Vector of individuals to select from pedigree}

\item{ngen}{Number of previous generations of parents to select starting from selectVector.}
}
\value{
Returns Subsetted pedigree as a DataFrame.
}
\description{
Subsets a pedigree for a specified vector of individuals up to a
specified number of previous generations using recursion.
}
\examples{
ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
                dam =  c(NA, NA, 2, NA, 3, 2),
                label = 1:6)
}
