% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pedigree.R
\name{ped2DF}
\alias{ped2DF}
\title{Convert a pedigree to a data frame}
\usage{
ped2DF(x)
}
\arguments{
\item{x}{\code{\link{pedigree}}}
}
\value{
a data frame
}
\description{
Express a pedigree as a data frame with \code{sire} and
  \code{dam} stored as factors. If the pedigree is an object of
  class \linkS4class{pedinbred} then the inbreeding coefficients are
  appended as the variable \code{F}
}
\examples{
ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5),
                dam =  c(NA, NA, 2, NA, 3, 2),
                label = 1:6)
ped2DF(ped)
}
