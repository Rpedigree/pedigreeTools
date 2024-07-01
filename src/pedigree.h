#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>

SEXP pedigree_chol(SEXP x, SEXP ans);
SEXP pedigree_inbreeding(SEXP x);
SEXP expand_pedigree_selfing(SEXP labels, SEXP sires, SEXP dams, SEXP selfing_generations, SEXP sep_char);


#endif /* PEDIGREE_H */
