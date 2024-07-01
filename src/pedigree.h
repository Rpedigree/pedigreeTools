#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <R.h>
/* Rdefines.h includes Rinternals.h (for SEXP, REAL, etc.) and defines
 * GET_SLOT, MAKE_CLASS, NEW_OBJECT, SET_SLOT, etc. */
#include <Rdefines.h>
#include <R_ext/Print.h>

SEXP pedigree_chol(SEXP x, SEXP ans);
SEXP pedigree_inbreeding(SEXP x);
SEXP get_generation(SEXP sire, SEXP dam, SEXP label);
SEXP expand_pedigree_selfing(SEXP labels, SEXP sires, SEXP dams, SEXP selfing_generations, SEXP sep_char, SEXP verbose);
void print_progress_bar(int current, int total, int bar_width);

#endif /* PEDIGREE_H */
