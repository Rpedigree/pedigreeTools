#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>
#include "pedigree.h"

static const R_CallMethodDef CallEntries[] = {
    {"pedigree_chol", (DL_FUNC) &pedigree_chol, 2},
    {"pedigree_inbreeding", (DL_FUNC) &pedigree_inbreeding, 1},
    {"get_generation", (DL_FUNC) &get_generation, 3},
    {"expand_pedigree_selfing", (DL_FUNC) &expand_pedigree_selfing, 6},
    {NULL, NULL, 0}
};

void R_init_pedigreeTools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
