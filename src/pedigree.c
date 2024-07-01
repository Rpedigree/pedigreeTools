#include "pedigree.h"
#include <string.h>

#ifdef ENABLE_NLS		/** Allow for translation of error messages */
#include <libintl.h>
#define _(String) dgettext ("pedigreemm", String)
#else
#define _(String) (String)
#endif

/* When appropriate, alloca is cleaner than malloc/free.  The storage
 * is freed automatically on return from a function. When using gcc the
 * builtin version is much faster. */
#ifdef __GNUC__
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
/* this is necessary (and sufficient) for Solaris 10: */
# ifdef __sun
#  include <alloca.h>
# endif
#endif

/** alloca n elements of type t */
#define Alloca(n, t)   (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

/* Functions for the pedigree class */

/**
 * Create the left Cholesky factor of the numerator relationship
 * matrix from a pedigree.
 *
 * @param x a pedigree object
 * @param ans T stored as a non-unit, lower dtCMatrix
 *
 * @return ans with elements modified to incorporate D
 */
SEXP pedigree_chol(SEXP x, SEXP ans)
{
    SEXP Sire = GET_SLOT(x, install("sire"));
    int *ai = INTEGER(GET_SLOT(ans, install("i"))),
	*ap = INTEGER(GET_SLOT(ans, install("p"))),
	*dam = INTEGER(GET_SLOT(x, install("dam"))),
	*sire = INTEGER(Sire), j, n = LENGTH(Sire);
    double *ax = REAL(GET_SLOT(ans, install("x"))), *F, Di, tmp;

    SEXP tmp2 = PROTECT(allocVector(REALSXP, n));
    setAttrib(ans, install("F"), tmp2);

    F = REAL(getAttrib(ans, install("F")));
    for (int i = 0; i < n; i++) {
	int p = sire[i] - 1, q = dam[i] - 1;
	if (sire[i] == NA_INTEGER) {
	    F[i] = 0;
	    Di = (dam[i] == NA_INTEGER) ? 1 : sqrt(0.75 - 0.25 * F[q]);
	} else {
	    if (dam[i] == NA_INTEGER) { /* sire only */
		F[i] = 0;
		Di = sqrt(0.75 - 0.25 * F[p]);
	    } else {		/* both parents in pedigree */
		Di = sqrt(0.5 - 0.25 * (F[p] + F[q]));
		F[i] = NA_REAL;
		if ((ap[i + 1] - ap[i]) > 1) {	  /* skip if no progeny */
		    if (p > q) {j = p; p = q; q = j;} /* ensure p <= q */
		    F[i] = 0;
		    for (int j = 0; j <= p; j++) {
			tmp = 0;
			for (int k = ap[j]; k < ap[j + 1] && ai[k] <= q; k++) {
			    int ii = ai[k];
			    if (ii == p) tmp = ax[k];
			    if (ii == q) F[i] += tmp * ax[k]/2;
			}
		    }
		}
	    }
	}
	for (int j = ap[i]; j < ap[i + 1]; j++) ax[j] *= Di;
    }
    UNPROTECT(1);
    return ans;
}

/**
 * Create the inbreeding coefficients according to the algorithm given
 * in "Comparison of four direct algorithms for computing inbreeding
 * coefficients" by Mehdi Sargolzaei and Hiroaki Iwaisaki, Animal
 * Science Journal (2005) 76, 401--406.  This function is a modified
 * version of the code published in an appendix to that paper.
 *
 * @param x a pedigree object
 *
 * @return a list of the inbreeding coefficients
 */

SEXP pedigree_inbreeding(SEXP x)
{
    SEXP ans, sp = GET_SLOT(x, install("sire"));
    int i, j, t, n = LENGTH(sp), S, D;
    int *SI, *MI,		/* start and minor */
	*sire = INTEGER(sp),
	*dam = INTEGER(GET_SLOT(x, install("dam")));
    double *F = Calloc(n + 1, double), /* inbreeding coefficients */
      *L = Calloc(n + 1, double), *B = Calloc(n + 1, double);
    int *Anc = Calloc(n + 1, int),	/* ancestor */
	*LAP = Calloc(n + 1, int); 	/* longest ancestoral path */
    R_CheckStack();
    
    for (i = 0; i < n; i++) {     /* Replace NA's by zeros */
	if (sire[i] == NA_INTEGER) sire[i] = 0;
	if (dam[i] == NA_INTEGER) dam[i] = 0;
    }
    F[0] =-1; LAP[0] =-1; /* set F and lap for unknown parents */
    for(i = 1, t = -1; i <= n; i++) { 	/* evaluate LAP and its maximum */
	S = sire[i-1]; D = dam[i-1]; /* parents of animal i-1 */
	LAP[i] = ((LAP[S] < LAP[D]) ? LAP[D] : LAP[S]) + 1;
	if (LAP[i] > t) t = LAP[i];
    }
    SI = Calloc(t + 1, int);
    MI = Calloc(t + 1, int);
    for(i = 0; i <= t ; ++i) SI[i] = MI[i] = 0; /* initialize start and minor */
    for(i = 1; i <= n; i++) { 	/* evaluate F */
	S = sire[i-1]; D = dam[i-1]; /* parents of animal i */
	B[i] = 0.5 - 0.25 * (F[S] + F[D]); 
				/* adjust start and minor */
	for (j = 0; j < LAP[i]; j++) {++SI[j]; ++MI[j];} 
	if (S == 0 || D == 0) { /* both parents unknown */
	    F[i] = L[i] = 0; continue;
	}
	if(S == sire[i-2] && D == dam[i-2]) { /* full-sib with last animal */
	    F[i] = F[i-1]; L[i] = L[i-1]; continue;
	}
    
	F[i] = -1; L[i] = 1; 
	t = LAP[i]; /* largest lap group number in the animal's pedigree */
	Anc[MI[t]++] = i; /* initialize Anc and increment MI[t] */
	while(t > -1) { /* from the largest lap group to zero */
	    j = Anc[--MI[t]]; /* next ancestor */
	    S = sire[j-1]; D = dam[j-1]; /* parents of the ancestor */
	    if (S) {
		if (!L[S]) Anc[MI[LAP[S]]++] = S; 
				/* add sire in its lap group in Anc
				 * array if it is not added yet and
				 * increment the minor index for the group */ 
		L[S] += 0.5 * L[j]; /* contribution to sire */
	    }
	    if (D) {
		if (!L[D]) Anc[MI[LAP[D]]++] = D;
		L[D] += 0.5 * L[j]; /* contribution to dam */
	    }
	    F[i] += L[j] * L[j] * B[j];
	    L[j] = 0; /*clear L[j] for the evaluation of the next animal */
	    if (MI[t] == SI[t]) --t; /* move to the next lap group when
				      * all ancestors in group t have been
				      * evaluated */
	} 
    }
    ans = PROTECT(allocVector(REALSXP, n));
    Memcpy(REAL(ans), F + 1, n);

    for (i = 0; i < n; i++) {     /* Restore the NA's */
	if (!sire[i]) sire[i] = NA_INTEGER;
	if (!dam[i]) dam[i] = NA_INTEGER;
    }
    Free(F); Free(L); Free(B); Free(Anc); Free(LAP); Free(SI); Free(MI);
    UNPROTECT(1);
    return ans;
}

/**
 * @title Expand Pedigree for Selfing Generations
 * 
 * @description This function expands a pedigree to account for selfing generations.
 * It creates new entries for each selfing generation of an individual.
 *
 * @param labels SEXP (STRSXP) Vector of individual labels
 * @param sires SEXP (STRSXP) Vector of sire labels
 * @param dams SEXP (STRSXP) Vector of dam labels
 * @param selfing_generations SEXP (INTSXP) Vector of selfing generation numbers
 * @param sep_char SEXP (STRSXP) Separator character for expanded pedigree IDs
 * @param verbose SEXP (STRSXP) print progress
 *
 * @return SEXP (VECSXP) A list containing:
 *         - labels: Expanded vector of labels
 *         - sires: Expanded vector of sire labels
 *         - dams: Expanded vector of dam labels
 *         - generations: Vector of generation numbers (all NA)
 *         - selfing_generations: Expanded vector of selfing generation numbers
 *         - expanded: Logical vector indicating if each entry is an expansion
 *
 */
SEXP expand_pedigree_selfing(SEXP labels, SEXP sires, SEXP dams, SEXP selfing_generations, SEXP sep_char, SEXP verbose) {
    int n = LENGTH(labels);
    int total_rows = 0;
    int *sg = INTEGER(selfing_generations);
    int show_progress = asLogical(verbose);
    
    // Calculate total number of rows in expanded pedigree
    for (int i = 0; i < n; i++) {
        total_rows += (sg[i] == 0) ? 1 : (sg[i] + 1);
    }
    
    // Allocate memory for result
    SEXP result = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(result, 0, allocVector(STRSXP, total_rows));  // label
    SET_VECTOR_ELT(result, 1, allocVector(STRSXP, total_rows));  // sire
    SET_VECTOR_ELT(result, 2, allocVector(STRSXP, total_rows));  // dam
    SET_VECTOR_ELT(result, 3, allocVector(INTSXP, total_rows));  // generation
    SET_VECTOR_ELT(result, 4, allocVector(INTSXP, total_rows));  // selfing_generation
    SET_VECTOR_ELT(result, 5, allocVector(LGLSXP, total_rows));  // expanded
    
    const char* sep = CHAR(STRING_ELT(sep_char, 0));
    int row = 0;
    char buffer[256];  // Adjust buffer size as needed
    
    if (show_progress) {
        R_FlushConsole();
    }
    
    for (int i = 0; i < n; i++) {
        // Update progress
        if (show_progress) print_progress_bar(i + 1, n, 50);
        const char* id = CHAR(STRING_ELT(labels, i));
        int cycles = sg[i];
        
        if (cycles == 0) {
            SET_STRING_ELT(VECTOR_ELT(result, 0), row, STRING_ELT(labels, i));
            SET_STRING_ELT(VECTOR_ELT(result, 1), row, STRING_ELT(sires, i));
            SET_STRING_ELT(VECTOR_ELT(result, 2), row, STRING_ELT(dams, i));
            INTEGER(VECTOR_ELT(result, 3))[row] = NA_INTEGER;
            INTEGER(VECTOR_ELT(result, 4))[row] = 0;
            LOGICAL(VECTOR_ELT(result, 5))[row] = FALSE;
            row++;
        } else {
            for (int j = 0; j <= cycles; j++) {
                if (j == cycles) {
                    SET_STRING_ELT(VECTOR_ELT(result, 0), row, STRING_ELT(labels, i));
                    LOGICAL(VECTOR_ELT(result, 5))[row] = FALSE;
                } else {
                    snprintf(buffer, sizeof(buffer), "%s%s%d", id, sep, j);
                    SET_STRING_ELT(VECTOR_ELT(result, 0), row, mkChar(buffer));
                    LOGICAL(VECTOR_ELT(result, 5))[row] = TRUE;
                }
                
                if (j == 0) {
                    SET_STRING_ELT(VECTOR_ELT(result, 1), row, STRING_ELT(sires, i));
                    SET_STRING_ELT(VECTOR_ELT(result, 2), row, STRING_ELT(dams, i));
                } else {
                    snprintf(buffer, sizeof(buffer), "%s%s%d", id, sep, j-1);
                    SET_STRING_ELT(VECTOR_ELT(result, 1), row, mkChar(buffer));
                    SET_STRING_ELT(VECTOR_ELT(result, 2), row, mkChar(buffer));
                }
                
                INTEGER(VECTOR_ELT(result, 3))[row] = NA_INTEGER;
                INTEGER(VECTOR_ELT(result, 4))[row] = j;
                row++;
            }
        }
        
    }
    
    if (show_progress) {
		Rprintf("\n");
        R_FlushConsole();
    }
    
    UNPROTECT(1);
    return result;
}

// for printing progress in c functions
void print_progress_bar(int current, int total, int bar_width) {
    float progress = (float)current / total;
    int filled_width = (int)(bar_width * progress);

    Rprintf("\r[");
    for (int i = 0; i < bar_width; ++i) {
        if (i < filled_width) {
            Rprintf("=");
        } else {
            Rprintf(" ");
        }
    }
    Rprintf("] %.1f%%", progress * 100);
    R_FlushConsole();
}
