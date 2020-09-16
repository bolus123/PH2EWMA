#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _PH2EWMA_integrandSteady(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PH2EWMA_integrandZero(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PH2EWMA_qintegrandSteady(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PH2EWMA_qintegrandZero(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_PH2EWMA_integrandSteady",  (DL_FUNC) &_PH2EWMA_integrandSteady,   9},
    {"_PH2EWMA_integrandZero",    (DL_FUNC) &_PH2EWMA_integrandZero,    10},
    {"_PH2EWMA_qintegrandSteady", (DL_FUNC) &_PH2EWMA_qintegrandSteady, 10},
    {"_PH2EWMA_qintegrandZero",   (DL_FUNC) &_PH2EWMA_qintegrandZero,   11},
    {NULL, NULL, 0}
};

void R_init_PH2EWMA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}