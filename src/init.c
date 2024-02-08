#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */

extern void dec2binC(void *, void *, void *);
extern void bin2decC(void *, void *, void *);


/* .Call calls */

extern SEXP simulate_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP simulate_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP simulate_async_return_states_R(SEXP);

static const R_CMethodDef cMethods[] = {
  {"dec2binC", (DL_FUNC) &dec2binC, 3},
  {"bin2decC", (DL_FUNC) &bin2decC, 3},
  {NULL, NULL, 0}
};


static const R_CallMethodDef callMethods[] = {
  // { "your_R_function_name", (DL_FUNC) &your_C_function_name, number_of_arguments },
  {"simulate_async_R", (DL_FUNC) &simulate_async_R,  14},
  {"simulate_sync_R", (DL_FUNC) &simulate_sync_R,  13},
  {"simulate_async_return_states_R", (DL_FUNC) &simulate_async_return_states_R, 12},
  { NULL, NULL, 0 }
};



void R_init_PARBONET(DllInfo *info) {
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
