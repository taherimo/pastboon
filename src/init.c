#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */

extern void dec2binC(void *, void *, void *);
extern void bin2decC(void *, void *, void *);


/* .Call calls */

extern SEXP get_node_activities_SDDS_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_node_activities_SDDS_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_SDDS_async_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_SDDS_sync_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_convergence_time_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_convergence_time_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_SDDS_async_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_SDDS_sync_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_transition_matrix_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_transition_matrix_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


extern SEXP get_node_activities_BNp_sync_R();


static const R_CMethodDef cMethods[] = {
  {"dec2binC", (DL_FUNC) &dec2binC, 3},
  {"bin2decC", (DL_FUNC) &bin2decC, 3},
  {NULL, NULL, 0}
};


static const R_CallMethodDef callMethods[] = {
  // { "your_R_function_name", (DL_FUNC) &your_C_function_name, number_of_arguments },
  {"get_node_activities_SDDS_async_R", (DL_FUNC) &get_node_activities_SDDS_async_R,  14},
  {"get_node_activities_SDDS_sync_R", (DL_FUNC) &get_node_activities_SDDS_sync_R,  13},
  {"get_reached_states_SDDS_async_batch_R", (DL_FUNC) &get_reached_states_SDDS_async_batch_R, 13},
  {"get_reached_states_SDDS_sync_batch_R", (DL_FUNC) &get_reached_states_SDDS_sync_batch_R, 12},
  {"get_convergence_time_async_R", (DL_FUNC) &get_convergence_time_async_R, 14},
  {"get_convergence_time_sync_R", (DL_FUNC) &get_convergence_time_sync_R, 13},
  {"get_reached_states_SDDS_async_single_R", (DL_FUNC) &get_reached_states_SDDS_async_single_R, 13},
  {"get_reached_states_SDDS_sync_single_R", (DL_FUNC) &get_reached_states_SDDS_sync_single_R, 12},
  {"get_transition_matrix_async_R", (DL_FUNC) &get_transition_matrix_async_R, 14},
  {"get_transition_matrix_sync_R", (DL_FUNC) &get_transition_matrix_sync_R, 13},
  { NULL, NULL, 0 }
};



void R_init_PARBONET(DllInfo *info) {
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
