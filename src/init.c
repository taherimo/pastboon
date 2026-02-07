<<<<<<< HEAD
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

extern SEXP get_reached_states_SDDS_async_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_SDDS_sync_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_SDDS_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_SDDS_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_SDDS_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_SDDS_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


extern SEXP get_node_activities_BNp_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_node_activities_BNp_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_BNp_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_BNp_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_BNp_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_BNp_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_BNp_async_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_BNp_async_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_BNp_sync_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_BNp_sync_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


extern SEXP get_node_activities_PEW_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_node_activities_PEW_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_PEW_async_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_PEW_sync_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_PEW_async_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_PEW_sync_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_PEW_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_PEW_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_PEW_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_PEW_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


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
  {"get_reached_states_SDDS_async_single_R", (DL_FUNC) &get_reached_states_SDDS_async_single_R, 13},
  {"get_reached_states_SDDS_sync_single_R", (DL_FUNC) &get_reached_states_SDDS_sync_single_R, 12},
  {"count_pairwise_transitions_SDDS_async_R", (DL_FUNC) &count_pairwise_transitions_SDDS_async_R, 14},
  {"count_pairwise_transitions_SDDS_sync_R", (DL_FUNC) &count_pairwise_transitions_SDDS_sync_R, 13},
  {"calc_pairwise_reachability_SDDS_async_R", (DL_FUNC) &calc_pairwise_reachability_SDDS_async_R, 14},
  {"calc_pairwise_reachability_SDDS_sync_R", (DL_FUNC) &calc_pairwise_reachability_SDDS_sync_R, 13},
  {"get_node_activities_BNp_async_R", (DL_FUNC) &get_node_activities_BNp_async_R, 11},
  {"get_node_activities_BNp_sync_R", (DL_FUNC) &get_node_activities_BNp_sync_R, 10},
  {"count_pairwise_transitions_BNp_async_R", (DL_FUNC) &count_pairwise_transitions_BNp_async_R, 11},
  {"count_pairwise_transitions_BNp_sync_R", (DL_FUNC) &count_pairwise_transitions_BNp_sync_R, 10},
  {"calc_pairwise_reachability_BNp_async_R", (DL_FUNC) &calc_pairwise_reachability_BNp_async_R, 11},
  {"calc_pairwise_reachability_BNp_sync_R", (DL_FUNC) &calc_pairwise_reachability_BNp_sync_R, 10},
  {"get_reached_states_BNp_async_single_R", (DL_FUNC) &get_reached_states_BNp_async_single_R, 10},
  {"get_reached_states_BNp_async_batch_R", (DL_FUNC) &get_reached_states_BNp_async_batch_R, 10},
  {"get_reached_states_BNp_sync_single_R", (DL_FUNC) &get_reached_states_BNp_sync_single_R, 9},
  {"get_reached_states_BNp_sync_batch_R", (DL_FUNC) &get_reached_states_BNp_sync_batch_R, 9},
  {"get_node_activities_PEW_async_R", (DL_FUNC) &get_node_activities_PEW_async_R, 12},
  {"get_node_activities_PEW_sync_R", (DL_FUNC) &get_node_activities_PEW_sync_R, 11},
  {"get_reached_states_PEW_async_batch_R", (DL_FUNC) &get_reached_states_PEW_async_batch_R, 11},
  {"get_reached_states_PEW_sync_batch_R", (DL_FUNC) &get_reached_states_PEW_sync_batch_R, 10},
  {"get_reached_states_PEW_async_single_R", (DL_FUNC) &get_reached_states_PEW_async_single_R, 11},
  {"get_reached_states_PEW_sync_single_R", (DL_FUNC) &get_reached_states_PEW_async_single_R, 10},
  {"count_pairwise_transitions_PEW_async_R", (DL_FUNC) &count_pairwise_transitions_PEW_async_R, 12},
  {"count_pairwise_transitions_PEW_sync_R", (DL_FUNC) &count_pairwise_transitions_PEW_sync_R, 11},
  {"calc_pairwise_reachability_PEW_async_R", (DL_FUNC) &calc_pairwise_reachability_PEW_async_R, 12},
  {"calc_pairwise_reachability_PEW_sync_R", (DL_FUNC) &calc_pairwise_reachability_PEW_sync_R, 11},
  { NULL, NULL, 0 }
};



void R_init_pastboon(DllInfo *info) {
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
=======
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

extern SEXP get_reached_states_SDDS_async_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_SDDS_sync_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_SDDS_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_SDDS_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_SDDS_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_SDDS_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


extern SEXP get_node_activities_BNp_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_node_activities_BNp_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_BNp_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_BNp_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_BNp_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_BNp_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_BNp_async_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_BNp_async_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_BNp_sync_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_BNp_sync_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


extern SEXP get_node_activities_PEW_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_node_activities_PEW_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_PEW_async_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_PEW_sync_batch_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_PEW_async_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP get_reached_states_PEW_sync_single_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_PEW_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP count_pairwise_transitions_PEW_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_PEW_async_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP calc_pairwise_reachability_PEW_sync_R(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


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
  {"get_reached_states_SDDS_async_single_R", (DL_FUNC) &get_reached_states_SDDS_async_single_R, 13},
  {"get_reached_states_SDDS_sync_single_R", (DL_FUNC) &get_reached_states_SDDS_sync_single_R, 12},
  {"count_pairwise_transitions_SDDS_async_R", (DL_FUNC) &count_pairwise_transitions_SDDS_async_R, 14},
  {"count_pairwise_transitions_SDDS_sync_R", (DL_FUNC) &count_pairwise_transitions_SDDS_sync_R, 13},
  {"calc_pairwise_reachability_SDDS_async_R", (DL_FUNC) &calc_pairwise_reachability_SDDS_async_R, 14},
  {"calc_pairwise_reachability_SDDS_sync_R", (DL_FUNC) &calc_pairwise_reachability_SDDS_sync_R, 13},
  {"get_node_activities_BNp_async_R", (DL_FUNC) &get_node_activities_BNp_async_R, 11},
  {"get_node_activities_BNp_sync_R", (DL_FUNC) &get_node_activities_BNp_sync_R, 10},
  {"count_pairwise_transitions_BNp_async_R", (DL_FUNC) &count_pairwise_transitions_BNp_async_R, 11},
  {"count_pairwise_transitions_BNp_sync_R", (DL_FUNC) &count_pairwise_transitions_BNp_sync_R, 10},
  {"calc_pairwise_reachability_BNp_async_R", (DL_FUNC) &calc_pairwise_reachability_BNp_async_R, 11},
  {"calc_pairwise_reachability_BNp_sync_R", (DL_FUNC) &calc_pairwise_reachability_BNp_sync_R, 10},
  {"get_reached_states_BNp_async_single_R", (DL_FUNC) &get_reached_states_BNp_async_single_R, 10},
  {"get_reached_states_BNp_async_batch_R", (DL_FUNC) &get_reached_states_BNp_async_batch_R, 10},
  {"get_reached_states_BNp_sync_single_R", (DL_FUNC) &get_reached_states_BNp_sync_single_R, 9},
  {"get_reached_states_BNp_sync_batch_R", (DL_FUNC) &get_reached_states_BNp_sync_batch_R, 9},
  {"get_node_activities_PEW_async_R", (DL_FUNC) &get_node_activities_PEW_async_R, 12},
  {"get_node_activities_PEW_sync_R", (DL_FUNC) &get_node_activities_PEW_sync_R, 11},
  {"get_reached_states_PEW_async_batch_R", (DL_FUNC) &get_reached_states_PEW_async_batch_R, 11},
  {"get_reached_states_PEW_sync_batch_R", (DL_FUNC) &get_reached_states_PEW_sync_batch_R, 10},
  {"get_reached_states_PEW_async_single_R", (DL_FUNC) &get_reached_states_PEW_async_single_R, 11},
  {"get_reached_states_PEW_sync_single_R", (DL_FUNC) &get_reached_states_PEW_async_single_R, 10},
  {"count_pairwise_transitions_PEW_async_R", (DL_FUNC) &count_pairwise_transitions_PEW_async_R, 12},
  {"count_pairwise_transitions_PEW_sync_R", (DL_FUNC) &count_pairwise_transitions_PEW_sync_R, 11},
  {"calc_pairwise_reachability_PEW_async_R", (DL_FUNC) &calc_pairwise_reachability_PEW_async_R, 12},
  {"calc_pairwise_reachability_PEW_sync_R", (DL_FUNC) &calc_pairwise_reachability_PEW_sync_R, 11},
  { NULL, NULL, 0 }
};



void R_init_pastboon(DllInfo *info) {
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
>>>>>>> b1a170373c94b8c09d5cd4181a69638be70554b2
