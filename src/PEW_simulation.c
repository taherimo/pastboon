/*
 * This code is derived from the BoolNet package.
 * Original files: BoolNet/src/statespace_search.h,
 *                 BoolNet/src/statespace_search.c,
 *                 BoolNet/src/attractor_search_interface.c
 */

#include "boolean_network.h"
#include "helper_functions.h"
#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <stdlib.h>

static inline void apply_single_function_PEW(unsigned int *currentState,
                                             unsigned int geneIdx,
                                             ProbabilisticEdgeWeight *net) {
  unsigned int k = 0;
  unsigned int gene, bit = 0;

  int previous_node_state = GET_BIT(currentState[geneIdx / BITS_PER_BLOCK_32],
                                    geneIdx % BITS_PER_BLOCK_32);

  if (net->fixed_nodes[geneIdx] == -1)
  // the gene is not fixed
  {
    unsigned long long inputdec = 0;

    for (k = net->input_positions[geneIdx];
         k < net->input_positions[geneIdx + 1]; k++) {
      if (net->inputs[k])
      // if the input of the function is not 0 (constant gene), take input bit
      {
        gene = net->inputs[k] - 1;
        bit = (GET_BIT(currentState[gene / BITS_PER_BLOCK_32],
                       gene % BITS_PER_BLOCK_32));

        if (bit == 0) {
          if (doublerand_1() > net->p_off[gene])
            bit = 1;
        } else {
          if (doublerand_1() > net->p_on[gene])
            bit = 0;
        }

        inputdec |= bit << (net->input_positions[geneIdx + 1] - k - 1);
      }
    }
    // determine transition function
    int transition = net->outputs[net->output_positions[geneIdx] + inputdec];

    currentState[geneIdx / BITS_PER_BLOCK_32] = CLEAR_BIT(
        currentState[geneIdx / BITS_PER_BLOCK_32], geneIdx % BITS_PER_BLOCK_32);

    if (transition != -1) {
      // apply transition function
      currentState[geneIdx / BITS_PER_BLOCK_32] |=
          (transition << (geneIdx % BITS_PER_BLOCK_32));

    } else
      // this is a dummy function for a constant gene
      // => value does not change
      currentState[geneIdx / BITS_PER_BLOCK_32] |=
          (previous_node_state << (geneIdx % BITS_PER_BLOCK_32));

  } else {

    // int transition = previous_node_state;
    int transition = net->fixed_nodes[geneIdx];

    currentState[geneIdx / BITS_PER_BLOCK_32] = CLEAR_BIT(
        currentState[geneIdx / BITS_PER_BLOCK_32], geneIdx % BITS_PER_BLOCK_32);

    currentState[geneIdx / BITS_PER_BLOCK_32] |=
        (transition << (geneIdx % BITS_PER_BLOCK_32));
  }
}

static inline void
state_transition_PEW_asynchronous(unsigned int *currentState,
                                  double *update_prob,
                                  ProbabilisticEdgeWeight *net) {
  unsigned int i;

  if (update_prob == NULL)
  // uniform gene selection
  {
    unsigned int r;

    r = intrand(net->num_nodes);

    // make a transition with the chosen gene
    apply_single_function_PEW(currentState, r, net);
  } else {
    double r = doublerand_1();

    // find the last index in the cumulative distribution that
    // is less than <r>
    for (i = 0; i < net->num_nodes; ++i) {
      if ((update_prob[i] < r) && (update_prob[i + 1] >= r))
        break;
    }
    // make a transition with the chosen gene
    apply_single_function_PEW(currentState, i, net);
  }
}

void state_transition_PEW_synchronous(unsigned int *currentState,
                                      ProbabilisticEdgeWeight *net,
                                      unsigned int elementsPerEntry) {

  unsigned int i = 0, k = 0;
  unsigned int gene, bit = 0;

  unsigned int nextState[elementsPerEntry];

  for (i = 0; i < elementsPerEntry; ++i)
    nextState[i] = 0;

  for (i = 1; i <= net->num_nodes; ++i) {

    int previous_node_state = GET_BIT(currentState[(i - 1) / BITS_PER_BLOCK_32],
                                      (i - 1) % BITS_PER_BLOCK_32);

    if (net->fixed_nodes[i - 1] == -1)
    // the gene is not fixed
    {

      unsigned long long inputdec = 0;

      for (k = net->input_positions[i - 1]; k < net->input_positions[i]; k++) {
        if (net->inputs[k])
        // if the input of the function is not 0 (constant gene), take input bit
        {
          gene = net->inputs[k] - 1;
          bit = (GET_BIT(currentState[gene / BITS_PER_BLOCK_32],
                         gene % BITS_PER_BLOCK_32));

          // bit = (GET_BIT(currentState[net->non_fixed_node_bits[gene] / BITS_PER_BLOCK_32],
          //                net->non_fixed_node_bits[gene] % BITS_PER_BLOCK_32));


          if (bit == 0) {
            if (doublerand_1() > net->p_off[gene]) // net->p_off[gene]
              bit = 1;
          } else {
            if (doublerand_1() > net->p_on[gene])
              bit = 0;
          }

          inputdec |= bit << (net->input_positions[i] - k - 1);
        }
      }
      // determine transition function
      int transition = net->outputs[net->output_positions[i - 1] + inputdec];

      if (transition != -1) {

        nextState[(i - 1) / BITS_PER_BLOCK_32] |=
            ((unsigned int)transition << ((i - 1) % BITS_PER_BLOCK_32));
      } else
        // this is a dummy function for a constant gene
        // => value does not change
        nextState[(i - 1) / BITS_PER_BLOCK_32] |=
            (previous_node_state << ((i - 1) % BITS_PER_BLOCK_32));

      //(GET_BIT(currentState[(i-1) / BITS_PER_BLOCK_32],
      //// (i-1) % BITS_PER_BLOCK_32) << (idx % BITS_PER_BLOCK_32));

    } else {

      // int transition = previous_node_state;
      int transition = net->fixed_nodes[i - 1];

      nextState[(i - 1) / BITS_PER_BLOCK_32] |=
          ((unsigned int)transition << ((i - 1) % BITS_PER_BLOCK_32));
    }
  }

  memcpy(currentState, &nextState, sizeof(unsigned int) * elementsPerEntry);
}


double **get_node_activities_PEW_async_traj(
    ProbabilisticEdgeWeight *net, double *update_prob, double *initial_prob,
    unsigned int num_repeats, int num_steps, unsigned int num_elements)

{

  // https://stackoverflow.com/questions/13761988/what-happens-to-memory-allocated-by-c-functions-in-r-language



  double *traj_vals = CALLOC(net->num_nodes * (num_steps + 1), sizeof(double));
  double **traj = CALLOC(net->num_nodes, sizeof(double *));

  double c = 1.0 / num_repeats;

  unsigned int current_state[num_elements];

  unsigned int i = 0, j = 0, k = 0;

  for (i = 0; i < net->num_nodes; i++) {
    traj[i] = traj_vals + i * (num_steps + 1);
  }

  for (i = 0; i < num_repeats; i++) {

    for (j = 0; j < num_elements; j++) {
      current_state[j] = 0;
    }

    for (k = 0; k < net->num_nodes; k++) {
      if (initial_prob == NULL) {
        if (doublerand_1() < 0.5) {
          current_state[k / BITS_PER_BLOCK_32] |=
              (1 << (k % BITS_PER_BLOCK_32));
        }
      } else if ((initial_prob[k] > 0) && (initial_prob[k] < 1)) {
        if (doublerand_1() < initial_prob[k]) {
          current_state[k / BITS_PER_BLOCK_32] |=
              (1 << (k % BITS_PER_BLOCK_32));
        }
      } else { // initial state probability is 0 or 1
        current_state[k / BITS_PER_BLOCK_32] |=
            (((unsigned int)initial_prob[k]) << (k % BITS_PER_BLOCK_32));
      }

      if (GET_BIT(current_state[k / BITS_PER_BLOCK_32],
                  k % BITS_PER_BLOCK_32)) {
        traj[k][0] += c;
      }
    }

    for (j = 1; j <= num_steps; j++) {

      state_transition_PEW_asynchronous(current_state, update_prob, net);

      for (k = 0; k < net->num_nodes; k++) {

        if (GET_BIT(current_state[k / BITS_PER_BLOCK_32],
                    k % BITS_PER_BLOCK_32)) {
          traj[k][j] += c;
        }
      }
    }
  }

  return traj;
}

double **get_node_activities_PEW_sync_traj(ProbabilisticEdgeWeight *net,
                                           double *initial_prob,
                                           unsigned int num_repeats,
                                           int num_steps,
                                           unsigned int num_elements) {
  double *traj_vals = CALLOC(net->num_nodes * (num_steps + 1), sizeof(double));
  double **traj = CALLOC(net->num_nodes, sizeof(double *));

  double c = 1.0 / num_repeats;

  unsigned int current_state[num_elements];

  for (unsigned int i = 0; i < net->num_nodes; i++) {
    traj[i] = traj_vals + i * (num_steps + 1);
  }

  unsigned int i = 0, j = 0, k = 0;

  for (i = 0; i < num_repeats; i++) {

    for (j = 0; j < num_elements; j++) {
      current_state[j] = 0;
    }

    for (k = 0; k < net->num_nodes; k++) {
      if (initial_prob == NULL) {
        if (doublerand_1() < 0.5) {
          current_state[k / BITS_PER_BLOCK_32] |=
            (1 << (k % BITS_PER_BLOCK_32));
        }
      } else if ((initial_prob[k] > 0) && (initial_prob[k] < 1)) {
        if (doublerand_1() < initial_prob[k]) {
          current_state[k / BITS_PER_BLOCK_32] |=
            (1 << (k % BITS_PER_BLOCK_32));
        }
      } else { // initial state probability is 0 or 1
        current_state[k / BITS_PER_BLOCK_32] |=
          (((unsigned int)initial_prob[k]) << (k % BITS_PER_BLOCK_32));
      }

      if (GET_BIT(current_state[k / BITS_PER_BLOCK_32],
                  k % BITS_PER_BLOCK_32)) {
        // if(GET_BIT_ARRAY(current_state, k)) {
        traj[k][0] += c;
      }
    }

    for (j = 1; j <= num_steps; j++) {

      state_transition_PEW_synchronous(current_state, net, num_elements);

      for (k = 0; k < net->num_nodes; k++) {
        if (GET_BIT(current_state[k / BITS_PER_BLOCK_32],
                    k % BITS_PER_BLOCK_32)) {

          traj[k][j] += c;
        }
      }
    }
  }

  return traj;
}

double *get_node_activities_PEW_async_last_step(
    ProbabilisticEdgeWeight *net, double *update_prob, double *initial_prob,
    unsigned int num_repeats, int num_steps, unsigned int num_elements)

{

  // https://stackoverflow.com/questions/13761988/what-happens-to-memory-allocated-by-c-functions-in-r-language



  double *traj = CALLOC(net->num_nodes, sizeof(double));

  double c = 1.0 / num_repeats;

  unsigned int current_state[num_elements];

  unsigned int i = 0, j = 0, k = 0;

  for (i = 0; i < num_repeats; i++) {

    for (j = 0; j < num_elements; j++) {
      current_state[j] = 0;
    }

    for (k = 0; k < net->num_nodes; k++) {
      if (initial_prob == NULL) {
        if (doublerand_1() < 0.5) {
          current_state[k / BITS_PER_BLOCK_32] |=
            (1 << (k % BITS_PER_BLOCK_32));
        }
      } else if ((initial_prob[k] > 0) && (initial_prob[k] < 1)) {
        if (doublerand_1() < initial_prob[k]) {
          current_state[k / BITS_PER_BLOCK_32] |=
            (1 << (k % BITS_PER_BLOCK_32));
        }
      } else { // initial state probability is 0 or 1
        current_state[k / BITS_PER_BLOCK_32] |=
          (((unsigned int)initial_prob[k]) << (k % BITS_PER_BLOCK_32));
      }
    }

    for (j = 0; j < num_steps; j++) {
      state_transition_PEW_asynchronous(current_state, update_prob, net);
    }

    for (k = 0; k < net->num_nodes; k++) {
      if (GET_BIT(current_state[k / BITS_PER_BLOCK_32],
                  k % BITS_PER_BLOCK_32)) {
        traj[k] += c;
      }
    }
  }

  return traj;
}

double *get_node_activities_PEW_sync_last_step(ProbabilisticEdgeWeight *net,
                                               double *initial_prob,
                                               unsigned int num_repeats,
                                               int num_steps,
                                               unsigned int num_elements) {

  double *traj = CALLOC(net->num_nodes, sizeof(double));

  double c = 1.0 / num_repeats;

  unsigned int current_state[num_elements];

  unsigned int i = 0, j = 0, k = 0;

  for (i = 0; i < num_repeats; i++) {

    for (j = 0; j < num_elements; j++) {
      current_state[j] = 0;
    }

    for (k = 0; k < net->num_nodes; k++) {
      if (initial_prob == NULL) {
        if (doublerand_1() < 0.5) {
          current_state[k / BITS_PER_BLOCK_32] |=
            (1 << (k % BITS_PER_BLOCK_32));
        }
      } else if ((initial_prob[k] > 0) && (initial_prob[k] < 1)) {
        if (doublerand_1() < initial_prob[k]) {
          current_state[k / BITS_PER_BLOCK_32] |=
            (1 << (k % BITS_PER_BLOCK_32));
        }
      } else { // initial state probability is 0 or 1
        current_state[k / BITS_PER_BLOCK_32] |=
          (((unsigned int)initial_prob[k]) << (k % BITS_PER_BLOCK_32));
      }
    }

    for (j = 0; j < num_steps; j++) {
      state_transition_PEW_synchronous(current_state, net, num_elements);
    }

    for (k = 0; k < net->num_nodes; k++) {
      if (GET_BIT(current_state[k / BITS_PER_BLOCK_32],
                  k % BITS_PER_BLOCK_32)) {
        traj[k] += c;
      }
    }
  }

  return traj;
}



double **get_pairwise_transitions_PEW_async(
    ProbabilisticEdgeWeight *net, double *update_prob, unsigned int **states,
    unsigned int num_states, unsigned int num_repeats, int num_steps,
    unsigned int num_elements) {

  double *trans_mat_vals = CALLOC(num_states * num_states, sizeof(double));
  double **trans_mat = CALLOC(num_states, sizeof(double *));

  unsigned int i, j, k, l = 0;

  // double c = 1.0 / num_repeats;
  double c = 1.0;

  for (i = 0; i < num_states; i++) {
    trans_mat[i] = trans_mat_vals + i * num_states;
  }

  unsigned int current_state[num_elements];

  for (i = 0; i < num_states; i++) {

    for (j = 0; j < num_repeats; j++) {

      for (k = 0; k < num_elements; k++) {
        current_state[k] = states[i][k];
      }

      for (k = 1; k <= num_steps; k++) {

        state_transition_PEW_asynchronous(current_state, update_prob, net);

        for (l = 0; l < num_states; l++) {

          if (areArraysEqual(current_state, states[l], num_elements)) {
            trans_mat[i][l] += c;
            break;
          }
        }
      }
    }
  }

  return (trans_mat);
}

double **get_pairwise_transitions_PEW_sync(ProbabilisticEdgeWeight *net,
                                           unsigned int **states,
                                           unsigned int num_states,
                                           unsigned int num_repeats,
                                           int num_steps,
                                           unsigned int num_elements) {

  double *trans_mat_vals = CALLOC(num_states * num_states, sizeof(double));
  double **trans_mat = CALLOC(num_states, sizeof(double *));

  unsigned int i, j, k, l = 0;

  // double c = 1.0 / num_repeats;
  double c = 1.0;

  for (i = 0; i < num_states; i++) {
    trans_mat[i] = trans_mat_vals + i * num_states;
  }

  unsigned int current_state[num_elements];

  for (i = 0; i < num_states; i++) {

    for (j = 0; j < num_repeats; j++) {

      for (k = 0; k < num_elements; k++) {
        current_state[k] = states[i][k];
      }

      for (k = 1; k <= num_steps; k++) {

        state_transition_PEW_synchronous(current_state, net, num_elements);

        for (l = 0; l < num_states; l++) {

          // if(areArraysEqual(reached_states[j], states[k], num_elements)) {
          if (areArraysEqual(current_state, states[l], num_elements)) {
            trans_mat[i][l] += c;
            break;
          }
        }
      }
    }
  }

  return (trans_mat);
}


unsigned int **get_reached_states_PEW_async_batch(
    ProbabilisticEdgeWeight *net, double *update_prob,
    unsigned int *initial_states, unsigned int num_initial_states,
    int num_steps, unsigned int num_elements)

{

  unsigned int i = 0, j = 0;

  unsigned int *reached_states_vals =
    CALLOC(num_initial_states * num_elements, sizeof(unsigned int));
  unsigned int **reached_states = CALLOC(num_initial_states, sizeof(int *));

  for (i = 0; i < num_initial_states; i++) {
    reached_states[i] = reached_states_vals + i * num_elements;
  }

  if (initial_states == NULL) {
    initial_states =
      CALLOC(num_initial_states * num_elements, sizeof(unsigned int));
    for (i = 0; i < num_initial_states; i++) {
      for (j = 0; j < num_elements; j++) {
        initial_states[i * num_elements + j] = intrand_fullrange();
      }
    }
  }

  unsigned int current_state[num_elements];

  for (i = 0; i < num_initial_states; i++) {

    for (j = 0; j < num_elements; j++) {
      current_state[j] = initial_states[i * num_elements + j];
    }

    for (j = 1; j <= num_steps; j++) {

      state_transition_PEW_asynchronous(current_state, update_prob, net);
    }

    for (j = 0; j < num_elements; j++) {
      reached_states[i][j] = current_state[j];
    }
  }

  return (reached_states);
}



unsigned int **get_reached_states_PEW_sync_batch(
    ProbabilisticEdgeWeight *net, unsigned int *initial_states,
    unsigned int num_initial_states, int num_steps, unsigned int num_elements)

{

  unsigned int i = 0, j = 0;

  unsigned int *reached_states_vals =
    CALLOC(num_initial_states * num_elements, sizeof(unsigned int));
  unsigned int **reached_states = CALLOC(num_initial_states, sizeof(int *));

  for (i = 0; i < num_initial_states; i++) {
    reached_states[i] = reached_states_vals + i * num_elements;
  }

  if (initial_states == NULL) {
    initial_states =
      CALLOC(num_initial_states * num_elements, sizeof(unsigned int));
    for (i = 0; i < num_initial_states; i++) {
      for (j = 0; j < num_elements; j++) {
        initial_states[i * num_elements + j] = intrand_fullrange();
      }
    }
  }

  unsigned int current_state[num_elements];

  for (i = 0; i < num_initial_states; i++) {

    for (j = 0; j < num_elements; j++) {
      current_state[j] = initial_states[i * num_elements + j];
    }

    for (j = 1; j <= num_steps; j++) {

      state_transition_PEW_synchronous(current_state, net, num_elements);
    }

    for (j = 0; j < num_elements; j++) {
      reached_states[i][j] = current_state[j];
    }
  }

  return (reached_states);
}

unsigned int **get_reached_states_PEW_async_single(ProbabilisticEdgeWeight *net,
                                                   double *update_prob,
                                                   unsigned int *initial_state,
                                                   unsigned int num_repeats,
                                                   int num_steps,
                                                   unsigned int num_elements)

{

  unsigned int i = 0, j = 0;

  unsigned int *reached_states_vals =
    CALLOC(num_repeats * num_elements, sizeof(unsigned int));
  unsigned int **reached_states = CALLOC(num_repeats, sizeof(int *));

  for (i = 0; i < num_repeats; i++) {
    reached_states[i] = reached_states_vals + i * num_elements;
  }

  if (initial_state == NULL) {
    initial_state = CALLOC(num_elements, sizeof(unsigned int));
    for (i = 0; i < num_elements; i++) {
      initial_state[i] = intrand_fullrange();
    }
  }

  unsigned int current_state[num_elements];

  for (i = 0; i < num_repeats; i++) {

    for (j = 0; j < num_elements; j++) {
      current_state[j] = initial_state[j];
    }

    for (j = 1; j <= num_steps; j++) {
      state_transition_PEW_asynchronous(current_state, update_prob, net);
    }

    for (j = 0; j < num_elements; j++) {
      reached_states[i][j] = current_state[j];
    }
  }

  return (reached_states);
}

unsigned int **get_reached_states_PEW_sync_single(ProbabilisticEdgeWeight *net,
                                                  unsigned int *initial_state,
                                                  unsigned int num_repeats,
                                                  int num_steps,
                                                  unsigned int num_elements)

{

  unsigned int i = 0, j = 0;

  unsigned int *reached_states_vals =
    CALLOC(num_repeats * num_elements, sizeof(unsigned int));
  unsigned int **reached_states = CALLOC(num_repeats, sizeof(int *));

  for (i = 0; i < num_repeats; i++) {
    reached_states[i] = reached_states_vals + i * num_elements;
  }

  if (initial_state == NULL) {
    initial_state = CALLOC(num_elements, sizeof(unsigned int));
    for (i = 0; i < num_elements; i++) {
      initial_state[i] = intrand_fullrange();
    }
  }

  unsigned int current_state[num_elements];

  for (i = 0; i < num_repeats; i++) {

    for (j = 0; j < num_elements; j++) {
      current_state[j] = initial_state[j];
    }

    for (j = 1; j <= num_steps; j++) {
      state_transition_PEW_synchronous(current_state, net, num_elements);
    }

    for (j = 0; j < num_elements; j++) {
      reached_states[i][j] = current_state[j];
    }
  }

  return (reached_states);
}

SEXP get_node_activities_PEW_async_R(SEXP inputs, SEXP input_positions,
                                     SEXP outputs, SEXP output_positions,
                                     SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                     SEXP initial_prob, SEXP update_prob,
                                     SEXP steps, SEXP repeats, SEXP last_step) {

  ProbabilisticEdgeWeight network;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  //network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  double *_update_prob = NULL;
  if ((!isNull(update_prob)) && (length(update_prob) > 0))
    _update_prob = REAL(update_prob);

  double *_initial_prob = NULL;
  if ((!isNull(initial_prob)) && (length(initial_prob) > 0))
    _initial_prob = REAL(initial_prob);

  // unsigned int num_non_fixed = 0, i;
  // for (i = 0; i < network.num_nodes; i++) {
  //   if (network.fixed_nodes[i] == -1) {
  //     network.non_fixed_node_bits[i] = num_non_fixed++;
  //   }
  // }

  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;

  unsigned int _num_steps = (unsigned int)*INTEGER(steps);
  unsigned int _num_repeats = (unsigned int)*INTEGER(repeats);

  int _last_step = (bool)(*INTEGER(last_step));



  GetRNGstate();

  SEXP result;

  if (_last_step) {

    double *traj = get_node_activities_PEW_async_last_step(
      &network, _update_prob, _initial_prob, _num_repeats, _num_steps,
      _num_elements);

    result = PROTECT(allocVector(REALSXP, network.num_nodes));

    memcpy(REAL(result), traj, network.num_nodes * sizeof(double));

    FREE(traj);

  } else {

    double **traj = get_node_activities_PEW_async_traj(
      &network, _update_prob, _initial_prob, _num_repeats, _num_steps,
      _num_elements);

    result =
      PROTECT(allocVector(REALSXP, network.num_nodes * (_num_steps + 1)));

    for (unsigned int i = 0; i < network.num_nodes; ++i) {
      memcpy(&REAL(result)[i * (_num_steps + 1)], traj[i],
             (_num_steps + 1) * sizeof(double));
    }

    FREE(traj);

  }

  PutRNGstate();

  UNPROTECT(1);

  //FREE(network.non_fixed_node_bits);

  return result;
}

SEXP get_node_activities_PEW_sync_R(SEXP inputs, SEXP input_positions,
                                    SEXP outputs, SEXP output_positions,
                                    SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                    SEXP initial_prob, SEXP steps, SEXP repeats,
                                    SEXP last_step) {

  ProbabilisticEdgeWeight network;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  //network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  double *_initial_prob = NULL;
  if ((!isNull(initial_prob)) && (length(initial_prob) > 0))
    _initial_prob = REAL(initial_prob);

  // unsigned int num_non_fixed = 0, i;
  // for (i = 0; i < network.num_nodes; i++) {
  //   if (network.fixed_nodes[i] == -1) {
  //     network.non_fixed_node_bits[i] = num_non_fixed++;
  //   }
  // }

  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;

  unsigned int _num_steps = *INTEGER(steps);
  unsigned int _num_repeats = *INTEGER(repeats);

  int _last_step = (bool)(*INTEGER(last_step));

  GetRNGstate();

  SEXP result;

  if (_last_step) {

    double *traj = get_node_activities_PEW_sync_last_step(
      &network, _initial_prob, _num_repeats, _num_steps, _num_elements);

    result = PROTECT(allocVector(REALSXP, network.num_nodes));

    memcpy(REAL(result), traj, network.num_nodes * sizeof(double));

    FREE(traj);

  } else {

    double **traj = get_node_activities_PEW_sync_traj(
      &network, _initial_prob, _num_repeats, _num_steps, _num_elements);

    result =
    PROTECT(allocVector(REALSXP, network.num_nodes * (_num_steps + 1)));

    for (unsigned int i = 0; i < network.num_nodes; ++i) {
      memcpy(&REAL(result)[i * (_num_steps + 1)], traj[i],
             (_num_steps + 1) * sizeof(double));
    }

    FREE(traj);
  }

  PutRNGstate();

  UNPROTECT(1);

  //FREE(network.non_fixed_node_bits);

  return result;
}

SEXP get_pairwise_transitions_PEW_async_R(SEXP inputs, SEXP input_positions,
                                          SEXP outputs, SEXP output_positions,
                                          SEXP fixed_nodes, SEXP p_on,
                                          SEXP p_off, SEXP update_prob,
                                          SEXP states, SEXP num_states,
                                          SEXP steps, SEXP repeats) {

  ProbabilisticEdgeWeight network;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  //network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  double *_update_prob = NULL;
  if ((!isNull(update_prob)) && (length(update_prob) > 0))
    _update_prob = REAL(update_prob);

  // unsigned int num_non_fixed = 0, i;
  // for (i = 0; i < network.num_nodes; i++) {
  //   if (network.fixed_nodes[i] == -1) {
  //     network.non_fixed_node_bits[i] = num_non_fixed++;
  //   }
  // }

  unsigned int i, j = 0;

  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;

  unsigned int *_states = (unsigned int *)INTEGER(states);

  unsigned int _num_states = (unsigned int)*INTEGER(num_states);

  unsigned int _num_steps = (unsigned int)*INTEGER(steps);
  unsigned int _num_repeats = (unsigned int)*INTEGER(repeats);

  unsigned int *_states_2d_vals =
      CALLOC(_num_states * _num_elements, sizeof(double));
  unsigned int **_states_2d = CALLOC(_num_states, sizeof(unsigned int *));

  for (i = 0; i < _num_states; i++) {
    _states_2d[i] = _states_2d_vals + i * _num_elements;
  }

  for (i = 0; i < _num_states; i++) {
    for (j = 0; j < _num_elements; j++) {
      _states_2d[i][j] = _states[i * _num_elements + j];
    }
  }

  GetRNGstate();

  double **transition_matrix = get_pairwise_transitions_PEW_async(
      &network, _update_prob, _states_2d, _num_states, _num_repeats, _num_steps,
      _num_elements);

  SEXP result = PROTECT(allocVector(REALSXP, _num_states * _num_states));

  for (unsigned int i = 0; i < _num_states; ++i) {
    memcpy(&REAL(result)[i * _num_states], transition_matrix[i],
           _num_states * sizeof(double));
  }

  PutRNGstate();

  UNPROTECT(1);

  //FREE(network.non_fixed_node_bits);
  FREE(transition_matrix);
  FREE(_states_2d);

  return result;
}

SEXP get_pairwise_transitions_PEW_sync_R(SEXP inputs, SEXP input_positions,
                                         SEXP outputs, SEXP output_positions,
                                         SEXP fixed_nodes, SEXP p_on,
                                         SEXP p_off, SEXP states,
                                         SEXP num_states, SEXP steps,
                                         SEXP repeats) {

  ProbabilisticEdgeWeight network;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  //network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  // unsigned int num_non_fixed = 0, i;
  // for (i = 0; i < network.num_nodes; i++) {
  //   if (network.fixed_nodes[i] == -1) {
  //     network.non_fixed_node_bits[i] = num_non_fixed++;
  //   }
  // }

  unsigned int i, j = 0;

  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;

  unsigned int *_states = (unsigned int *)INTEGER(states);

  unsigned int _num_states = (unsigned int)*INTEGER(num_states);

  unsigned int _num_steps = (unsigned int)*INTEGER(steps);
  unsigned int _num_repeats = (unsigned int)*INTEGER(repeats);

  unsigned int *_states_2d_vals =
      CALLOC(_num_states * _num_elements, sizeof(double));
  unsigned int **_states_2d = CALLOC(_num_states, sizeof(unsigned int *));

  for (i = 0; i < _num_states; i++) {
    _states_2d[i] = _states_2d_vals + i * _num_elements;
  }

  for (i = 0; i < _num_states; i++) {
    for (j = 0; j < _num_elements; j++) {
      _states_2d[i][j] = _states[i * _num_elements + j];
    }
  }

  GetRNGstate();

  double **transition_matrix = get_pairwise_transitions_PEW_sync(
      &network, _states_2d, _num_states, _num_repeats, _num_steps,
      _num_elements);

  SEXP result = PROTECT(allocVector(REALSXP, _num_states * _num_states));

  for (unsigned int i = 0; i < _num_states; ++i) {
    memcpy(&REAL(result)[i * _num_states], transition_matrix[i],
           _num_states * sizeof(double));
  }

  PutRNGstate();

  UNPROTECT(1);

  //FREE(network.non_fixed_node_bits);
  FREE(transition_matrix);
  FREE(_states_2d);

  return result;
}



SEXP get_reached_states_PEW_async_batch_R(SEXP inputs, SEXP input_positions,
                                          SEXP outputs, SEXP output_positions,
                                          SEXP fixed_nodes, SEXP p_on,
                                          SEXP p_off, SEXP initial_states,
                                          SEXP num_initial_states,
                                          SEXP update_prob, SEXP steps) {

  ProbabilisticEdgeWeight network;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  //network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  double *_update_prob = NULL;
  if ((!isNull(update_prob)) && (length(update_prob) > 0))
    _update_prob = REAL(update_prob);

  unsigned int _num_initial_states = INTEGER(num_initial_states)[0];

  unsigned int *_initial_states = NULL;
  if (!isNull(initial_states) && length(initial_states) > 0)
    _initial_states = (unsigned int *)INTEGER(initial_states);

  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;

  // unsigned int num_non_fixed = 0, i;
  // for (i = 0; i < network.num_nodes; i++) {
  //   if (network.fixed_nodes[i] == -1) {
  //     network.non_fixed_node_bits[i] = num_non_fixed++;
  //   }
  // }

  int _num_steps = *INTEGER(steps);



  GetRNGstate();

  unsigned int **reached_states = get_reached_states_PEW_async_batch(
      &network, _update_prob, _initial_states, _num_initial_states, _num_steps,
      _num_elements);

  SEXP result =
      PROTECT(allocVector(INTSXP, _num_initial_states * _num_elements));

  for (unsigned int i = 0; i < _num_initial_states; ++i) {
    memcpy(&INTEGER(result)[i * _num_elements], reached_states[i],
           _num_elements * sizeof(unsigned int));
  }

  PutRNGstate();

  UNPROTECT(1);

  //FREE(network.non_fixed_node_bits);
  FREE(reached_states);

  return result;
}

SEXP get_reached_states_PEW_sync_batch_R(SEXP inputs, SEXP input_positions,
                                         SEXP outputs, SEXP output_positions,
                                         SEXP fixed_nodes, SEXP p_on,
                                         SEXP p_off, SEXP initial_states,
                                         SEXP num_initial_states, SEXP steps) {

  ProbabilisticEdgeWeight network;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  //network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  unsigned int _num_initial_states = INTEGER(num_initial_states)[0];

  unsigned int *_initial_states = NULL;
  if ((!isNull(initial_states)) && (length(initial_states) > 0))
    _initial_states = (unsigned int *)INTEGER(initial_states);

  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;

  // unsigned int num_non_fixed = 0, i;
  // for (i = 0; i < network.num_nodes; i++) {
  //   if (network.fixed_nodes[i] == -1) {
  //     network.non_fixed_node_bits[i] = num_non_fixed++;
  //   }
  // }

  int _num_steps = *INTEGER(steps);

  GetRNGstate();

  unsigned int **reached_states = get_reached_states_PEW_sync_batch(
      &network, _initial_states, _num_initial_states, _num_steps, _num_elements);

  SEXP result =
      PROTECT(allocVector(INTSXP, _num_initial_states * _num_elements));

  for (unsigned int i = 0; i < _num_initial_states; ++i) {
    memcpy(&INTEGER(result)[i * _num_elements], reached_states[i],
           _num_elements * sizeof(unsigned int));
  }

  PutRNGstate();

  UNPROTECT(1);

  //FREE(network.non_fixed_node_bits);
  FREE(reached_states);

  return result;
}

SEXP get_reached_states_PEW_async_single_R(SEXP inputs, SEXP input_positions,
                                           SEXP outputs, SEXP output_positions,
                                           SEXP fixed_nodes, SEXP p_on,
                                           SEXP p_off, SEXP update_prob,
                                           SEXP initial_state, SEXP repeats,
                                           SEXP steps) {

  ProbabilisticEdgeWeight network;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  //network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  double *_update_prob = NULL;
  if ((!isNull(update_prob)) && (length(update_prob) > 0))
    _update_prob = REAL(update_prob);

  // unsigned int num_non_fixed = 0, i;
  // for (i = 0; i < network.num_nodes; i++) {
  //   if (network.fixed_nodes[i] == -1) {
  //     network.non_fixed_node_bits[i] = num_non_fixed++;
  //   }
  // }

  unsigned int *_initial_state = NULL;
  if ((!isNull(initial_state)) && (length(initial_state) > 0))
    _initial_state = (unsigned int *)INTEGER(initial_state);

  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;

  unsigned int _num_repeats = (unsigned int)*INTEGER(repeats);
  unsigned int _num_steps = (unsigned int)*INTEGER(steps);



  GetRNGstate();

  unsigned int **reached_states = get_reached_states_PEW_async_single(
    &network, _update_prob, _initial_state, _num_repeats, _num_steps,
    _num_elements);


  SEXP result = PROTECT(allocVector(INTSXP, _num_repeats * _num_elements));

  for (unsigned int i = 0; i < _num_repeats; ++i) {
    memcpy(&INTEGER(result)[i * _num_elements], reached_states[i],
           _num_elements * sizeof(unsigned int));
  }

  PutRNGstate();

  UNPROTECT(1);

  //FREE(network.non_fixed_node_bits);
  FREE(reached_states);

  return result;
}

SEXP get_reached_states_PEW_sync_single_R(SEXP inputs, SEXP input_positions,
                                          SEXP outputs, SEXP output_positions,
                                          SEXP fixed_nodes, SEXP p_on,
                                          SEXP p_off, SEXP initial_state,
                                          SEXP repeats, SEXP steps) {

  ProbabilisticEdgeWeight network;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  //network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  // unsigned int num_non_fixed = 0, i;
  // for (i = 0; i < network.num_nodes; i++) {
  //   if (network.fixed_nodes[i] == -1) {
  //     network.non_fixed_node_bits[i] = num_non_fixed++;
  //   }
  // }

  unsigned int *_initial_state = NULL;
  if ((!isNull(initial_state)) && (length(initial_state) > 0))
    _initial_state = (unsigned int *)INTEGER(initial_state);

  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;

  unsigned int _num_repeats = (unsigned int)*INTEGER(repeats);
  unsigned int _num_steps = (unsigned int)*INTEGER(steps);



  GetRNGstate();

  unsigned int **reached_states = get_reached_states_PEW_sync_single(
    &network, _initial_state, _num_repeats, _num_steps, _num_elements);

  SEXP result = PROTECT(allocVector(INTSXP, _num_repeats * _num_elements));

  for (unsigned int i = 0; i < _num_repeats; ++i) {
    memcpy(&INTEGER(result)[i * _num_elements], reached_states[i],
           _num_elements * sizeof(unsigned int));
  }

  PutRNGstate();

  UNPROTECT(1);

  //FREE(network.non_fixed_node_bits);
  FREE(reached_states);

  return result;
}
