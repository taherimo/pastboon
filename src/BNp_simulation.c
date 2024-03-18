#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NUL
#include <stdbool.h>
#include "common.h"
#include "boolean_network.h"
#include "random.h"


void state_transition_BNp_synchronous(unsigned int * currentState, BooleanNetworkWithPerturbations * net, unsigned int elementsPerEntry)
{

  unsigned int i = 0, k = 0;

  unsigned int nextState[elementsPerEntry];

  for (i = 0; i < elementsPerEntry; ++i)
    nextState[i] = 0;

  for (i = 1; i <= net->num_nodes; ++i)
  {

    int previous_node_state = GET_BIT(currentState[(i - 1) / BITS_PER_BLOCK_32],
                                      (i - 1) % BITS_PER_BLOCK_32);

    if (doublerand_1() > net->p[i-1]) {

      if (net->fixed_nodes[i-1] == -1)
        // the gene is not fixed
      {


        unsigned long long inputdec = 0;

        for (k = net->input_positions[i-1]; k < net->input_positions[i]; k++)
        {
          if (net->inputs[k])
            // if the input of the function is not 0 (constant gene), take input bit
          {
            unsigned int gene = net->inputs[k] - 1;
            unsigned int bit = (GET_BIT(currentState[net->non_fixed_node_bits[gene] / BITS_PER_BLOCK_32],
                           net->non_fixed_node_bits[gene] % BITS_PER_BLOCK_32));

            inputdec |= bit	<< (net->input_positions[i] - k - 1);
          }
        }
        // determine transition function
        int transition = net->outputs[net->output_positions[i-1] + inputdec];

        if(transition != -1) {

          nextState[(i - 1) / BITS_PER_BLOCK_32] |= ((unsigned int)transition << ((i - 1) % BITS_PER_BLOCK_32));
        }
        else
          // this is a dummy function for a constant gene
          // => value does not change
          nextState[(i - 1) / BITS_PER_BLOCK_32] |= (previous_node_state << ((i - 1) % BITS_PER_BLOCK_32));

        //(GET_BIT(currentState[(i-1) / BITS_PER_BLOCK_32],												             // (i-1) % BITS_PER_BLOCK_32) << (idx % BITS_PER_BLOCK_32));




      }
      else {

        //int transition = previous_node_state;
        int transition = net->fixed_nodes[i-1];


        nextState[(i - 1) / BITS_PER_BLOCK_32] |= ((unsigned int)transition << ((i - 1) % BITS_PER_BLOCK_32));
      }
    }

    else {
      int transition = 1 - previous_node_state;
      nextState[(i - 1) / BITS_PER_BLOCK_32] |= ((unsigned int)transition << ((i - 1) % BITS_PER_BLOCK_32));

    }



  }

  //printf("stateTransition %u %u %d\n", currentState[0], nextState[0], elementsPerEntry);
  memcpy(currentState,&nextState,sizeof(unsigned int) * elementsPerEntry);
}


double ** get_node_activities_BNp_sync_traj(BooleanNetworkWithPerturbations * net, double * initial_prob, unsigned int num_repeats, int num_steps, unsigned int num_elements) {
  double* traj_vals = CALLOC(net->num_nodes * (num_steps+1), sizeof(double));
  double** traj = CALLOC(net->num_nodes, sizeof(double*));

  double c = 1.0 / num_repeats;

  unsigned int current_state[num_elements];


  for (unsigned int i=0;i<net->num_nodes;i++){
    //traj[i] = (unsigned int *)malloc(net->numElements*sizeof(int));
    traj[i] = traj_vals + i*(num_steps+1);
  }

  unsigned int i = 0, j = 0, k = 0;

  // for(i = 0; i < net->num_nodes; i++) {
  //   for(j = 0; j < num_steps; j++) {
  //     printf("traj[%u][%u]= %.00f\n", i, j, traj[i*num_steps  + j]);
  //   }
  // }


  for (i = 0; i < num_repeats; i++) {
    //unsigned int currentState[net->numElements];

    for(j = 0; j < num_elements; j++) {
      //current_state[j] = (rand() << (BITS_PER_BLOCK_32 - 1)) ^ rand();
      //printf("current state in for block %d:  %u\n", j, current_state[j]);
      current_state[j] = 0;
    }

    for(k = 0; k < net->num_nodes; k++) {
      if (initial_prob == NULL) {
        if (doublerand_1() < 0.5) {
          current_state[k / BITS_PER_BLOCK_32] |= (1 << (k % BITS_PER_BLOCK_32));
        }
      }
      else if(initial_prob[k] > 0 & initial_prob[k] < 1) {
        //double r = doublerand_1();
        if (doublerand_1() < initial_prob[k]) {
          //printf("%d ----- %f ----- %f ------ %f\n", k, net->initial_prob[k], r, unif_rand());
          current_state[k / BITS_PER_BLOCK_32] |= (1 << (k % BITS_PER_BLOCK_32));
        }
      }
      else { // initial state probability is 0 or 1
        current_state[k / BITS_PER_BLOCK_32] |= (((unsigned int)initial_prob[k]) << (k % BITS_PER_BLOCK_32));
      }
      //printf("inside loop\n");
      //printf("Bit %d = %d\n", k, GET_BIT(current_state[k / BITS_PER_BLOCK_32], k % BITS_PER_BLOCK_32));
      //printf("block = %lu, bit = %lu, value = %d\n", k / BITS_PER_BLOCK_32, k % BITS_PER_BLOCK_32, GET_BIT(current_state[k / BITS_PER_BLOCK_32], k % BITS_PER_BLOCK_32));
      if(GET_BIT(current_state[k / BITS_PER_BLOCK_32], k % BITS_PER_BLOCK_32)) {
        //if(GET_BIT_ARRAY(current_state, k)) {
        traj[k][0] += c;
      }
      //traj[k * num_steps] += c;
      //printf("traj[%d][0] = %f\n", k, traj[k][0]);
    }


    for (j = 1; j <= num_steps; j++)
    {

      state_transition_BNp_synchronous(current_state, net, num_elements);
      //stateTransition(current_state,net,num_elements);

      //printf("current state in for block 0 in step %d:  %u\n", j, current_state[0]);

      for(k = 0; k < net->num_nodes; k++) {
        //printf("inside loop\n");
        if(GET_BIT(current_state[k / BITS_PER_BLOCK_32], k % BITS_PER_BLOCK_32)) {
          // printf("yes, %f\n", c);
          //traj[i*num_steps  + j] += c;
          //printf("traj before = %0.0000f\n", traj[i][j]);
          traj[k][j] += c;
          //traj[k * num_steps + j] += c;
          //printf("traj after = %0.0000f\n", traj[i][j]);
        }
      }

    }


    // for(k = 0; k < net->num_nodes; k++) {
    //   printf("traj[%d][0] = %f\n", k, traj[k][0]);
    // }

    //printf("\n");


    //printf("%d %u %u\n", i, currentStates[i][0], currentStates[i][1]);
    //currentStates[i]=stateTransition(currentStates[i], net);

  }


  return traj;
}



double * get_node_activities_BNp_sync_last_step(BooleanNetworkWithPerturbations * net, double * initial_prob, unsigned int num_repeats, int num_steps, unsigned int num_elements) {

  double* traj = CALLOC(net->num_nodes, sizeof(double));

  double c = 1.0 / num_repeats;

  unsigned int current_state[num_elements];

  unsigned int i = 0, j = 0, k = 0;


  for (i = 0; i < num_repeats; i++) {

    for(j = 0; j < num_elements; j++) {
      current_state[j] = 0;
    }

    for(k = 0; k < net->num_nodes; k++) {
      if(initial_prob[k] > 0 & initial_prob[k] < 1) {
        if (doublerand_1() <= initial_prob[k]) {
          current_state[k / BITS_PER_BLOCK_32] |= (1 << (k % BITS_PER_BLOCK_32));
        }
      }
      else { // initial state probability is 0 or 1
        current_state[k / BITS_PER_BLOCK_32] |= (((unsigned int)initial_prob[k]) << (k % BITS_PER_BLOCK_32));
      }

    }

    for (j = 1; j <= num_steps; j++)
    {
      state_transition_BNp_synchronous(current_state, net, num_elements);
    }

    for(k = 0; k < net->num_nodes; k++) {
      if(GET_BIT(current_state[k / BITS_PER_BLOCK_32], k % BITS_PER_BLOCK_32)) {
        traj[k] += c;
      }
    }

  }

  return traj;
}


unsigned int ** get_reached_states_BNp_sync_batch(BooleanNetworkWithPerturbations * net, unsigned int * initial_states, unsigned int num_initial_states, int num_steps, unsigned int num_elements)

{

  unsigned int i = 0, j = 0;

  printf("num_initial_states=%d, num_steps=%d, num_elements=%u\n",num_initial_states,num_steps,num_elements);

  unsigned int * reached_states_vals = CALLOC(num_initial_states * num_elements, sizeof(unsigned int));
  unsigned int ** reached_states = CALLOC(num_initial_states, sizeof(int*));

  for (i=0;i<num_initial_states;i++){
    //traj[i] = (unsigned int *)malloc(net->numElements*sizeof(int));
    reached_states[i] = reached_states_vals + i*num_elements;
  }


  if(initial_states==NULL) {
    initial_states = CALLOC(num_initial_states * num_elements, sizeof(unsigned int));
    for (i=0;i<num_initial_states;i++){
      for(j=0;j<num_elements;j++) {
        initial_states[i*num_elements + j] = uintrand();
      }
    }
  }



  unsigned int current_state[num_elements];

  for (i = 0; i < num_initial_states; i++) {

    for(j = 0; j < num_elements; j++) {
      current_state[j] = initial_states[i * num_elements + j];
      printf("initial_states[%u]=%u\n",j,initial_states[i * num_elements + j]);
    }

    for (j = 1; j <= num_steps; j++)
    {

      state_transition_BNp_synchronous(current_state, net, num_elements);
      //stateTransition(current_state,net,num_elements);
      //printf("current state in for block 0 in step %d:  %u\n", j, current_state[0]);

    }

    for(j = 0; j < num_elements; j++) {
      reached_states[i][j] = current_state[j];
      printf("reached_states[%u][%u]=%u\n",i,j,reached_states[i][j]);
    }

  }


  return(reached_states);


}

unsigned int ** get_reached_states_BNp_sync_single(BooleanNetworkWithPerturbations * net, unsigned int * initial_state, unsigned int num_repeats, int num_steps, unsigned int num_elements)

{

  unsigned int i = 0, j = 0;

  printf("num_initial_states=%d, num_steps=%d, num_elements=%u\n",num_repeats,num_steps,num_elements);

  unsigned int * reached_states_vals = CALLOC(num_repeats * num_elements, sizeof(unsigned int));
  unsigned int ** reached_states = CALLOC(num_repeats, sizeof(int*));

  for (i=0;i<num_repeats;i++){
    //traj[i] = (unsigned int *)malloc(net->numElements*sizeof(int));
    reached_states[i] = reached_states_vals + i*num_elements;
  }


  if(initial_state==NULL) {
    for(i=0;i<num_elements;i++) {
      initial_state[i] = uintrand();
    }
  }



  unsigned int current_state[num_elements];

  for (i = 0; i < num_repeats; i++) {

    for(j = 0; j < num_elements; j++) {
      current_state[j] = initial_state[j];
      printf("initial_states[%u]=%u\n",j,initial_state[j]);
    }

    for (j = 1; j <= num_steps; j++)
    {

      state_transition_BNp_synchronous(current_state, net, num_elements);
      //stateTransition(current_state,net,num_elements);
      //printf("current state in for block 0 in step %d:  %u\n", j, current_state[0]);

    }

    for(j = 0; j < num_elements; j++) {
      reached_states[i][j] = current_state[j];
      printf("reached_states[%u][%u]=%u\n",i,j,reached_states[i][j]);
    }

  }


  return(reached_states);


}


SEXP get_reached_states_BNp_sync_single_R(SEXP inputs, SEXP input_positions,
                               SEXP outputs, SEXP output_positions,
                               SEXP fixed_nodes, SEXP p, SEXP initial_state,
                               SEXP repeats, SEXP steps) {



  BooleanNetworkWithPerturbations network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p = REAL(p);



  unsigned int numNonFixed = 0, i;
  for (i = 0; i < network.num_nodes; i++)
  {
    if (network.fixed_nodes[i] == -1)
    {
      network.non_fixed_node_bits[i] = numNonFixed++;
    }
  }


  unsigned int * _initial_state = NULL;
  if (!isNull(initial_state) && length(initial_state) > 0)
    _initial_state = (unsigned int *) INTEGER(initial_state);


  unsigned int _numElements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _numElements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _numElements = network.num_nodes / BITS_PER_BLOCK_32 + 1;



  unsigned int _num_repeats = (unsigned int) *INTEGER(repeats);
  unsigned int _num_steps = (unsigned int) *INTEGER(steps);


  //srand(INTEGER(seed)[0]);


  GetRNGstate();  // Activate R's random number generator



  unsigned int ** reached_states = get_reached_states_BNp_sync_single(&network, _initial_state, _num_repeats, _num_steps, _numElements);


  for(unsigned int j = 0; j < _numElements; j++) {
    printf("reached_states[%u][%u]=%u\n",0,j,reached_states[0][j]);
  }


  SEXP result = PROTECT(allocVector(INTSXP, _num_repeats * _numElements));
  //memcpy(&REAL(result)[0], traj, network.num_nodes * (_num_steps + 1) * sizeof(double));

  //memcpy(INTEGER(result), reached_states, _num_initial_states * sizeof(unsigned int));


  for (unsigned int i = 0; i < _num_repeats; ++i) {
    memcpy(&INTEGER(result)[i * _numElements], reached_states[i], _numElements * sizeof(unsigned int));
  }

  PutRNGstate();  // Deactivate R's random number generator


  UNPROTECT(1);

  FREE(network.non_fixed_node_bits);

  return result;


}


SEXP get_reached_states_BNp_sync_batch_R(SEXP inputs, SEXP input_positions,
                                   SEXP outputs, SEXP output_positions,
                                   SEXP fixed_nodes, SEXP p, SEXP initial_states,
                                   SEXP num_initial_states, SEXP steps) {



  BooleanNetworkWithPerturbations network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p = REAL(p);


  unsigned int _num_initial_states = INTEGER(num_initial_states)[0];

  unsigned int * _initial_states = NULL;
  if (!isNull(initial_states) && length(initial_states) > 0)
    _initial_states = (unsigned int *) INTEGER(initial_states);


  unsigned int _numElements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _numElements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _numElements = network.num_nodes / BITS_PER_BLOCK_32 + 1;


  //unsigned int _num_initial_states = length(initial_states) / _numElements;


  //int * _initial_states = INTEGER(initial_states);
  //unsigned int * _initial_states = (unsigned int *) INTEGER(initial_states);



  unsigned int numNonFixed = 0, i;
  for (i = 0; i < network.num_nodes; i++)
  {
    if (network.fixed_nodes[i] == -1)
    {
      network.non_fixed_node_bits[i] = numNonFixed++;
    }
  }



  int _num_steps = *INTEGER(steps);



  //srand(INTEGER(seed)[0]);




  GetRNGstate();  // Activate R's random number generator



  unsigned int ** reached_states =get_reached_states_BNp_sync_batch(&network, _initial_states, _num_initial_states, _num_steps, _numElements);


  for(unsigned int j = 0; j < _numElements; j++) {
    printf("reached_states[%u][%u]=%u\n",0,j,reached_states[0][j]);
  }


  SEXP result = PROTECT(allocVector(INTSXP, _num_initial_states * _numElements));
  //memcpy(&REAL(result)[0], traj, network.num_nodes * (_num_steps + 1) * sizeof(double));

  //memcpy(INTEGER(result), reached_states, _num_initial_states * sizeof(unsigned int));


  for (unsigned int i = 0; i < _num_initial_states; ++i) {
    memcpy(&INTEGER(result)[i * _numElements], reached_states[i], _numElements * sizeof(unsigned int));
  }

  PutRNGstate();  // Deactivate R's random number generator


  UNPROTECT(1);

  FREE(network.non_fixed_node_bits);

  return result;


}


SEXP get_node_activities_sync_R(SEXP inputs, SEXP input_positions,
                     SEXP outputs, SEXP output_positions,
                     SEXP fixed_nodes, SEXP p, SEXP initial_prob,
                     SEXP steps, SEXP repeats, SEXP last_step) {


  BooleanNetworkWithPerturbations network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p = REAL(p);

  double * _initial_prob = NULL;
  if (!isNull(initial_prob) && length(initial_prob) > 0)
    _initial_prob = REAL(initial_prob);


  unsigned int numNonFixed = 0, i;
  for (i = 0; i < network.num_nodes; i++)
  {
    if (network.fixed_nodes[i] == -1)
    {
      network.non_fixed_node_bits[i] = numNonFixed++;
    }
  }


  unsigned int _numElements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _numElements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _numElements = network.num_nodes / BITS_PER_BLOCK_32 + 1;


  //unsigned int* _startStates = (unsigned int*) INTEGER(startStates);
  //unsigned long long * _startStates = (unsigned long long *) INTEGER(startStates);

  //printf("start state in simulate_R: %u\n", _startStates[0]);
  //unsigned int _numStartStates = length(startStates)/_numElements;  // max 32 bits
  //unsigned long long _numStartStates = length(startStates);

  unsigned int _num_steps = *INTEGER(steps);
  unsigned int _numRepeats = *INTEGER(repeats);

  int _last_step = (bool) (*INTEGER(last_step));



  GetRNGstate();  // Activate R's random number generator

  SEXP result;

  if(_last_step) {

    double * traj = get_node_activities_BNp_sync_last_step(&network, _initial_prob, _numRepeats, _num_steps, _numElements);


    result = PROTECT(allocVector(REALSXP, network.num_nodes));
    //memcpy(&REAL(result)[0], traj, network.num_nodes * (_num_steps + 1) * sizeof(double));

    memcpy(REAL(result), traj, network.num_nodes * sizeof(double));

  }
  else {

    double ** traj = get_node_activities_BNp_sync_traj(&network, _initial_prob, _numRepeats, _num_steps, _numElements);
    //unsigned long long * reachedStates = simulate_singleInt(&network, (long long *)_startStates, (long long)_numStartStates, _steps);


    result = PROTECT(allocVector(REALSXP, network.num_nodes * (_num_steps + 1)));
    //memcpy(&REAL(result)[0], traj, network.num_nodes * (_num_steps + 1) * sizeof(double));

    for (unsigned int i = 0; i < network.num_nodes; ++i) {
      memcpy(&REAL(result)[i * (_num_steps+1)], traj[i], (_num_steps + 1) * sizeof(double));
    }


    //memcpy(INTEGER(result), reachedStates, _numStartStates * network.numElements * sizeof(int));

    //free(reachedStates);

  }

  PutRNGstate();  // Deactivate R's random number generator



  UNPROTECT(1);

  FREE(network.non_fixed_node_bits);

  return result;


}



