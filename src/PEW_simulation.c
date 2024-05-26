#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NUL
#include <stdbool.h>
#include "common.h"
#include "boolean_network.h"
#include "random.h"


static inline void apply_single_function_PEW(unsigned int * currentState, unsigned int geneIdx, ProbabilisticEdgeWeight * net)
{
  unsigned int k = 0;
  unsigned int gene, bit = 0;


  int previous_node_state = GET_BIT(currentState[geneIdx / BITS_PER_BLOCK_32],
                                    geneIdx % BITS_PER_BLOCK_32);


  if (net->fixed_nodes[geneIdx] == -1)
    // the gene is not fixed
  {
    unsigned long long inputdec = 0;

    for (k = net->input_positions[geneIdx]; k < net->input_positions[geneIdx+1]; k++)
    {
      if (net->inputs[k])
        // if the input of the function is not 0 (constant gene), take input bit
      {
        gene = net->inputs[k] - 1;
        bit = (GET_BIT(currentState[gene / BITS_PER_BLOCK_32], gene % BITS_PER_BLOCK_32));

        if(bit==0) {
          if(doublerand_1() > net->p_off[gene])
            bit = 1;
        }
        else {
          if(doublerand_1() > net->p_on[gene])
            bit = 0;
        }

        inputdec |= bit	<< (net->input_positions[geneIdx+1] - k - 1);
      }
    }
    // determine transition function
    int transition = net->outputs[net->output_positions[geneIdx] + inputdec];

    currentState[geneIdx / BITS_PER_BLOCK_32] = CLEAR_BIT(currentState[geneIdx / BITS_PER_BLOCK_32],
                                                          geneIdx % BITS_PER_BLOCK_32);


    if(transition != -1) {
      // apply transition function
      currentState[geneIdx / BITS_PER_BLOCK_32] |= (transition << (geneIdx % BITS_PER_BLOCK_32));

    }
    else
      // this is a dummy function for a constant gene
      // => value does not change
      currentState[geneIdx / BITS_PER_BLOCK_32] |= (previous_node_state << (geneIdx % BITS_PER_BLOCK_32));

    // if(geneIdx==80) {
    //   printf("previous state = %d\n", previous_node_state);
    //   printf("transition = %d\n", transition);
    //   printf("input position start = %d\n", net->input_positions[geneIdx]);
    //   printf("input position end = %d\n", net->input_positions[geneIdx+1]);
    // }

  }
  else {


    //int transition = previous_node_state;
    int transition = net->fixed_nodes[geneIdx];


    currentState[geneIdx / BITS_PER_BLOCK_32] = CLEAR_BIT(currentState[geneIdx / BITS_PER_BLOCK_32],
                                                          geneIdx % BITS_PER_BLOCK_32);

    currentState[geneIdx / BITS_PER_BLOCK_32] |= (transition << (geneIdx % BITS_PER_BLOCK_32));



  }


}


static inline void state_transition_PEW_asynchronous(unsigned int * currentState, double* update_prob, ProbabilisticEdgeWeight * net)
{
  unsigned int i;

  if (update_prob == NULL)
    // uniform gene selection
  {
    unsigned int r;

    r = intrand(net->num_nodes);

    // make a transition with the chosen gene
    apply_single_function_PEW(currentState,r,net);
  }
  else
  {
    double r = doublerand_1();

    // find the last index in the cumulative distribution that
    // is less than <r>
    for (i = 0; i < net->num_nodes; ++i)
    {
      if (update_prob[i] < r && update_prob[i+1] >= r)
        break;
    }
    // make a transition with the chosen gene
    apply_single_function_PEW(currentState,i,net);
  }
}

void state_transition_PEW_synchronous(unsigned int * currentState, ProbabilisticEdgeWeight * net, unsigned int elementsPerEntry)
{

  unsigned int i = 0, k = 0;
  unsigned int gene, bit = 0;

  unsigned int nextState[elementsPerEntry];

  for (i = 0; i < elementsPerEntry; ++i)
    nextState[i] = 0;

  for (i = 1; i <= net->num_nodes; ++i)
  {

    int previous_node_state = GET_BIT(currentState[(i - 1) / BITS_PER_BLOCK_32],
                                      (i - 1) % BITS_PER_BLOCK_32);

    if (net->fixed_nodes[i-1] == -1)
      // the gene is not fixed
    {


      unsigned long long inputdec = 0;

      for (k = net->input_positions[i-1]; k < net->input_positions[i]; k++)
      {
        if (net->inputs[k])
          // if the input of the function is not 0 (constant gene), take input bit
        {
          gene = net->inputs[k] - 1;
          bit = (GET_BIT(currentState[net->non_fixed_node_bits[gene] / BITS_PER_BLOCK_32],
                         net->non_fixed_node_bits[gene] % BITS_PER_BLOCK_32));

          if(bit==0) {
            if(doublerand_1() > net->p_off[gene]) // net->p_off[gene]
              bit = 1;
          }
          else {
            if(doublerand_1() > net->p_on[gene])
              bit = 0;
          }

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

  //printf("stateTransition %u %u %d\n", currentState[0], nextState[0], elementsPerEntry);
  memcpy(currentState,&nextState,sizeof(unsigned int) * elementsPerEntry);
}


double * get_node_activities_PEW_async_last_step(ProbabilisticEdgeWeight * net, double * update_prob, double * initial_prob, unsigned int num_repeats, int num_steps, unsigned int num_elements)

{

  // https://stackoverflow.com/questions/13761988/what-happens-to-memory-allocated-by-c-functions-in-r-language

  // double * traj = CALLOC((num_steps + 1) * net->num_nodes, sizeof(double));

  double* traj = CALLOC(net->num_nodes, sizeof(double));

  double c = 1.0 / num_repeats;

  unsigned int current_state[num_elements];



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
      if(initial_prob[k] > 0 & initial_prob[k] < 1) {
        //double r = doublerand_1();
        if (doublerand_1() <= initial_prob[k]) {
          //printf("%d ----- %f ----- %f ------ %f\n", k, net->initial_prob[k], r, unif_rand());
          current_state[k / BITS_PER_BLOCK_32] |= (1 << (k % BITS_PER_BLOCK_32));
        }
      }
      else { // initial state probability is 0 or 1
        current_state[k / BITS_PER_BLOCK_32] |= (((unsigned int)initial_prob[k]) << (k % BITS_PER_BLOCK_32));
      }

    }


    for (j = 1; j <= num_steps; j++)
    {

      state_transition_PEW_asynchronous(current_state, update_prob, net);
      //stateTransition(current_state,net,num_elements);
      //printf("current state in for block 0 in step %d:  %u\n", j, current_state[0]);
    }

    for(k = 0; k < net->num_nodes; k++) {
      //printf("inside loop\n");
      if(GET_BIT(current_state[k / BITS_PER_BLOCK_32], k % BITS_PER_BLOCK_32)) {
        // printf("yes, %f\n", c);
        //traj[i*num_steps  + j] += c;
        //printf("traj before = %0.0000f\n", traj[i][j]);
        traj[k] += c;
        //traj[k * num_steps + j] += c;
        //printf("traj after = %0.0000f\n", traj[i][j]);
      }

    }


  }


  return traj;

}


double ** get_node_activities_PEW_async_traj(ProbabilisticEdgeWeight * net, double * update_prob, double * initial_prob, unsigned int num_repeats, int num_steps, unsigned int num_elements)

{

  // https://stackoverflow.com/questions/13761988/what-happens-to-memory-allocated-by-c-functions-in-r-language

  // double * traj = CALLOC((num_steps + 1) * net->num_nodes, sizeof(double));

  double* traj_vals = CALLOC(net->num_nodes * (num_steps+1), sizeof(double));
  double** traj = CALLOC(net->num_nodes, sizeof(double*));

  double c = 1.0 / num_repeats;

  unsigned int current_state[num_elements];

  unsigned int i = 0, j = 0, k = 0;


  for (i=0;i<net->num_nodes;i++){
    //traj[i] = (unsigned int *)malloc(net->numElements*sizeof(int));
    traj[i] = traj_vals + i*(num_steps+1);
  }


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

      state_transition_PEW_asynchronous(current_state, update_prob, net);
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

double ** get_node_activities_PEW_sync_traj(ProbabilisticEdgeWeight * net, double * initial_prob, unsigned int num_repeats, int num_steps, unsigned int num_elements) {
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

      state_transition_PEW_synchronous(current_state, net, num_elements);
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


double * get_node_activities_PEW_sync_last_step(ProbabilisticEdgeWeight * net, double * initial_prob, unsigned int num_repeats, int num_steps, unsigned int num_elements) {

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
      state_transition_PEW_synchronous(current_state, net, num_elements);
    }

    for(k = 0; k < net->num_nodes; k++) {
      if(GET_BIT(current_state[k / BITS_PER_BLOCK_32], k % BITS_PER_BLOCK_32)) {
        traj[k] += c;
      }
    }

  }

  return traj;
}

unsigned int ** get_reached_states_PEW_async_batch(ProbabilisticEdgeWeight * net, double* update_prob, unsigned int * initial_states, unsigned int num_initial_states, int num_steps, unsigned int num_elements)

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
      //printf("initial_states[%u]=%u\n",j,initial_states[i * num_elements + j]);
    }

    for (j = 1; j <= num_steps; j++)
    {

      state_transition_PEW_asynchronous(current_state, update_prob, net);
      //stateTransition(current_state,net,num_elements);
      //printf("current state in for block 0 in step %d:  %u\n", j, current_state[0]);

    }

    for(j = 0; j < num_elements; j++) {
      reached_states[i][j] = current_state[j];
      //printf("reached_states[%u][%u]=%u\n",i,j,reached_states[i][j]);
    }

  }


  return(reached_states);


}

unsigned int ** get_reached_states_PEW_async_single(ProbabilisticEdgeWeight * net, double * update_prob, unsigned int * initial_state, unsigned int num_repeats, int num_steps, unsigned int num_elements)

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
    initial_state = CALLOC(num_elements, sizeof(unsigned int));
    for(i=0;i<num_elements;i++) {
      initial_state[i] = uintrand();
    }
  }

  unsigned int current_state[num_elements];


  for (i = 0; i < num_repeats; i++) {

    for(j = 0; j < num_elements; j++) {
      current_state[j] = initial_state[j];
      //printf("initial_states[%u]=%u\n",j,initial_state[j]);
    }

    for (j = 1; j <= num_steps; j++)
    {

      state_transition_PEW_asynchronous(current_state, update_prob, net);
      //stateTransition(current_state,net,num_elements);
      //printf("current state in for block 0 in step %d:  %u\n", j, current_state[0]);

    }

    for(j = 0; j < num_elements; j++) {
      reached_states[i][j] = current_state[j];
      //printf("reached_states[%u][%u]=%u\n",i,j,reached_states[i][j]);
    }

  }


  return(reached_states);


}

unsigned int ** get_reached_states_PEW_sync_batch(ProbabilisticEdgeWeight * net, unsigned int * initial_states, unsigned int num_initial_states, int num_steps, unsigned int num_elements)

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

      state_transition_PEW_synchronous(current_state, net, num_elements);
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


unsigned int ** get_reached_states_PEW_sync_single(ProbabilisticEdgeWeight * net, unsigned int * initial_state, unsigned int num_repeats, int num_steps, unsigned int num_elements)

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
    initial_state = CALLOC(num_elements, sizeof(unsigned int));
    for(i=0;i<num_elements;i++) {
      initial_state[i] = uintrand();
    }
  }



  unsigned int current_state[num_elements];

  for (i = 0; i < num_repeats; i++) {

    for(j = 0; j < num_elements; j++) {
      current_state[j] = initial_state[j];
      //printf("initial_states[%u]=%u\n",j,initial_state[j]);
    }

    for (j = 1; j <= num_steps; j++)
    {

      state_transition_PEW_synchronous(current_state, net, num_elements);
      //stateTransition(current_state,net,num_elements);
      //printf("current state in for block 0 in step %d:  %u\n", j, current_state[0]);

    }

    for(j = 0; j < num_elements; j++) {
      reached_states[i][j] = current_state[j];
      //printf("reached_states[%u][%u]=%u\n",i,j,reached_states[i][j]);
    }

  }


  return(reached_states);


}


double ** get_cumulative_transition_matrix_PEW_async(ProbabilisticEdgeWeight * net, double * update_prob, unsigned int ** states, unsigned int num_states, unsigned int num_repeats, int num_steps, unsigned int num_elements) {

  double* trans_mat_vals = CALLOC(num_states * num_states, sizeof(double));
  double** trans_mat = CALLOC(num_states, sizeof(double*));

  unsigned int i , j, k, l = 0;
  unsigned int ** reached_states;

  //double c = 1.0 / num_repeats;
  double c = 1.0;

  for (i=0;i<num_states;i++){
    //traj[i] = (unsigned int *)malloc(net->numElements*sizeof(int));
    trans_mat[i] = trans_mat_vals + i*num_states;
  }

  unsigned int current_state[num_elements];

  for (i=0;i<num_states;i++){


    //reached_states = get_reached_states_async(net,current_state,num_repeats,num_steps, num_elements);

    for (j=0;j<num_repeats;j++){

      for(k = 0; k < num_elements; k++) {
        current_state[k] = states[i][k];
        //current_state[j] = states[i * num_elements + j];
        //printf("states[%u]=%u\n",i,states[i][j]);
      }

      for (k = 1; k <= num_steps; k++)
      {

        state_transition_PEW_asynchronous(current_state, update_prob, net);
        //stateTransition(current_state,net,num_elements);
        //printf("current state in for block 0 in step %d:  %u\n", j, current_state[0]);

        for(l=0; l<num_states;l++) {

          //if(areArraysEqual(reached_states[j], states[k], num_elements)) {
          if(areArraysEqual(current_state, states[l], num_elements)) {
            trans_mat[i][l] += c;
            //printf("%u --- %u\n", j,k);
            break;
          }
        }

      }

    }
  }

  return(trans_mat);

}

double ** get_cumulative_transition_matrix_PEW_sync(ProbabilisticEdgeWeight * net, unsigned int ** states, unsigned int num_states, unsigned int num_repeats, int num_steps, unsigned int num_elements) {

  double* trans_mat_vals = CALLOC(num_states * num_states, sizeof(double));
  double** trans_mat = CALLOC(num_states, sizeof(double*));

  unsigned int i , j, k, l = 0;
  unsigned int ** reached_states;

  //double c = 1.0 / num_repeats;
  double c = 1.0;

  for (i=0;i<num_states;i++){
    //traj[i] = (unsigned int *)malloc(net->numElements*sizeof(int));
    trans_mat[i] = trans_mat_vals + i*num_states;
  }

  unsigned int current_state[num_elements];

  for (i=0;i<num_states;i++){


    //reached_states = get_reached_states_async(net,current_state,num_repeats,num_steps, num_elements);

    for (j=0;j<num_repeats;j++){

      for(k = 0; k < num_elements; k++) {
        current_state[k] = states[i][k];
        //current_state[j] = states[i * num_elements + j];
        //printf("states[%u]=%u\n",i,states[i][j]);
      }

      for (k = 1; k <= num_steps; k++)
      {

        state_transition_PEW_synchronous(current_state, net, num_elements);
        //stateTransition(current_state,net,num_elements);
        //printf("current state in for block 0 in step %d:  %u\n", j, current_state[0]);


        for(l=0; l<num_states;l++) {

          //if(areArraysEqual(reached_states[j], states[k], num_elements)) {
          if(areArraysEqual(current_state, states[l], num_elements)) {
            trans_mat[i][l] += c;
            //printf("%u --- %u\n", j,k);
            break;
          }
        }

      }


    }
  }

  return(trans_mat);

}

SEXP get_reached_states_PEW_async_single_R(SEXP inputs, SEXP input_positions,
                                            SEXP outputs, SEXP output_positions,
                                            SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                            SEXP update_prob, SEXP initial_state,
                                            SEXP repeats, SEXP steps) {



  ProbabilisticEdgeWeight network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  double * _update_prob = NULL;
  if (!isNull(update_prob) && length(update_prob) > 0)
    _update_prob = REAL(update_prob);


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



  unsigned int ** reached_states = get_reached_states_PEW_async_single(&network, _update_prob, _initial_state, _num_repeats, _num_steps, _numElements);


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


SEXP get_reached_states_PEW_sync_single_R(SEXP inputs, SEXP input_positions,
                                           SEXP outputs, SEXP output_positions,
                                           SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                           SEXP initial_state, SEXP repeats,
                                           SEXP steps) {



  ProbabilisticEdgeWeight network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);


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



  unsigned int ** reached_states = get_reached_states_PEW_sync_single(&network, _initial_state, _num_repeats, _num_steps, _numElements);


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



SEXP get_cumulative_transition_matrix_PEW_async_R(SEXP inputs, SEXP input_positions,
                                                   SEXP outputs, SEXP output_positions,
                                                   SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                                   SEXP update_prob, SEXP states, SEXP num_states,
                                                   SEXP steps, SEXP repeats) {

  ProbabilisticEdgeWeight network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  double * _update_prob = NULL;
  if (!isNull(update_prob) && length(update_prob) > 0)
    _update_prob = REAL(update_prob);

  unsigned int numNonFixed = 0, i;
  for (i = 0; i < network.num_nodes; i++)
  {
    if (network.fixed_nodes[i] == -1)
    {
      network.non_fixed_node_bits[i] = numNonFixed++;
    }
  }

  unsigned int j = 0;


  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;


  unsigned int * _states = (unsigned int *) INTEGER(states);

  unsigned int _num_states = (unsigned int) *INTEGER(num_states);


  unsigned int _num_steps = (unsigned int) *INTEGER(steps);
  unsigned int _num_repeats = (unsigned int) *INTEGER(repeats);

  unsigned int * _states_2d_vals = CALLOC(_num_states * _num_elements, sizeof(double));
  unsigned int ** _states_2d = CALLOC(_num_states, sizeof(unsigned int *));

  for (i=0;i<_num_states;i++){
    //traj[i] = (unsigned int *)malloc(net->numElements*sizeof(int));
    _states_2d[i] = _states_2d_vals + i*_num_elements;
  }

  for(i = 0; i < _num_states; i++) {
    for(j = 0; j < _num_elements; j++) {
      _states_2d[i][j] = _states[i*_num_elements + j];
    }
  }



  GetRNGstate();  // Activate R's random number generator



  double ** transition_matrix = get_cumulative_transition_matrix_PEW_async(&network, _update_prob, _states_2d, _num_states, _num_repeats, _num_steps, _num_elements);


  SEXP result = PROTECT(allocVector(REALSXP, _num_states * _num_states));

  for (unsigned int i = 0; i < _num_states; ++i) {
    memcpy(&REAL(result)[i * _num_states], transition_matrix[i], _num_states * sizeof(double));
  }

  PutRNGstate();  // Deactivate R's random number generator


  UNPROTECT(1);

  FREE(network.non_fixed_node_bits);
  // FREE(_states);
  // FREE(_states_2d);
  // FREE(_states_2d_vals);

  return result;


}


SEXP get_cumulative_transition_matrix_PEW_sync_R(SEXP inputs, SEXP input_positions,
                                                  SEXP outputs, SEXP output_positions,
                                                  SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                                  SEXP states, SEXP num_states, SEXP steps,
                                                  SEXP repeats) {

  ProbabilisticEdgeWeight network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);


  unsigned int numNonFixed = 0, i;
  for (i = 0; i < network.num_nodes; i++)
  {
    if (network.fixed_nodes[i] == -1)
    {
      network.non_fixed_node_bits[i] = numNonFixed++;
    }
  }

  unsigned int j = 0;


  unsigned int _num_elements;

  if (network.num_nodes % BITS_PER_BLOCK_32 == 0)
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32;
  else
    _num_elements = network.num_nodes / BITS_PER_BLOCK_32 + 1;


  unsigned int * _states = (unsigned int *) INTEGER(states);

  unsigned int _num_states = (unsigned int) *INTEGER(num_states);


  unsigned int _num_steps = (unsigned int) *INTEGER(steps);
  unsigned int _num_repeats = (unsigned int) *INTEGER(repeats);

  unsigned int * _states_2d_vals = CALLOC(_num_states * _num_elements, sizeof(double));
  unsigned int ** _states_2d = CALLOC(_num_states, sizeof(unsigned int *));

  for (i=0;i<_num_states;i++){
    //traj[i] = (unsigned int *)malloc(net->numElements*sizeof(int));
    _states_2d[i] = _states_2d_vals + i*_num_elements;
  }

  for(i = 0; i < _num_states; i++) {
    for(j = 0; j < _num_elements; j++) {
      _states_2d[i][j] = _states[i*_num_elements + j];
    }
  }



  GetRNGstate();  // Activate R's random number generator



  double ** transition_matrix = get_cumulative_transition_matrix_PEW_sync(&network, _states_2d, _num_states, _num_repeats, _num_steps, _num_elements);


  SEXP result = PROTECT(allocVector(REALSXP, _num_states * _num_states));

  for (unsigned int i = 0; i < _num_states; ++i) {
    memcpy(&REAL(result)[i * _num_states], transition_matrix[i], _num_states * sizeof(double));
  }

  PutRNGstate();  // Deactivate R's random number generator


  UNPROTECT(1);

  FREE(network.non_fixed_node_bits);
  // FREE(_states);
  // FREE(_states_2d);
  // FREE(_states_2d_vals);

  return result;


}


SEXP get_node_activities_PEW_async_R(SEXP inputs, SEXP input_positions,
                                      SEXP outputs, SEXP output_positions,
                                      SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                      SEXP initial_prob, SEXP update_prob,
                                      SEXP steps, SEXP repeats, SEXP last_step) {
  //int * inputs = INTEGER(inputs_R);
  //SEXP result;
  //result = PROTECT(allocVector(REALSXP, 2));
  //REAL(result)[0] = 123.45;
  //REAL(result)[1] = 67.89;
  //UNPROTECT(1);

  ProbabilisticEdgeWeight network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  double * _update_prob = NULL;
  if (!isNull(update_prob) && length(update_prob) > 0)
    _update_prob = REAL(update_prob);

  double * _initial_prob = NULL;
  if (!isNull(initial_prob) && length(initial_prob) > 0)
    _initial_prob = REAL(initial_prob);
  //network.epsilon = REAL(epsilon)[0];

  // count fixed genes, and create an index array for non-fixed genes:
  // <network.non_fixed_node_bits[i]> contains the bit positions in a state
  // at which the <i>-th gene is stored - this is different from <i>
  // as fixed genes are not stored
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

  unsigned int _num_steps = (unsigned int) *INTEGER(steps);
  unsigned int _numRepeats = (unsigned int) *INTEGER(repeats);

  int _last_step = (bool) (*INTEGER(last_step));

  //srand(INTEGER(seed)[0]);

  GetRNGstate();  // Activate R's random number generator

  SEXP result;

  if(_last_step) {

    double * traj = get_node_activities_PEW_async_last_step(&network, _update_prob, _initial_prob, _numRepeats, _num_steps, _numElements);


    result = PROTECT(allocVector(REALSXP, network.num_nodes));
    //memcpy(&REAL(result)[0], traj, network.num_nodes * (_num_steps + 1) * sizeof(double));

    memcpy(REAL(result), traj, network.num_nodes * sizeof(double));

  }
  else {

    double ** traj = get_node_activities_PEW_async_traj(&network, _update_prob, _initial_prob, _numRepeats, _num_steps, _numElements);


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



SEXP get_node_activities_PEW_sync_R(SEXP inputs, SEXP input_positions,
                                     SEXP outputs, SEXP output_positions,
                                     SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                     SEXP initial_prob, SEXP steps,
                                     SEXP repeats, SEXP last_step) {


  ProbabilisticEdgeWeight network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

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

    double * traj = get_node_activities_PEW_sync_last_step(&network, _initial_prob, _numRepeats, _num_steps, _numElements);


    result = PROTECT(allocVector(REALSXP, network.num_nodes));
    //memcpy(&REAL(result)[0], traj, network.num_nodes * (_num_steps + 1) * sizeof(double));

    memcpy(REAL(result), traj, network.num_nodes * sizeof(double));

  }
  else {

    double ** traj = get_node_activities_PEW_sync_traj(&network, _initial_prob, _numRepeats, _num_steps, _numElements);
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


SEXP get_reached_states_PEW_async_batch_R(SEXP inputs, SEXP input_positions,
                                           SEXP outputs, SEXP output_positions,
                                           SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                           SEXP initial_states, SEXP num_initial_states,
                                           SEXP update_prob, SEXP steps) {



  ProbabilisticEdgeWeight network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

  double * _update_prob = NULL;
  if (!isNull(update_prob) && length(update_prob) > 0)
    _update_prob = REAL(update_prob);


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



  unsigned int ** reached_states = get_reached_states_PEW_async_batch(&network, _update_prob, _initial_states, _num_initial_states, _num_steps, _numElements);


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


SEXP get_reached_states_PEW_sync_batch_R(SEXP inputs, SEXP input_positions,
                                          SEXP outputs, SEXP output_positions,
                                          SEXP fixed_nodes, SEXP p_on, SEXP p_off,
                                          SEXP initial_states, SEXP num_initial_states,
                                          SEXP steps) {



  ProbabilisticEdgeWeight network;
  //network.type = TRUTHTABLE_BOOLEAN_NETWORK;
  network.num_nodes = length(fixed_nodes);
  network.inputs = INTEGER(inputs);
  network.input_positions = INTEGER(input_positions);
  network.outputs = INTEGER(outputs);
  network.output_positions = INTEGER(output_positions);
  network.fixed_nodes = INTEGER(fixed_nodes);
  network.non_fixed_node_bits = CALLOC(network.num_nodes, sizeof(unsigned int));
  network.p_on = REAL(p_on);
  network.p_off = REAL(p_off);

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



  unsigned int ** reached_states = get_reached_states_PEW_sync_batch(&network, _initial_states, _num_initial_states, _num_steps, _numElements);


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



