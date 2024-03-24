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

        if(bit==1) {
          if(doublerand_1() < net->p_off[gene])
            bit = 0;
        }
        else {
          if(doublerand_1() < net->p_on[gene])
            bit = 1;
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

          if(bit==1) {
            if(doublerand_1() < net->p_off[gene])
              bit = 0;
          }
          else {
            if(doublerand_1() < net->p_on[gene])
              bit = 1;
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
