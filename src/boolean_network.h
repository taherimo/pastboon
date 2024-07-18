#ifndef BOOLEAN_NETWORK_H
#define BOOLEAN_NETWORK_H

/*
 * This code is derived from the BoolNet package.
 * Original file: BoolNet/src/boolean_network.h
 */

typedef struct
{
  unsigned char type;

  unsigned int num_nodes;

  int * fixed_nodes;
  unsigned int * non_fixed_node_bits;

  //int * node_types;

  int * inputs;
  int * input_positions;
  int * outputs;
  int * output_positions;
  double * p;

} BooleanNetworkWithPerturbations;


typedef struct
{
  unsigned char type;

  unsigned int num_nodes;

  int * fixed_nodes;
  unsigned int * non_fixed_node_bits;

  //int * node_types;

  int * inputs;
  int * input_positions;
  int * outputs;
  int * output_positions;

  double * p_on;
  double * p_off;



} ProbabilisticEdgeWeight;


typedef struct
{
  unsigned char type;

  unsigned int num_nodes;

  int * fixed_nodes;
  unsigned int * non_fixed_node_bits;

  //int * node_types;

  int * inputs;
  int * input_positions;
  int * outputs;
  int * output_positions;
  double * p00;
  double * p01;
  double * p10;
  double * p11;

} StochasticDiscreteDynamicalSystem;


#endif
