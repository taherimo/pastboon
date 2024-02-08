simulate_bn_return_states <- function(net, p00, p01, p10, p11,
                                      steps, asynchronous=T, initial_states=NULL,
                                      num_random_initial_states=1000, update_prob=NULL) {


  # assert ncols(initial_states) == length(net$genes)

  if(is.null(initial_states)) {
    initial_states_bin <- sample(0:1, num_random_initial_states * length(bn$genes), rep = T)
    if(num_random_initial_states>1) {
      initial_states_bin <- matrix(initial_states_bin, nrow = num_random_initial_states, byrow = TRUE)
      colnames(initial_states_bin) <- net$genes
    } else {
      names(initial_states_bin) <- net$genes
    }
  }

  if(is.vector(initial_states_bin)) {
    initial_states_dec <- bin2dec(initial_states_bin, length(net$genes))
  }
  else {
    initial_states_dec <- as.vector(apply(initial_states, 1, bin2dec, len=length(net$genes)))
    print(t(apply(apply(initial_states_bin, 1, bin2dec, len=length(net$genes)),2,dec2bin,len=length(net$genes))))
  }

  # the C code requires all interactions to be coded into one vector:
  # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
  inputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$input)))
  input_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$input)))))

  # Do the same for the transition functions.
  outputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$func)))
  output_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$func)))))

  if(asynchronous) {

    #reached_states <- .Call("simulate_async_return_states_R", inputs, input_positions,
    #      outputs, output_positions, as.integer(net$fixed), p00, p01, p10, p11,
    #      as.integer(initial_states), update_prob, as.integer(steps))

    reached_states <- .Call("simulate_async_return_states_R", inputs, input_positions,
                            outputs, output_positions, as.integer(net$fixed), p00, p01, p10, p11,
                            as.integer(initial_states_dec), update_prob, as.integer(steps))

  }


  print(reached_states)


  if(is.vector(initial_states_bin)) {

    reached_states_bin <- dec2bin(reached_states,len=length(net$genes))

  }
  else {

    reached_states <- matrix(reached_states, nrow = nrow(initial_states), byrow = TRUE)

    reached_states_bin <- apply(reached_states, 1, dec2bin, len=length(net$genes))


    reached_states_bin <- t(reached_states_bin)[,1:length(net$genes)]
    colnames(reached_states_bin) <- net$genes

  }


  return(reached_states_bin)


}
