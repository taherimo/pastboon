get_reached_states_single <- function(net, p00, p01, p10, p11,
                               repeats, steps=1, initial_state=NULL,
                               update_prob=NULL,asynchronous=T) {


  # Odd behavour when initial_state is a decimal number like 1000


  # the C code requires all interactions to be coded into one vector:
  # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
  inputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$input)))
  input_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$input)))))

  # Do the same for the transition functions.
  outputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$func)))
  output_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$func)))))

  # if (any(initial_state != 0 & initial_state != 1)) {
  #   cat("Error!")
  #   return(NA)
  # }

  # if(is.vector(initial_state)) {
  #   initial_state_dec <- bin2dec(initial_state, length(net$genes))
  # }
  # else {
  #   cat("error!")
  #   return(NA)
  #   #initial_state_dec <- as.vector(apply(initial_state, 1, bin2dec, len=length(net$genes)))
  #   #print(t(apply(apply(initial_states, 1, bin2dec, len=length(net$genes)),2,dec2bin,len=length(net$genes))))
  #
  # }


  initial_state_dec <- bin2dec(initial_state, length(net$genes))


  if(asynchronous) {

    reached_states <- .Call("get_reached_states_SDDS_async_single_R", inputs, input_positions,
                             outputs, output_positions,
                             as.integer(net$fixed),
                             p00, p01, p10, p11,
                             update_prob, as.integer(initial_state_dec),
                             as.integer(repeats), as.integer(steps),
                             PACKAGE = "PARBONET")


  } else {

    reached_states <- .Call("get_reached_states_SDDS_sync_single_R", inputs, input_positions,
                            outputs, output_positions,
                            as.integer(net$fixed),
                            p00, p01, p10, p11,
                            as.integer(initial_state_dec),
                            as.integer(repeats), as.integer(steps),
                            PACKAGE = "PARBONET")




  }


  if(repeats==1) {

    reached_states_bin <- dec2bin(reached_states,len=length(net$genes))

  } else {

    reached_states <- matrix(reached_states, nrow = repeats, byrow = TRUE)

    reached_states_bin <- apply(reached_states, 1, dec2bin, len=length(net$genes))


    reached_states_bin <- t(reached_states_bin)[,1:length(net$genes)]
    colnames(reached_states_bin) <- net$genes

  }


  return(reached_states_bin)



}
