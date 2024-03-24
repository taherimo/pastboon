get_cumulative_transition_matrix <- function(net, method=c("SDDS","BNp","PEW"), params, states,
                                  time_step=1, repeats=1000,
                                  asynchronous=T, update_prob=NULL)
{

  # the C code requires all interactions to be coded into one vector:
  # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
  inputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$input)))
  input_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$input)))))

  # Do the same for the transition functions.
  outputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$func)))
  output_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$func)))))

  if (any(states != 0 & states != 1)) {
    cat("Error!")
    return(NA)
  }

  if(is.vector(states)) {
    states_dec <- bin2dec(states, length(net$genes))
    num_states <- 1
  }
  else {
    states_dec <- as.vector(apply(states, 1, bin2dec, len=length(net$genes)))
    num_states <- nrow(states)
    #print(t(apply(apply(initial_states, 1, bin2dec, len=length(net$genes)),2,dec2bin,len=length(net$genes))))

  }


  switch(match.arg(method), SDDS={

    p00 <- params$p00
    p01 <- params$p01
    p10 <- params$p10
    p11 <- params$p11


  if(asynchronous) {

    transition_matrix <- .Call("get_cumulative_transition_matrix_SDDS_async_R", inputs, input_positions,
                             outputs, output_positions,
                             as.integer(net$fixed),
                             p00, p01, p10, p11,
                             update_prob, states_dec, num_states,
                             as.integer(time_step), as.integer(repeats),
                             PACKAGE = "PARBONET")


    # SEXP inputs, SEXP input_positions,
    # SEXP outputs, SEXP output_positions,
    # SEXP fixed_nodes, SEXP p00, SEXP p01,
    # SEXP p10, SEXP p11, SEXP update_prob,
    # SEXP states, SEXP num_states,
    # SEXP steps, SEXP repeats


  } else {

    transition_matrix <- .Call("get_cumulative_transition_matrix_SDDS_sync_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed),
                               p00, p01, p10, p11,
                               states_dec, num_states,
                               as.integer(time_step), as.integer(repeats),
                               PACKAGE = "PARBONET")




  }

  },
  BNp={

    if(asynchronous) {

      transition_matrix <- .Call("get_cumulative_transition_matrix_BNp_async_R", inputs, input_positions,
                                 outputs, output_positions,
                                 as.integer(net$fixed),
                                 params,
                                 update_prob, states_dec, num_states,
                                 as.integer(time_step), as.integer(repeats),
                                 PACKAGE = "PARBONET")


      # SEXP inputs, SEXP input_positions,
      # SEXP outputs, SEXP output_positions,
      # SEXP fixed_nodes, SEXP p00, SEXP p01,
      # SEXP p10, SEXP p11, SEXP update_prob,
      # SEXP states, SEXP num_states,
      # SEXP steps, SEXP repeats


    } else {

      transition_matrix <- .Call("get_cumulative_transition_matrix_BNp_sync_R", inputs, input_positions,
                                 outputs, output_positions,
                                 as.integer(net$fixed),
                                 params,
                                 states_dec, num_states,
                                 as.integer(time_step), as.integer(repeats),
                                 PACKAGE = "PARBONET")


    }


  },
  PEW={

    p_on <- params$p_on
    p_off <- params$p_off


    if(asynchronous) {

      transition_matrix <- .Call("get_cumulative_transition_matrix_PEW_async_R", inputs, input_positions,
                                 outputs, output_positions,
                                 as.integer(net$fixed),
                                 p_on, p_off, update_prob, states_dec, num_states,
                                 as.integer(time_step), as.integer(repeats),
                                 PACKAGE = "PARBONET")


      # SEXP inputs, SEXP input_positions,
      # SEXP outputs, SEXP output_positions,
      # SEXP fixed_nodes, SEXP p00, SEXP p01,
      # SEXP p10, SEXP p11, SEXP update_prob,
      # SEXP states, SEXP num_states,
      # SEXP steps, SEXP repeats


    } else {

      transition_matrix <- .Call("get_cumulative_transition_matrix_PEW_sync_R", inputs, input_positions,
                                 outputs, output_positions,
                                 as.integer(net$fixed),
                                 p_on, p_off, states_dec, num_states,
                                 as.integer(time_step), as.integer(repeats),
                                 PACKAGE = "PARBONET")




    }
  },
  stop("'method' must be one of \"SDDS\",\"BNp\",\"PEW\"")
  )

  # print(num_states)
  #
  # print(states_dec)

  transition_matrix <- matrix(transition_matrix, nrow = num_states, byrow = TRUE)
  #colnames(transition_matrix) <- net$genes

}
