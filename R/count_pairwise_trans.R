
count_pairwise_trans <- function(net, method = c("SDDS","BNp","PEW"), params,
                                 states, steps = 1, repeats = 1000,
                                 asynchronous = TRUE, update_prob = NULL)
{

  if(!is.vector(states) & !is.matrix(states))
    stop("The value of the argument \"initial_states\" must be either a vector or a matrix.")

  if(!all(states == 0 | states == 1))
    stop("All the elements in \"states\" must be either zero or one.")


  if(is.vector(states)) {
    if(length(states)!=length(net$genes)) {
      stop("The number of variables (elements) in \"states\" must be equal to the number of network nodes.")
    }
    states_dec <- bin2dec(states, length(net$genes))
    num_states <- 1
  } else if(nrow(states)==1) {
    if(length(states)!=length(net$genes)) {
      stop("The number of variables (elements) in \"states\" must be equal to the number of network nodes.")
    }
    states_dec <- bin2dec(states[1,], length(net$genes))
    num_states <- 1
  }
  else {
    if(ncol(states)!=length(net$genes)) {
      stop("The number of variables (columns) in \"states\" must be equal to the number of network nodes.")
    }
    num_states <- nrow(states)

    states_dec <- as.vector(apply(states, 1, bin2dec, len=length(net$genes)))

  }



  if(!is.positive.integer(steps))
    stop("The value of the argument \"steps\" must be an integer.")


  if(!is.positive.integer(repeats))
    stop("The value of the argument \"repeats\" must be an integer.")


  if (!is.logical_value(asynchronous))
    stop("The value of the argument \"asynchronous\" must be logical (TRUE or FALSE).")

  if(!is.null(update_prob)) {
    if (asynchronous) {
      if (!is.all_non_negative_float(update_prob)) {
        if(is.vector(update_prob)) {
          if(length(update_prob)!=length(net$genes)) {
            stop("The length of \"update_prob\" must be a equal to the number of network nodes.")
          } else if (sum(update_prob) != 1) {
            stop("The sum of the \"update_prob\" values must be one.")
          }
        } else {
          stop("The value of the argument \"update_prob\" must be a vector.")
        }
      } else {
        stop("All \"update prob\" values must be non-negative and non-NA.")
      }
    } else {
      warning("Since \"asynchronous = FALSE\", ignoring \"update_prob\".")
    }
  }


  # the C code requires all interactions to be coded into one vector:
  # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
  inputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$input)))
  input_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$input)))))

  # Do the same for the transition functions.
  outputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$func)))
  output_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$func)))))

  # if (any(states != 0 & states != 1)) {
  #   cat("Error!")
  #   return(NA)
  # }
  #
  # if(is.vector(states)) {
  #   states_dec <- bin2dec(states, length(net$genes))
  #   num_states <- 1
  # }
  # else {
  #   states_dec <- as.vector(apply(states, 1, bin2dec, len=length(net$genes)))
  #   num_states <- nrow(states)
  #   #print(t(apply(apply(initial_states, 1, bin2dec, len=length(net$genes)),2,dec2bin,len=length(net$genes))))
  #
  # }


  switch(match.arg(method), SDDS={

    if (!is.list(params) || is.null(names(params)))
      stop("The value of the argument \"params\" must be a named list.")

    if (!all(c("p00", "p01", "p10", "p11") %in% names(params)))
      stop("The value of the argument \"params\" must be a named list consisting of \"p00\", \"p01\", \"p10\", and \"p11\".")


    if(length(params$p00) != length(net$genes) |
       length(params$p01) != length(net$genes) |
       length(params$p10) != length(net$genes) |
       length(params$p11) != length(net$genes))
      stop("The lengths of \"p00\", \"p01\", \"p10\", and \"p11\" must be equal to the number of network nodes.")


    if(!is.nonNA.numeric(params$p00) | !is.nonNA.numeric(params$p01) | !is.nonNA.numeric(params$p10) | !is.nonNA.numeric(params$p11))
      stop("The vectors\"p00\", \"p01\", \"p10\", and \"p11\" must be numeric without NA values.")



  if(asynchronous) {

    pairwise_transitions <- .Call("get_pairwise_transitions_SDDS_async_R", inputs, input_positions,
                             outputs, output_positions,
                             as.integer(net$fixed),
                             params$p00, params$p01, params$p10, params$p11,
                             update_prob, states_dec, num_states,
                             as.integer(steps), as.integer(repeats),
                             PACKAGE = "PARBONET")


    # SEXP inputs, SEXP input_positions,
    # SEXP outputs, SEXP output_positions,
    # SEXP fixed_nodes, SEXP p00, SEXP p01,
    # SEXP p10, SEXP p11, SEXP update_prob,
    # SEXP states, SEXP num_states,
    # SEXP steps, SEXP repeats


  } else {

    pairwise_transitions <- .Call("get_pairwise_transitions_SDDS_sync_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed),
                               params$p00, params$p01, params$p10, params$p11,
                               states_dec, num_states,
                               as.integer(steps), as.integer(repeats),
                               PACKAGE = "PARBONET")




  }

  },
  BNp={

    if(!is.nonNA.numeric(params))
      stop("The value of the argument \"params\" must be numeric vector without NA values.")


    if(length(params) != length(net$genes))
      stop("The length of \"params\" must be equal to the number of network nodes.")


    if(asynchronous) {

      pairwise_transitions <- .Call("get_pairwise_transitions_BNp_async_R", inputs, input_positions,
                                 outputs, output_positions,
                                 as.integer(net$fixed),
                                 params,
                                 update_prob, states_dec, num_states,
                                 as.integer(steps), as.integer(repeats),
                                 PACKAGE = "PARBONET")


      # SEXP inputs, SEXP input_positions,
      # SEXP outputs, SEXP output_positions,
      # SEXP fixed_nodes, SEXP p00, SEXP p01,
      # SEXP p10, SEXP p11, SEXP update_prob,
      # SEXP states, SEXP num_states,
      # SEXP steps, SEXP repeats


    } else {

      pairwise_transitions <- .Call("get_pairwise_transitions_BNp_sync_R", inputs, input_positions,
                                 outputs, output_positions,
                                 as.integer(net$fixed),
                                 params,
                                 states_dec, num_states,
                                 as.integer(steps), as.integer(repeats),
                                 PACKAGE = "PARBONET")


    }


  },
  PEW={

    if (!is.list(params) || is.null(names(params)))
      stop("The value of the argument \"params\" must be a named list.")

    if (!all(c("p_on", "p_off") %in% names(params)))
      stop("The value of the argument \"params\" must be a named list consisting of \"p_on\" and \"p_off\".")

    if(length(params$p_on) != nrow(extract_edges(net)) |
       length(params$p_off) != nrow(extract_edges(net)))
      stop("The lengths of \"p_on\" and \"p_off\" must be equal to the number of network edges.")

    if(!is.nonNA.numeric(params$p_on) | !is.nonNA.numeric(params$p_off))
      stop("The vectors \"p_on\" and \"p_off\" must be numeric without NA values.")


    if(asynchronous) {

      pairwise_transitions <- .Call("get_pairwise_transitions_PEW_async_R", inputs, input_positions,
                                 outputs, output_positions,
                                 as.integer(net$fixed),
                                 params$p_on, params$p_off, update_prob, states_dec, num_states,
                                 as.integer(steps), as.integer(repeats),
                                 PACKAGE = "PARBONET")


    } else {

      pairwise_transitions <- .Call("get_pairwise_transitions_PEW_sync_R", inputs, input_positions,
                                 outputs, output_positions,
                                 as.integer(net$fixed),
                                 params$p_on, params$p_off, states_dec, num_states,
                                 as.integer(steps), as.integer(repeats),
                                 PACKAGE = "PARBONET")




    }
  },
  stop("The value of the argument \"method\" must be one of \"SDDS\",\"BNp\",\"PEW\"")
  )

  # print(num_states)
  #
  # print(states_dec)

  pairwise_transitions <- matrix(pairwise_transitions, nrow = num_states, byrow = TRUE)
  #colnames(transition_matrix) <- net$genes

  return(pairwise_transitions)

}
