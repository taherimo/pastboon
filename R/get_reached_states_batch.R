get_reached_states_batch <- function(net, method=c("SDDS","BNp","PEW"), params,
                                      steps, repeats=1, num_initial_states, initial_states=NULL,
                                      update_prob=NULL,asynchronous=T) {

  if(!is.positive.integer(reoeats)) {
   stop("The value of the argument repeats is not integer.")
  }

  if(!is.positive.integer(steps)) {
    stop("The value of the argument steps is not integer.")
  }

  initial_states_dec <- NULL

  # assert ncols(initial_states) == length(net$genes)
  if(!is.null(initial_states)) {
    if(!all(initial_states == 0 | initial_states == 1)) {
      stop("Non-binary value(s) in initial_states.")
    }
    if(is.vector(initial_states)) {
      if(length(initial_states)!=length(net$genes)) {
        stop("The number of variables in initial_states doesn't match the number of network nodes.")
      }
      initial_states_dec <- bin2dec(initial_states, length(net$genes))
    } else {
      if(ncol(initial_states)!=length(net$genes)) {
        stop("The number of variables in initial_states doesn't match the number of network nodes.")
      }
      if(nrow(initial_states)>1) {
        if(repeats>1) {
          stop("in the case of repeats>1, a single initial state or no initial state can be given")
        }
        initial_states_dec <- as.vector(apply(initial_states, 1, bin2dec, len=length(net$genes)))
      }

    }
  }

  if(!is.null(update_prob)) {
    if(is.vector(update_prob)) {
      if(length(update_prob)!=length(net$genes)) {
        stop("The length of update_prob should be a equal to the number of network nodes.")
      }
    } else {
      stop("The argument update_prob should be a vector.")
    }
  }

  # if(is.null(initial_states)) {
  #   initial_states_bin <- sample(0:1, num_random_initial_states * length(bn$genes), rep = T)
  #   if(num_random_initial_states>1) {
  #     initial_states_bin <- matrix(initial_states_bin, nrow = num_random_initial_states, byrow = TRUE)
  #     colnames(initial_states_bin) <- net$genes
  #   } else {
  #     names(initial_states_bin) <- net$genes
  #   }
  # }

  # initial_states_dec <- NULL
  #
  # if(!is.null(initial_states)) {
  #
  #   # assert ncols(initial_states) == length(net$genes)
  #   # assert nrow(initial_states) == num_initial_states --> warning
  #   # at least one of the args initial_states or num_initial_states should be passed
  #
  #   if (any(initial_states != 0 & initial_states != 1)) {
  #     cat("Error!") # there are values other than zero and one
  #     return(NA)
  #   }
  #
  #   if(is.vector(initial_states)) {
  #     if(num_initial_states!=1) {
  #       cat("Warning!")
  #       num_initial_states <- 1
  #     }
  #     initial_states_dec <- bin2dec(initial_states, length(net$genes))
  #   }
  #   else {
  #     initial_states_dec <- as.vector(apply(initial_states, 1, bin2dec, len=length(net$genes)))
  #     #print(t(apply(apply(initial_states, 1, bin2dec, len=length(net$genes)),2,dec2bin,len=length(net$genes))))
  #
  #     if(num_initial_states!=nrow(initial_states)) {
  #       cat("Warning!")
  #       num_initial_states <- nrow(initial_states)
  #     }
  #   }
  #
  # } else if(num_initial_states <=0) {
  #     cat("Error!") # guarantee that num_initial_states is positive integer
  #     return(NA)
  # }



  # the C code requires all interactions to be coded into one vector:
  # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
  inputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$input)))
  input_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$input)))))

  # Do the same for the transition functions.
  outputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$func)))
  output_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$func)))))


  if (!is.list(params) || is.null(names(params))) {
    stop("The params argument must be a named list.")
  }

  switch(match.arg(method), SDDS={


    if (!all(c("p00", "p01", "p10", "p11") %in% names(params))) {
      stop("Input list must contain vectors named 'p00', 'p01', 'p10', and 'p11'.")
    }

    if(length(params$p00) != length(net$genes) |
       length(params$p01) != length(net$genes) |
       length(params$p10) != length(net$genes) |
       length(params$p11) != length(net$genes)) {
      stop("Length of p00, p01, p10, and p11 must be equal to the number of network nodes.")
    }

    if(!is.nonNA.numeric(params$p00) | !is.nonNA.numeric(params$p01) | !is.nonNA.numeric(params$p10) | !is.nonNA.numeric(params$p11)) {
      stop("The vectors p00, p01, p10, and p11 must be numeric without NA values.")
    }


    if(num_initial_states==1) {

      if(asynchronous) {


        reached_states <- .Call("get_reached_states_SDDS_async_single_R", inputs, input_positions,
                                outputs, output_positions,
                                as.integer(net$fixed),
                                params$p00, params$p01, params$p10, params$p11,
                                update_prob, as.integer(initial_state_dec),
                                as.integer(repeats), as.integer(steps),
                                PACKAGE = "PARBONET")
      } else {


        reached_states <- .Call("get_reached_states_SDDS_sync_single_R", inputs, input_positions,
                                outputs, output_positions,
                                as.integer(net$fixed),
                                params$p00, params$p01, params$p10, params$p11,
                                as.integer(initial_state_dec),
                                as.integer(repeats), as.integer(steps),
                                PACKAGE = "PARBONET")

      }

    } else {

      if(asynchronous) {

        #reached_states <- .Call("simulate_async_return_states_R", inputs, input_positions,
        #      outputs, output_positions, as.integer(net$fixed), p00, p01, p10, p11,
        #      as.integer(initial_states), update_prob, as.integer(steps))

        reached_states <- .Call("get_reached_states_SDDS_async_batch_R", inputs, input_positions,
                                outputs, output_positions, as.integer(net$fixed), params$p00, params$p01, params$p10, params$p11,
                                as.integer(initial_states_dec), as.integer(num_initial_states), update_prob, as.integer(steps))

      } else {
        reached_states <- .Call("get_reached_states_SDDS_sync_batch_R", inputs, input_positions,
                                outputs, output_positions, as.integer(net$fixed), params$p00, params$p01, params$p10, params$p11,
                                as.integer(initial_states_dec), as.integer(num_initial_states), as.integer(steps))

      }


    }

  },
  BNp={

    if(length(params) != length(net$genes)) {
      stop("Length of params must be equal to the number of network nodes.")
    }

    if(!is.nonNA.numeric(params)) {
      stop("The vector params must be numeric without NA values.")
    }


    if(num_initial_states==1) {

      if(asynchronous) {


        reached_states <- .Call("get_reached_states_BNp_async_single_R", inputs, input_positions,
                                outputs, output_positions,
                                as.integer(net$fixed),
                                params,
                                update_prob, as.integer(initial_state_dec),
                                as.integer(repeats), as.integer(steps),
                                PACKAGE = "PARBONET")
      } else {


        reached_states <- .Call("get_reached_states_BNp_sync_single_R", inputs, input_positions,
                                outputs, output_positions,
                                as.integer(net$fixed),
                                params,
                                as.integer(initial_state_dec),
                                as.integer(repeats), as.integer(steps),
                                PACKAGE = "PARBONET")

      }

    } else {

      if(asynchronous) {

        #reached_states <- .Call("simulate_async_return_states_R", inputs, input_positions,
        #      outputs, output_positions, as.integer(net$fixed), p00, p01, p10, p11,
        #      as.integer(initial_states), update_prob, as.integer(steps))

        reached_states <- .Call("get_reached_states_BNp_async_batch_R", inputs, input_positions,
                                outputs, output_positions, as.integer(net$fixed), params,
                                as.integer(initial_states_dec), as.integer(num_initial_states), update_prob, as.integer(steps))

      } else {
        reached_states <- .Call("get_reached_states_BNp_sync_batch_R", inputs, input_positions,
                                outputs, output_positions, as.integer(net$fixed), params,
                                as.integer(initial_states_dec), as.integer(num_initial_states), as.integer(steps))

      }


    }

  },
  PEW={


    if (!all(c("p_on", "p_off") %in% names(params))) {
      stop("Input list must contain vectors named 'p_on' and 'p_off'.")
    }

    if(length(params$p_on) != length(net$genes) |
       length(params$p_off) != length(net$genes)) {
      stop("Length of p_on and p_off must be equal to the number of network nodes.")
    }

    if(!is.nonNA.numeric(params$p_on) | !is.nonNA.numeric(params$p_off)) {
      stop("The vectors p_on and p_off must be numeric without NA values.")
    }


    if(num_initial_states==1) {

      if(asynchronous) {


        reached_states <- .Call("get_reached_states_PEW_async_single_R", inputs, input_positions,
                                outputs, output_positions,
                                as.integer(net$fixed),
                                params$p_on, params$p_off, update_prob,
                                as.integer(initial_state_dec),
                                as.integer(repeats), as.integer(steps),
                                PACKAGE = "PARBONET")
      } else {


        reached_states <- .Call("get_reached_states_PEW_sync_single_R", inputs, input_positions,
                                outputs, output_positions,
                                as.integer(net$fixed),
                                params$p_on, params$p_off, as.integer(initial_state_dec),
                                as.integer(repeats), as.integer(steps),
                                PACKAGE = "PARBONET")

      }

    } else {

      if(asynchronous) {

        #reached_states <- .Call("simulate_async_return_states_R", inputs, input_positions,
        #      outputs, output_positions, as.integer(net$fixed), p00, p01, p10, p11,
        #      as.integer(initial_states), update_prob, as.integer(steps))

        reached_states <- .Call("get_reached_states_PEW_async_batch_R", inputs, input_positions,
                                outputs, output_positions, as.integer(net$fixed), params$p_on, params$p_off,
                                as.integer(initial_states_dec), as.integer(num_initial_states), update_prob, as.integer(steps))

      } else {
        reached_states <- .Call("get_reached_states_PEW_sync_batch_R", inputs, input_positions,
                                outputs, output_positions, as.integer(net$fixed), params$p_on, params$p_off,
                                as.integer(initial_states_dec), as.integer(num_initial_states), as.integer(steps))

      }


    }
  },
  stop("'method' must be one of \"SDDS\",\"BNp\",\"PEW\"")

  )





  #print(reached_states)


  if(num_initial_states==1) {

    reached_states_bin <- dec2bin(reached_states,len=length(net$genes))

  } else {

    reached_states <- matrix(reached_states, nrow = num_initial_states, byrow = TRUE)

    reached_states_bin <- apply(reached_states, 1, dec2bin, len=length(net$genes))


    reached_states_bin <- t(reached_states_bin)[,1:length(net$genes)]
    colnames(reached_states_bin) <- net$genes

  }


  return(reached_states_bin)


}
