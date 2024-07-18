get_reached_states <- function(net, method = c("BNp", "SDDS", "PEW"), params,
                               steps, repeats = NULL, initial_states = NULL,
                               asynchronous = TRUE, update_prob = NULL) {
  if (!is.BooleanNetwork(net)) {
    stop("The value of the argument \"net\" must accord to the \"BooleanNetwork\" definition in \"BoolNet\".")
  }


  if (!is.scalar(steps)) {
    stop("The value of the argument \"steps\" must be a scalar.")
  }

  if (!is.nonnegative.integer(steps)) {
    stop("The value of the argument \"steps\" must be a non-negative integer.")
  }

  if (!is.null(repeats)) {
    if (!is.scalar(repeats)) {
      stop("The value of the argument \"repeats\" must be a scalar.")
    }

    if (!is.positive.integer(repeats)) {
      stop("The value of the argument \"repeats\" must be a positive integer.")
    }

  }


  if (!is.vector(initial_states) & !is.matrix(initial_states) & !is.null(initial_states)) {
    stop("The value of the argument \"initial_states\" must be either a vector, a matrix or NULL.")
  }


  initial_states_dec <- NULL

  # assert ncols(initial_states) == length(net$genes)
  if (!is.null(initial_states)) {
    if (is.null(repeats)) {
      repeats <- 1
    }
    if (!all(initial_states == 0 | initial_states == 1)) {
      stop("All the elements in \"states\" must be either zero or one.")
    }

    if (is.vector(initial_states)) {
      if (length(initial_states) != length(net$genes)) {
        stop("The number of variables in initial_states doesn't match the number of network nodes.")
      }
      initial_states_dec <- bin2dec(initial_states, length(net$genes))
      num_initial_states <- 1
    } else if (nrow(initial_states) == 1) {
      if (length(initial_states) != length(net$genes)) {
        stop("The number of variables in initial_states doesn't match the number of network nodes.")
      }
      initial_states_dec <- bin2dec(initial_states[1, ], length(net$genes))
      num_initial_states <- 1
    } else {
      if (ncol(initial_states) != length(net$genes)) {
        stop("The number of variables in initial_states doesn't match the number of network nodes.")
      }
      num_initial_states <- nrow(initial_states)
      if (repeats > 1) {
        stop("In the case of \"repeats > 1\", a single initial state or no initial state can be given.")
      }

      # repeats is NULL or 1
      initial_states_dec <- as.vector(apply(initial_states, 1, bin2dec, len = length(net$genes)))
    }
  } else if (!is.null(repeats)) {
    num_initial_states <- repeats
  } else {
    stop("The values of the arguments \"repeats\" and \"initial_states\" cannot be NULL at the same time!")
  }

  if (!is.logical_value(asynchronous)) {
    stop("The value of the argument \"asynchronous\" must be logical (TRUE or FALSE).")
  }

  if (!is.null(update_prob)) {
    if (asynchronous) {
      if (is.vector(update_prob)) {

        if (length(update_prob) != length(net$genes)) {
          stop("The length of \"update_prob\" must be a equal to the number of network nodes.")
        }

        if (!is.all_in_range_0_1(update_prob)) {
          stop("All values in \"initial_prob\" must be in the range [0,1].")
        }

        if (sum(update_prob) != 1) {
          stop("The sum of the \"update_prob\" values must be one.")
        }

      } else {
        stop("The value of the argument \"update_prob\" must be a vector.")
      }

    } else {
      warning("Since \"asynchronous = FALSE\", ignoring \"update_prob\".")
    }
  }

  # This code is derived from the BoolNet package.
  # Original file: BoolNet/R/getAttractors.R
  # Original function: getAttractors
  # Modifications: Changed the naming of the variables.

  # The start of the derived code from BoolNet

  inputs <- as.integer(unlist(lapply(net$interactions, function(interaction) interaction$input)))
  input_positions <- as.integer(cumsum(c(0, sapply(net$interactions, function(interaction) length(interaction$input)))))

  outputs <- as.integer(unlist(lapply(net$interactions, function(interaction) interaction$func)))
  output_positions <- as.integer(cumsum(c(0, sapply(net$interactions, function(interaction) length(interaction$func)))))

  # The end of the derived code from BoolNet

  switch(match.arg(method),
     BNp = {

       if (length(params) != length(net$genes)) {
         stop("The length of \"params\" must be equal to the number of network nodes.")
       }

       if (!is.all_in_range_0_1(params)) {
         stop("The value of the argument \"params\" must be a vector consisting of values in the range [0,1].")
       }


       if (num_initial_states == 1) {
         if (asynchronous) {
           reached_states <- .Call("get_reached_states_BNp_async_single_R", inputs, input_positions,
                                   outputs, output_positions,
                                   as.integer(net$fixed),
                                   params,
                                   update_prob, as.integer(initial_states_dec),
                                   as.integer(repeats), as.integer(steps),
                                   PACKAGE = "pastboon"
           )
         } else {
           reached_states <- .Call("get_reached_states_BNp_sync_single_R", inputs, input_positions,
                                   outputs, output_positions,
                                   as.integer(net$fixed),
                                   params,
                                   as.integer(initial_states_dec),
                                   as.integer(repeats), as.integer(steps),
                                   PACKAGE = "pastboon"
           )
         }
       } else {
         if (asynchronous) {
           reached_states <- .Call(
             "get_reached_states_BNp_async_batch_R", inputs, input_positions,
             outputs, output_positions, as.integer(net$fixed), params,
             as.integer(initial_states_dec), as.integer(num_initial_states), update_prob, as.integer(steps)
           )
         } else {
           reached_states <- .Call(
             "get_reached_states_BNp_sync_batch_R", inputs, input_positions,
             outputs, output_positions, as.integer(net$fixed), params,
             as.integer(initial_states_dec), as.integer(num_initial_states), as.integer(steps)
           )
         }
       }
     },
    SDDS = {
      if (!is.list(params) || is.null(names(params))) {
        stop("The value of the argument \"params\" must be a named list.")
      }


      if (!all(c("p00", "p01", "p10", "p11") %in% names(params))) {
        stop("The value of the argument \"params\" must be a named list consisting of \"p00\", \"p01\", \"p10\", and \"p11\".")
      }

      if (length(params$p00) != length(net$genes) |
        length(params$p01) != length(net$genes) |
        length(params$p10) != length(net$genes) |
        length(params$p11) != length(net$genes)) {
        stop("The lengths of \"p00\", \"p01\", \"p10\", and \"p11\" must be equal to the number of network nodes.")
      }

      if (!is.all_in_range_0_1(params$p00) | !is.all_in_range_0_1(params$p01) | !is.all_in_range_0_1(params$p10) | !is.all_in_range_0_1(params$p11)) {
        stop("The vectors\"p00\", \"p01\", \"p10\", and \"p11\" must consist of values in the range [0,1].")
      }

      if (num_initial_states == 1) {
        if (asynchronous) {
          reached_states <- .Call("get_reached_states_SDDS_async_single_R", inputs, input_positions,
            outputs, output_positions,
            as.integer(net$fixed),
            params$p00, params$p01, params$p10, params$p11,
            update_prob, as.integer(initial_states_dec),
            as.integer(repeats), as.integer(steps),
            PACKAGE = "pastboon"
          )
        } else {
          reached_states <- .Call("get_reached_states_SDDS_sync_single_R", inputs, input_positions,
            outputs, output_positions,
            as.integer(net$fixed),
            params$p00, params$p01, params$p10, params$p11,
            as.integer(initial_states_dec),
            as.integer(repeats), as.integer(steps),
            PACKAGE = "pastboon"
          )
        }
      } else {
        if (asynchronous) {
          reached_states <- .Call(
            "get_reached_states_SDDS_async_batch_R", inputs, input_positions,
            outputs, output_positions, as.integer(net$fixed), params$p00, params$p01, params$p10, params$p11,
            as.integer(initial_states_dec), as.integer(num_initial_states), update_prob, as.integer(steps)
          )
        } else {
          reached_states <- .Call(
            "get_reached_states_SDDS_sync_batch_R", inputs, input_positions,
            outputs, output_positions, as.integer(net$fixed), params$p00, params$p01, params$p10, params$p11,
            as.integer(initial_states_dec), as.integer(num_initial_states), as.integer(steps)
          )
        }
      }
    },
    PEW = {
      if (!is.list(params) || is.null(names(params))) {
        stop("The value of the argument \"params\" must be a named list.")
      }

      if (!all(c("p_on", "p_off") %in% names(params))) {
        stop("The value of the argument \"params\" must be a named list consisting of \"p_on\" and \"p_off\".")
      }

      if (length(params$p_on) != nrow(extract_edges(net)) |
        length(params$p_off) != nrow(extract_edges(net))) {
        stop("The lengths of \"p_on\" and \"p_off\" must be equal to the number of network edges.")
      }

      if (!is.all_in_range_0_1(params$p_on) | !is.all_in_range_0_1(params$p_off)) {
        stop("The vectors \"p_on\" and \"p_off\" must consist of values in the range [0,1].")
      }

      if (num_initial_states == 1) {
        if (asynchronous) {
          reached_states <- .Call("get_reached_states_PEW_async_single_R", inputs, input_positions,
            outputs, output_positions,
            as.integer(net$fixed),
            params$p_on, params$p_off, update_prob,
            as.integer(initial_states_dec),
            as.integer(repeats), as.integer(steps),
            PACKAGE = "pastboon"
          )
        } else {
          reached_states <- .Call("get_reached_states_PEW_sync_single_R", inputs, input_positions,
            outputs, output_positions,
            as.integer(net$fixed),
            params$p_on, params$p_off, as.integer(initial_states_dec),
            as.integer(repeats), as.integer(steps),
            PACKAGE = "pastboon"
          )
        }
      } else {
        if (asynchronous) {
          reached_states <- .Call(
            "get_reached_states_PEW_async_batch_R", inputs, input_positions,
            outputs, output_positions, as.integer(net$fixed), params$p_on, params$p_off,
            as.integer(initial_states_dec), as.integer(num_initial_states), update_prob, as.integer(steps)
          )
        } else {
          reached_states <- .Call(
            "get_reached_states_PEW_sync_batch_R", inputs, input_positions,
            outputs, output_positions, as.integer(net$fixed), params$p_on, params$p_off,
            as.integer(initial_states_dec), as.integer(num_initial_states), as.integer(steps)
          )
        }
      }
    },
    stop("The value of the argument \"method\" must be one of \"SDDS\",\"BNp\",\"PEW\"")
  )

  if (num_initial_states == 1) {
    reached_states_bin <- dec2bin(reached_states, len = length(net$genes))
    names(reached_states_bin) <- net$genes
  } else {
    reached_states <- matrix(reached_states, nrow = num_initial_states, byrow = TRUE)

    reached_states_bin <- apply(reached_states, 1, dec2bin, len = length(net$genes))


    # reached_states_bin <- t(reached_states_bin)[,1:length(net$genes)]
    reached_states_bin <- t(reached_states_bin)
    colnames(reached_states_bin) <- net$genes
  }


  return(reached_states_bin)
}
