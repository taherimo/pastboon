
calc_node_activities <- function(net, method=c("SDDS","BNp","PEW"), params,
                     steps, repeats = 1000, initial_prob= NULL, last_step = FALSE,
                     asynchronous = TRUE, update_prob = NULL) {

  if(!is.positive.integer(steps)) {
    stop("The value of the argument \"steps\" is not integer.")
  }

  if(!is.positive.integer(repeats)) {
    stop("The value of the argument \"repeats\" is not integer.")
  }

  if(!is.null(initial_prob)) {
    if(is.vector(initial_prob)) {
      if(length(initial_prob)!=length(net$genes)) {
        stop("The length of initial_prob should be a equal to the number of network nodes.")
      }
    } else {
      stop("The argument initial_prob should be a vector.")
    }
  }

  if(!is.null(update_prob)) {
    if (!is.all_non_negative_float(update_prob)) {
      if(is.vector(update_prob)) {
        if(length(update_prob)!=length(net$genes)) {
          stop("The length of update_prob should be a equal to the number of network nodes.")
        } else if (sum(update_prob) != 1) {
          stop("The sum of the update_prob values should be 1.")
        }
      } else {
        stop("The argument update_prob should be a vector.")
      }
    } else {
      stop("All update prob values should be non-negative and non-NA.")
    }
  }


  if (!is.logical_value(last_step))
    stop("The value of the last_step argument should be logical (TRUE or FALSE).")

  if (!is.logical_value(asynchronous))
    stop("The value of the asybchronous argument should be logical (TRUE or FALSE).")



  # the C code requires all interactions to be coded into one vector:
  # Assemble all input gene lists in one list <inputs>, and remember the split positions in <input_positions>.
  inputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$input)))
  input_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$input)))))

  # Do the same for the transition functions.
  outputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$func)))
  output_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$func)))))


  switch(match.arg(method), SDDS={

    p00 <- params$p00
    p01 <- params$p01
    p10 <- params$p10
    p11 <- params$p11

    if(asynchronous) {

      node_activities <- .Call("get_node_activities_SDDS_async_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed),
                               p00, p01, p10, p11,
                               initial_prob, update_prob,
                               as.integer(steps), as.integer(repeats),
                               as.integer(last_step), PACKAGE = "PARBONET")


    } else {



      node_activities <- .Call("get_node_activities_SDDS_sync_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed),
                               p00, p01, p10, p11,
                               initial_prob,
                               as.integer(steps),
                               as.integer(repeats),
                               as.integer(last_step), PACKAGE = "PARBONET")

    }


  },
  BNp={

    if(asynchronous) {
      node_activities <- .Call("get_node_activities_BNp_async_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed), params,
                               initial_prob, update_prob,
                               as.integer(steps), as.integer(repeats),
                               as.integer(last_step), PACKAGE = "PARBONET")

    } else {

      node_activities <- .Call("get_node_activities_BNp_sync_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed), params,
                               initial_prob, as.integer(steps), as.integer(repeats),
                               as.integer(last_step), PACKAGE = "PARBONET")
    }

  },
  PEW={

    p_on <- params$p_on
    p_off <- params$p_off

    if(length(p_on) != length(inputs))
      stop("The length of \"p_on\" should be equal to the number of edges!")

    if(length(p_off) != length(inputs))
      stop("The length of \"p_off\" should be equal to the number of edges!")


    if(asynchronous) {

      node_activities <- .Call("get_node_activities_PEW_async_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed),
                               p_on, p_off, initial_prob, update_prob,
                               as.integer(steps), as.integer(repeats),
                               as.integer(last_step), PACKAGE = "PARBONET")


    } else {



      node_activities <- .Call("get_node_activities_PEW_sync_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed),
                               p_on, p_off, initial_prob,
                               as.integer(steps),
                               as.integer(repeats),
                               as.integer(last_step), PACKAGE = "PARBONET")

    }

  },
  stop("'method' must be one of \"SDDS\",\"BNp\",\"PEW\"")
  )



  if(last_step) {
    names(node_activities) <- net$genes
  }
  else {
    node_activities <- matrix(node_activities, nrow = steps+1, byrow = FALSE)
    colnames(node_activities) <- net$genes
    rownames(node_activities) <- 1:(steps+1)
  }

  return(node_activities)


}

