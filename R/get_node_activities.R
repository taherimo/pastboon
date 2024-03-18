
get_node_activities <- function(net, method=c("SDDS","BNp"), params,
                     steps, repeats=1000, return_last_step=F,
                     asynchronous=T, initial_prob=NULL, update_prob=NULL) {



  # the C code requires all interactions to be coded into one vector:
  # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
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
                               as.integer(return_last_step), PACKAGE = "PARBONET")


    } else {



      node_activities <- .Call("get_node_activities_SDDS_sync_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed),
                               p00, p01, p10, p11,
                               initial_prob,
                               as.integer(steps),
                               as.integer(repeats),
                               as.integer(return_last_step), PACKAGE = "PARBONET")

    }


  },
  BNp={

    if(asynchronous) {

    } else {

      node_activities <- .Call("get_node_activities_BNp_sync_R", inputs, input_positions,
                               outputs, output_positions,
                               as.integer(net$fixed),
                               p00, p01, p10, p11,
                               initial_prob, update_prob,
                               as.integer(steps), as.integer(repeats),
                               as.integer(return_last_step), PACKAGE = "PARBONET")
    }

  },
  stop("'method' must be one of \"SDDS\",\"BNp\"")
  )



  if(return_last_step) {
    names(node_activities) <- net$genes
  }
  else {
    node_activities <- matrix(node_activities, nrow = steps+1, byrow = FALSE)
    colnames(node_activities) <- net$genes
    rownames(node_activities) <- 1:(steps+1)
  }

  return(node_activities)


}

