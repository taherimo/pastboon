
simulate <- function(net, p00, p01, p10, p11, initial_prob,
                     steps, repeats, return_last_step=F,
                     update_scheme="asynchronous") {

  if (!is.loaded("simulate_R")) dyn.load("package.so")

  # the C code requires all interactions to be coded into one vector:
  # Assemble all input gene lists in one list <inputGenes>, and remember the split positions in <inputGenePositions>.
  inputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$input)))
  input_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$input)))))

  # Do the same for the transition functions.
  outputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$func)))
  output_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$func)))))

  if(update_scheme=="asynchronous") {
    traj <- .Call("simulate_async_R", inputs, input_positions,
                  outputs, output_positions,
                  as.integer(net$fixed),
                  up_prob, down_prob,
                  initial_prob,
                  as.integer(steps),
                  as.integer(repeats),
                  as.double(noise),
                  as.integer(return_last_step))
  } else {
    traj <- .Call("simulate_sync_R", inputs, input_positions,
                  outputs, output_positions,
                  as.integer(net$fixed),
                  up_prob, down_prob,
                  initial_prob,
                  as.integer(steps),
                  as.integer(repeats),
                  as.double(noise),
                  as.integer(return_last_step))
  }




  if(return_last_step) {
    names(traj) <- net$genes
  }
  else {
    traj <- matrix(traj, nrow = steps+1, byrow = FALSE)
    colnames(traj) <- net$genes
  }

  return(traj)


}
