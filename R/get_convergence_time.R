
get_convergence_time <- function(net, p00, p01, p10, p11,
                                 max_steps, threshold, repeats,
                                 asynchronous=T, initial_prob=NULL, update_prob=NULL) {


  inputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$input)))
  input_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$input)))))

  outputs <- as.integer(unlist(lapply(net$interactions,function(interaction)interaction$func)))
  output_positions <- as.integer(cumsum(c(0,sapply(net$interactions,function(interaction)length(interaction$func)))))

  if(asynchronous) {

    convergence_time <- .Call("get_convergence_time_async_R", inputs, input_positions,
                             outputs, output_positions,
                             as.integer(net$fixed),
                             p00, p01, p10, p11,
                             initial_prob, update_prob,
                             as.integer(repeats), threshold,
                             as.integer(max_steps), PACKAGE = "PARBONET")


  } else {

    convergence_time <- .Call("get_convergence_time_sync_R", inputs, input_positions,
                             outputs, output_positions,
                             as.integer(net$fixed),
                             p00, p01, p10, p11,
                             initial_prob,
                             as.integer(repeats),
                             threshold, as.integer(max_steps),
                             PACKAGE = "PARBONET")

  }


  return(convergence_time)

}
