extract_edges <- function(net, node_names = TRUE) {
  if (!is.BooleanNetwork(net)) {
    stop("The value of the argument \"net\" must accord to the \"BooleanNetwork\" definition in \"BoolNet\".")
  }


  if (!is.logical_value(node_names)) {
    stop("The value of the argument \"node_names\" must be logical (TRUE or FALSE).")
  }


  inputs <- sapply(net$interactions, FUN = function(x) {
    x$input
  })

  # edges <- list(source= unname(unlist(inputs)), destination= rep(1:length(network$genes),sapply(inputs, length)))
  # if (node_names) {
  #   edges$source <- network$genes[edges$source]
  #   edges$destination <- network$genes[edges$destination]
  # }

  if (node_names) {
    edges <- data.frame(source = net$genes[unname(unlist(inputs))], destination = net$genes[rep(1:length(net$genes), sapply(inputs, length))])
  } else {
    edges <- data.frame(source = unname(unlist(inputs)), destination = rep(1:length(net$genes), sapply(inputs, length)))
  }

  rownames(edges) <- 1:nrow(edges)
  colnames(edges) <- c("source", "destination")


  return(edges)
}
