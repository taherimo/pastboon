extract_edges <- function(network, node_names= TRUE) {

  inputs <- as.integer(sapply(network$interactions, FUN = function(x) {x$input}))
  #edges <- list(source= unname(unlist(inputs)), destination= rep(1:length(network$genes),sapply(inputs, length)))
  if (node_names)
    edges <- data.frame(source= unname(unlist(inputs)), destination= rep(1:length(network$genes),sapply(inputs, length)))
  else
    edges <- data.frame(source= network$genes[unname(unlist(inputs))], destination= network$genes[rep(1:length(network$genes),sapply(inputs, length))])


  rownames(edges) <- 1:nrow(edges)
  colnames(edges) <- c("source", "destination")

  # if (node_names) {
  #   edges$source <- network$genes[edges$source]
  #   edges$destination <- network$genes[edges$destination]
  # }

  return(edges)

}
