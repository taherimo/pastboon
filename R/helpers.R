
# Encode a vector of binary values <bin> with <len> bits
# to a decimal number
bin2dec <- function(bin,len)
{
  if (len %% 32 == 0)
    numElts <- as.integer(len / 32)
  else
    numElts <- as.integer(len / 32) + 1

  dec = rep(0,numElts)

  dec = .C("bin2decC",as.integer(dec),as.integer(bin),as.integer(len))[[1]]
}

# Decode the <len> low-order bits of <dec> to a vector of binary values,
# where the first entry is the least significant bit
dec2bin <- function(dec,len)
{
  bin = rep(0,len)

  bin = .C("dec2binC",as.integer(bin),as.integer(dec),as.integer(len),NAOK=TRUE)[[1]]
}

get_edges <- function(network) {

  inputs <- sapply(network$interactions, FUN = function(x) {x$input})
  edges <- cbind(unlist(inputs),rep(1:length(network$genes),sapply(inputs, length)))
  rownames(edges) <- 1:nrow(edges)
  colnames(edges) <- c("source", "detination")

  return(edges)

}


is.positive.integer <- function(x) {

  if (!is.integer(x) || !is.numeric(x)) {
    return(FALSE)
  }

  if (x <= 0) {
    return(FALSE)
  }

  return(TRUE)

}

