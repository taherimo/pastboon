get_pseudopotential_landscape <- function(bn, p00, p01, p10, p11, x_states) {


  x_states_bin <- t(apply(x_states,1,function(state) dec2bin(state,length(bn$genes))))



  return(x_states_bin)


}

get_attractor_states <- function(bn) {

  attractors <- getAttractors(bn, type = "synchronous", method="sat.exhaustive")


  attractor_states <- sapply(attractors$attractors, FUN = function(x) {

    t(apply(x$involvedStates,2,function(state)
      dec2bin(state,numGenes)))

  })

  return(attractor_states)

}


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

