
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


is.positive.integer <- function(x) {

  if (x != as.integer(x) || !is.numeric(x)) {
    return(FALSE)
  }

  if (x <= 0) {
    return(FALSE)
  }

  return(TRUE)

}


is.nonNA.numeric <- function(vector) {
  if (!is.numeric(vector)) {
    return(FALSE)
  }
  return(all(!is.na(vector)))
}

is.all_non_negative_float <- function(x) {
  is.numeric(x) && all(!is.na(x)) && all(x >= 0)
}

is.logical_value <- function(x) {
  is.logical(x) && length(x) == 1
}



is.scalar <- function(x) {
  is.atomic(x) && length(x) == 1
}



trim <- function(string)
{
  string <- gsub("^[ \t]+", "", string)
  string <- gsub("[ \t]+$", "", string)
  return(string)
}

