# This function is taken from the BoolNet package.
# Original file: BoolNet/R/helpers.R
# Original function: bin2dec
# Encode a vector of binary values <bin> with <len> bits
# to a decimal number
bin2dec <- function(bin, len) {
  if (len %% 32 == 0) {
    numElts <- as.integer(len / 32)
  } else {
    numElts <- as.integer(len / 32) + 1
  }

  dec <- rep(0, numElts)

  dec <- .C("bin2decC", as.integer(dec), as.integer(bin), as.integer(len))[[1]]
}

# This function is taken from the BoolNet package.
# Original file: BoolNet/R/helpers.R
# Original function: dec2bin
# Decode the <len> low-order bits of <dec> to a vector of binary values,
# where the first entry is the least significant bit
dec2bin <- function(dec, len) {
  bin <- rep(0, len)

  bin <- .C("dec2binC", as.integer(bin), as.integer(dec), as.integer(len), NAOK = TRUE)[[1]]
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

is.nonnegative.integer <- function(x) {
  if (x != as.integer(x) || !is.numeric(x)) {
    return(FALSE)
  }

  if (x < 0) {
    return(FALSE)
  }

  return(TRUE)
}

is.non_negative_real <- function(x) {
  is.numeric(x) && !is.na(x) && x >= 0
}

is.all_in_range_0_1 <- function(x) {
  is.numeric(x) && all(!is.na(x)) && all(x >= 0) && all(x <= 1)
}


is.logical_value <- function(x) {
  is.logical(x) && length(x) == 1
}



is.scalar <- function(x) {
  is.atomic(x) && length(x) == 1
}



trim <- function(string) {
  string <- gsub("^[ \t]+", "", string)
  string <- gsub("[ \t]+$", "", string)
  return(string)
}

is.BooleanNetwork <- function(net) {

  # Check if 'net' is a list
  if (!is.list(net)) {
    stop("The value of the argument \"net\" must be a named list consisting of \"interactions\", \"genes\", and \"fixed\".")
  } else if (!all(c("interactions", "genes", "fixed") %in% names(net))) {
    stop("The value of the argument \"net\" must be a named list consisting of \"interactions\", \"genes\", and \"fixed\".")
  }

  # Check if 'interactions' is a list of lists
  if (!is.list(net$interactions) || !all(sapply(net$interactions, is.list))) {
    stop("\"interactions\" must be a list of lists.")
  }

  if (!is.character(net$genes)) {
    stop("\"genes\" must be a character vector.")
  }

  if (!all(net$fixed == as.integer(net$fixed)) || !is.numeric(net$fixed)) {
    if (!all(net$fixed %in% c(0, 1, -1))) {
      stop("\"fixed\" must be numeric vector consisting of the values: 0, 1, or -1.")
    }
  }


  # Check if the lengths of 'genes', 'fixed', and 'interactions' are equal
  len_genes <- length(net$genes)
  len_fixed <- length(net$fixed)
  len_interactions <- length(net$interactions)

  if (!(len_genes == len_fixed && len_fixed == len_interactions)) {
    stop("The lengths of \"genes\", \"fixed\", and \"interactions\" must be equal.")
  }


  # Check if 'fixed' contains only 0, 1, or -1
  if (!all(net$fixed %in% c(0, 1, -1))) {
    stop("\"fixed\" must only contain the values 0, 1, or -1.")
  }


  return(TRUE)
}
