calc_convergence_time <- function(node_act, threshold, window_size = 1) {


  if (!is.matrix(node_act)) {
    stop("The value of the argument \"node_act\" must be a matrix.")
  }

  if (!is.scalar(threshold)) {
    stop("The value of the argument \"threshold\" must be a scalar.")
  }

  if (!is.non_negative_real(threshold)) {
    stop("The value of the argument \"threshold\" must be a non-negative real.")
  }

  if (!is.scalar(window_size)) {
    stop("The value of the argument \"window_size\" must be a scalar.")
  }

  if (!is.positive.integer(window_size)) {
    stop("The value of the argument \"window_size\" must be a positive integer.")
  }

  if (nrow(node_act) < window_size) {
    stop("The rows in \"node_act\" (time-steps) must be equal to or greater than \"window_size\".")
  }

  differences <- diff(node_act, 1, along = 1)

  rows_below_threshold <- apply(differences, 1, function(row) all(row < threshold))

  convergence_time <- find_consecutive_true(rows_below_threshold, window_size)

  if (is.na(convergence_time)) {
    # No convergence was detected!
    return(NA)
  }

  return(convergence_time + 1)
}

find_consecutive_true <- function(vec, n) {
  for (i in 1:(length(vec) - n + 1)) {
    if (all(vec[i:(i + n - 1)])) {
      return(i)
    }
  }
  return(NA)
}
