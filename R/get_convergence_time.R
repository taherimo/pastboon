
get_convergence_time <- function(node_activities, threshold, window_size=1) {
  # the first row is consideres as time 0

  # check dimensions of the matrix node_activities

  differences <- diff(node_activities, 1, along = 1)

  #first_row_below_threshold <- which(apply(matrix_data, 1, function(row) all(row < threshold)), arr.ind = TRUE)[1, 1]

  rows_below_threshold <- apply(differences, 1, function(row) all(row < threshold))

  convergence_time <- find_consecutive_true(rows_below_threshold, window_size)

  if(is.na(convergence_time)) {
    cat("No convergence time found! Try reducing threshold and/or window_size")
  }

  return(convergence_time + 1)

}

find_consecutive_true <- function(vec, n) {
  for (i in 1:(length(vec) - n + 1)) {
    if (all(vec[i:(i+n-1)])) {
      return(i)
    }
  }
  return(NA)  # Return NA if no consecutive TRUE values are found
}
