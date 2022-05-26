##### Check for nearest value in a vector
nearest_position <- function(vector, x) {

  y <- which.min(abs(vector - x))

  ##### Clean the space and return output
  rm(vector, x)
  return( y )
}
