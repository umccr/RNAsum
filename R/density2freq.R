##### Convert density to counts
density2freq <- function(density) {
  freq = length(density)/base::sum(density) * density
  return(freq)
}
