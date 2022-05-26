##### Perform range standardization between 0 and 1 (for the cumulative sums)
standardization <- function(x) c(x-min(x))/(max(x)-min(x))
