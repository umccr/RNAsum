##### Convert a vector of numbers into corresponding vector of their percentiles
perc_rank <- function(x) trunc(rank(x))*100/length(x)
