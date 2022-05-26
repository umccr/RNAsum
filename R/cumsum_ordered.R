##### Calculating cumulative sum for while keeping the original data order
cumsum_ordered <- function(x) {

  ##### Perform range standardization between 0 and 1, otherwise the negative values are summed up
  standarised <- standardization(x)

  ##### Sort and cumsum values
  sorted_cumsum <- cumsum(sort(standarised))

  ##### Restore the original elements order
  ordered_cumsum <- sorted_cumsum[ names(standarised) ]

  ##### Perform range standardization between 0 and 1, otherwise the negative values are summed up
  standarised_cumsum <- standardization(ordered_cumsum)

  ##### Clean the space and return output
  rm(x, standarised, sorted_cumsum, ordered_cumsum)
  return( standarised_cumsum )
}
