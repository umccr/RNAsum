##### Calculating cumulative sum for while keeping the original data order
#' Calculate cumulative sum
#'
#' @param x Input
#'
#' @return Cumulative sum while keeping the order of input data
#' @export
#'

cumsum_ordered <- function(x) {

  ##### Perform range standardization between 0 and 1, otherwise the negative values are summed up
  standarised <- standardization(x)

  ##### Sort and cumsum values
  sorted_cumsum <- base::cumsum(base::sort(standarised))

  ##### Restore the original elements order
  ordered_cumsum <- sorted_cumsum[ base::names(standarised) ]

  ##### Perform range standardization between 0 and 1, otherwise the negative values are summed up
  standarised_cumsum <- standardization(ordered_cumsum)

  ##### Clean the space and return output
  rm(x, standarised, sorted_cumsum, ordered_cumsum)
  return( standarised_cumsum )
}

##### Perform range standardization between 0 and 1 (for the cumulative sums)
standardization <- function(x) c(x-rapportools::min(x))/(rapportools::max(x)-rapportools::min(x))
