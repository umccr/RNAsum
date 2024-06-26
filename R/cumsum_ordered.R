#' Calculate cumulative sum
#'
#' Calculating cumulative sum while keeping the original data order.
#'
#' @param x Input.
#' @return Cumulative sum while keeping the order of input data.
#' @export
cumsum_ordered <- function(x) {

  ##### Perform range standardization between 0 and 1, otherwise the negative values are summed up
  standarised <- standardization(x)

  ##### Sort and cumsum values
  sorted_cumsum <- base::cumsum(base::sort(standarised))

  ##### Restore the original elements order
  ordered_cumsum <- sorted_cumsum[ base::names(standarised) ]

  ##### Perform range standardization between 0 and 1, otherwise the negative values are summed up
  standarised_cumsum <- standardization(ordered_cumsum)

  return( standarised_cumsum )
}

##### Perform range standardization between 0 and 1 (for the cumulative sums)
standardization <- function(x) c(x-base::min(x))/(base::max(x)-base::min(x))
