##### Calculate group-wise median, sd, quantiles and cumulative fractions for expression data from specific sample group
#' Calculate group-wise median, standard deviation, quantiles and cumulative fractions.
#'
#' @param data Input data.
#' @param targets Target groups.
#' @param target Defined biological group.
#'
#' @return Group-wise median, standard deviation, quantiles and cumulative fractions for expression data from specific sample group.
#' @export
#'

exprGroupStats_groupWise <- function(data, targets, target) {

  ##### Subset data for defined biological group
  data.group <- data[, targets$Target %in% target ]

  ##### For groups with > 1 sample get the median and standard deviation for each gene
  if ( !is.null(ncol(data.group)) )  {

    data.group.median <- matrixStats::rowMedians(data.group)
    names(data.group.median) <- rownames(data.group)
    data.group.median <- sort(data.group.median)
    data.group.sd <- matrixStats::rowSds(data.group)

  } else {
    data.group.median <- sort(data.group)
    data.group.sd <- rep( NA, length(data.group))
  }

  ##### Make sure the median and sd vectors have the same gene order
  names(data.group.sd) <- rownames(data.group)
  data.group.sd <- data.group.sd[names(data.group.median)]

  ##### Convert a expression values into corresponding percentiles
  data.group.q <- perc_rank(data.group.median)

  ##### Perform range standardization between 0 and 1 (for the cumulative sums), otherwise the negative values are summed up
  data.group.s <- sort(standardization(data.group.median))

  ##### Calculate cumulative sums and perform range standardization between 0 and 1
  data.group.cum <- standardization(cumsum(data.group.s))

  ##### Perform Z-score transformation of the median expression values
  data.group.z <- scale(data.group.median, scale = FALSE)

  ##### Organise the data into data frame
  data.group.df <- as.data.frame(cbind( data.group.median, data.group.sd, data.group.z, data.group.q, data.group.cum))
  names(data.group.df) <- c("median", "sd", "z", "quantile", "cum")

  ##### Clean the space and return output
  rm(data, targets, target, data.group, data.group.median, data.group.sd, data.group.q, data.group.s, data.group.cum, data.group.z)
  return( data.group.df )
}
