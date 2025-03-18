#' Calculate gene-wise median, sd, quantiles and cumulative fractions for expression data
#'
#' Calculates gene-wise median, sd, quantiles and cumulative fractions for expression data.
#'
#' @param data input data.
#' @param targets Target groups
#'
#' @return Gene-wise median, standard deviation, quantiles and cumulative fractions for expression data.
#' @export
exprGroupsStats_geneWise <- function(data, targets) {

  ##### Perform Z-score transformation of the expression values
  data.z <- t(apply(data, 1, scale, scale = TRUE))
  colnames(data.z) <- colnames(data)

  ##### Remove rows with potential NA's, which is due to SD = 0 across all samples
  data.z <- data.z[rowSums(!is.na(data.z)) > 0, , drop = FALSE]
  data <- data[ rownames(data) %in% rownames(data.z), , drop = FALSE]

  ##### Perform the gene-wise calculations across all groups
  ##### Convert a expression values into corresponding percentiles
  data.q <- t(apply(data, 1, perc_rank))

  ##### Calculate cumulative sums and perform range standardization between 0 and 1
  data.cum <- t(apply(data, 1, cumsum_ordered))
  colnames(data.cum) <- colnames(data.q)

  ##### Create lists with stats for each group and gene
  targets.list <- unique(targets$Target)
  group_stats.list <- vector("list", length(targets.list))
  names(group_stats.list) <- targets.list

  #### For each group...
  for ( group in targets.list ) {

    ##### For groups with > 1 sample get the median values for each gene
    if ( base::sum(c(targets$Target %in% group), na.rm = TRUE) > 1 && nrow(data) > 1 )  {

      ##### Extract the median expression values
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], matrixStats::rowMedians(data[ , colnames(data)[ targets$Target %in% group ] ]))

      ##### Extract the expression sd values
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], matrixStats::rowSds(data[ , colnames(data)[ targets$Target %in% group ] ]))

      ##### Extract the median Z-scores
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], matrixStats::rowMedians(data.z[ , colnames(data)[ targets$Target %in% group ] ]))

      ##### Extract the median percentiles
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], matrixStats::rowMedians(data.q[ , colnames(data)[ targets$Target %in% group ] ]))

      ##### Extract the cumulative fraction corresponding to the median Z-score
      ##### First, need to get the position of the Z-score nearest to the median Z-score, and then extract the cumulative value at this position
      data.z.median_pos <- apply(data.z, 1, nearest_position, stats::median(data.z[ , colnames(data)[ targets$Target %in% group, drop = FALSE ] ]))
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], data.cum[ data.z.median_pos ] )

      group_stats.list[[group]] <- as.data.frame(group_stats.list[[group]])
      names( group_stats.list[[group]] ) <- c("median", "sd", "z", "quantile", "cum")
      rownames( group_stats.list[[group]] ) <- rownames(data)

    } else if ( base::sum(c(targets$Target %in% group), na.rm = TRUE) > 1 && nrow(data) == 1 ) {

      ##### Extract the median expression values
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], stats::median(data[ , colnames(data)[ targets$Target %in% group, drop = FALSE ] ]))

      ##### Extract the expression sd values
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], stats::sd(data[ , colnames(data)[ targets$Target %in% group, drop = FALSE ] ]))

      ##### Extract the median Z-scores
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], stats::median(data.z[ , colnames(data)[ targets$Target %in% group, drop = FALSE ] ]))

      ##### Extract the median percentiles
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], stats::median(data.q[ , colnames(data)[ targets$Target %in% group, drop = FALSE ] ]))

      ##### Extract the cumulative fraction corresponding to the median Z-score
      ##### First, need to get the position of the Z-score nearest to the median Z-score, and then extract the cumulative value at this position
      data.z.median_pos <- nearest_position( data.z, stats::median(data.z[ , colnames(data)[ targets$Target %in% group, drop = FALSE ] ]))
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], data.cum[ data.z.median_pos ] )

      group_stats.list[[group]] <- as.data.frame(group_stats.list[[group]])
      names( group_stats.list[[group]] ) <- c("median", "sd", "z", "quantile", "cum")
      rownames( group_stats.list[[group]] ) <- rownames(data)

    } else {

      ##### Extract the median expression values
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], data[ , colnames(data)[ targets$Target %in% group ] ])

      ##### Extract the expression sd values
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], rep( NA, nrow(data)))

      ##### Extract the median Z-scores
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], data.z[ , colnames(data)[ targets$Target %in% group ] ])

      ##### Extract the median percentiles
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], data.q[ , colnames(data)[ targets$Target %in% group ] ])

      ##### Extract the median cumulative fraction
      group_stats.list[[group]] <- cbind(group_stats.list[[group]], data.cum[ , colnames(data)[ targets$Target %in% group ] ])

      group_stats.list[[group]] <- as.data.frame(group_stats.list[[group]])
      names( group_stats.list[[group]] ) <- c("median", "sd", "z", "quantile","cum")
      rownames( group_stats.list[[group]] ) <- rownames(data)
    }
  }

  ##### Finally, extract cumulative values for each gene within individual groups
  gene_stats.list <- vector("list", length(targets.list))
  names(gene_stats.list) <- targets.list

  #### For each group...
  for ( group in targets.list ) {

    ##### Extract per-gene expression values
    gene_stats.list[[group]]$median <- data[ , colnames(data)[ targets$Target %in% group ], drop = FALSE ]

    ##### Extract per-gene z-score values
    gene_stats.list[[group]]$z <- data.z[ , colnames(data.z)[ targets$Target %in% group ], drop = FALSE ]

    ##### Extract per-gene percentile values
    gene_stats.list[[group]]$q <- data.q[ , colnames(data.q)[ targets$Target %in% group ], drop = FALSE ]

    ##### Extract per-gene cumulative values
    gene_stats.list[[group]]$cum <- data.cum[ , colnames(data.cum)[ targets$Target %in% group ], drop = FALSE ]
  }

  ##### Clean the space and return output
  rm(data, targets, data.z, data.q, data.cum, targets.list, data.z.median_pos)
  return( list( group_stats.list, gene_stats.list) )
}
