#' Generate CDF plot for selected gene
#'
#' Generates cumulative distribution function (CDF) plot for selected gene.
#' If option "addBoxPlot" = TRUE, then generate additional boxplot
#' below to show the data variance for selected gene in individual groups
#'
#' @param gene Gene of interest.
#' @param data Input data.
#' @param targets Target data set.
#' @param sampleName Sample name.
#' @param int_cancer Internal cancer group.
#' @param ext_cancer External cancer group.
#' @param comp_cancer Complete cancer group.
#' @param add_cancer Used for reordering groups.
#' @param addBoxPlot Add boc plot (boolean).
#' @param scaling Gene-wise or sample-wise.
#' @param report_dir Output directory.
#'
#' @importFrom magrittr %>%
#' @return Cumulative distribution function (CDF) plot for selected gene.
#' @export
#'
cdfPlot <- function(gene, data, targets, sampleName, int_cancer, ext_cancer, comp_cancer, add_cancer = NULL, addBoxPlot = FALSE, scaling = "gene-wise", report_dir) {
  ##### Remove the internal reference cohort data if the patient samples origins from other tissue. Of note, the internal reference cohort was only used to process the in-house data (including the investigated patient sample) and to correct batch-effects
  if (comp_cancer != int_cancer) {
    targets <- targets[!targets$Target %in% int_cancer, ]
    data <- data[, rownames(targets)]
  }

  ##### Initiate lists with stats for each group
  targets.list <- unique(targets$Target)
  group.z <- vector("list", length(targets.list))
  names(group.z) <- targets.list

  ##### .... and for selected gene
  group.z.gene <- vector("list", length(targets.list))
  names(group.z.gene) <- targets.list

  ##### Get expression-related stats for each group
  ##### ... from gene-wise approach
  if (scaling == "gene-wise") {
    ##### Get stats for each group
    gene.data <- data[gene, , drop = FALSE]
    group.z.gene <- exprGroupsStats_geneWise(gene.data, targets)[[1]]

    ##### ... and for each sample in individual groups
    gene.stats <- exprGroupsStats_geneWise(gene.data, targets)[[2]]

    for (group in targets.list) {
      group.z[[group]] <- cbind(t(gene.stats[[group]]$median), t(gene.stats[[group]]$z), t(gene.stats[[group]]$q), t(gene.stats[[group]]$cum))
      group.z[[group]] <- as.data.frame(group.z[[group]])
      colnames(group.z[[group]]) <- c("median", "z", "quantile", "cum")
    }

    group.z[[sampleName]] <- do.call("rbind", group.z)

    ##### ... or from group-wise approach
  } else {
    group.z[[sampleName]] <- exprGroupStats_groupWise(data, targets, sampleName)
    group.z[[ext_cancer]] <- exprGroupStats_groupWise(data, targets, ext_cancer)

    ##### Extract expression for selected genes
    group.z.gene[[sampleName]] <- group.z[[sampleName]][rownames(group.z[[sampleName]]) %in% gene, ]
    group.z.gene[[ext_cancer]] <- group.z[[ext_cancer]][rownames(group.z[[ext_cancer]]) %in% gene, ]

    ##### Add info for internal cohort
    if (comp_cancer == int_cancer) {
      group.z[[int_cancer]] <- exprGroupStats_groupWise(data, targets, int_cancer)
      group.z.gene[[int_cancer]] <- group.z[[int_cancer]][rownames(group.z[[int_cancer]]) %in% gene, ]
    }

    ##### Add info for additional cancer type is specified
    if (!is.null(add_cancer)) {
      group.z[[add_cancer]] <- exprGroupStats_groupWise(data, targets, add_cancer)
      group.z.gene[[add_cancer]] <- group.z[[add_cancer]][rownames(group.z[[add_cancer]]) %in% gene, ]
    }
  }

  ##### Generate box-plot for selected gene
  if (addBoxPlot) {
    ##### Perform Z-score transformation of the median expression values
    if (scaling == "gene-wise") {
      data.z <- t(scale(t(data)))
    } else {
      data.z <- scale(data, scale = FALSE)
    }

    targets$Target[targets$Target == sampleName] <- "Patient"
    gene.expr.df <- data.frame(targets$Target, data.z[gene, ])
    colnames(gene.expr.df) <- c("Group", "Expression")

    ##### Reorder groups
    if (!is.null(add_cancer)) {
      gene.expr.df$Group <- factor(gene.expr.df$Group, levels = c(add_cancer, ext_cancer, int_cancer, "Patient"))
      group.colours <- c("forestgreen", "cornflowerblue", "red", "black")
    } else {
      gene.expr.df$Group <- factor(gene.expr.df$Group, levels = c(ext_cancer, int_cancer, "Patient"))
      group.colours <- c("cornflowerblue", "red", "black")
    }

    p2 <- plotly::plot_ly(gene.expr.df, x = ~Expression, color = ~Group, type = "box", jitter = 0.3, pointpos = 0, boxpoints = "all", colors = group.colours, opacity = 0.5, orientation = "h", width = 800, height = 400, showlegend = FALSE)
  }

  ##### Generate interactive CDF plot with plotly
  ##### Include the internal reference cohort in the plot
  if (comp_cancer == int_cancer) {
    p1 <- plotly::plot_ly(group.z[[sampleName]], x = ~z, color = I("black"), width = 700, height = 200) %>%
      ##### Add sample data
      plotly::add_markers(
        y = group.z.gene[[sampleName]]$quantile, x = group.z.gene[[sampleName]]$z,
        text = rownames(group.z.gene[[sampleName]]),
        name = "Patient",
        marker = list(size = 12, color = "black"),
        showlegend = TRUE
      ) %>%
      plotly::add_lines(
        y = group.z[[sampleName]]$quantile, x = group.z[[sampleName]]$z,
        line = list(color = "grey"),
        text = rownames(group.z[[sampleName]]),
        name = "Patient", showlegend = FALSE
      ) %>%
      ##### Add int_cancer data
      plotly::add_markers(
        y = group.z.gene[[int_cancer]]$quantile, x = group.z.gene[[int_cancer]]$z,
        text = rownames(group.z.gene[[int_cancer]]),
        name = int_cancer,
        marker = list(size = 12, opacity = 0.5, color = "red"),
        showlegend = TRUE
      ) %>%
      plotly::add_lines(
        y = group.z[[int_cancer]]$quantile, x = group.z[[int_cancer]]$z, opacity = 0.5,
        line = list(color = "red", dash = "dash"),
        text = rownames(group.z[[int_cancer]]),
        name = int_cancer, showlegend = FALSE
      ) %>%
      ##### Add ext_cancer data
      plotly::add_markers(
        y = group.z.gene[[ext_cancer]]$quantile, x = group.z.gene[[ext_cancer]]$z,
        text = rownames(group.z.gene[[ext_cancer]]),
        name = ext_cancer,
        marker = list(size = 12, opacity = 0.5, color = "cornflowerblue"),
        showlegend = TRUE
      ) %>%
      plotly::add_lines(
        y = group.z[[ext_cancer]]$quantile, x = group.z[[ext_cancer]]$z, opacity = 0.5,
        line = list(color = "cornflowerblue", dash = "dash"),
        text = rownames(group.z[[ext_cancer]]),
        name = ext_cancer, showlegend = FALSE
      ) %>%
      ##### Add quantile lines
      plotly::add_lines(
        y = seq(0, 100, 10), x = rep(stats::quantile(group.z[[sampleName]]$z)[2], 11), opacity = 0.5,
        line = list(color = "gray", dash = "dash"),
        name = "Q1", showlegend = FALSE
      ) %>%
      plotly::add_lines(
        y = seq(0, 100, 10), x = rep(stats::quantile(group.z[[sampleName]]$z)[3], 11), opacity = 0.5,
        line = list(color = "gray", dash = "dash"),
        name = "Q2", showlegend = FALSE
      ) %>%
      plotly::add_lines(
        y = seq(0, 100, 10), x = rep(stats::quantile(group.z[[sampleName]]$z)[4], 11), opacity = 0.5,
        line = list(color = "gray", dash = "dash"),
        name = "Q3", showlegend = FALSE
      ) %>%
      plotly::layout(
        title = gene, xaxis = list(title = "mRNA expression (Z-score)", zeroline = FALSE, range = c(min(group.z[[sampleName]]$z) - 1.5, max(group.z[[sampleName]]$z) + 1.5)),
        yaxis = list(title = "Percentile"),
        legend = list(orientation = "v", x = 0.02, y = 1, bgcolor = "white")
      )

    ##### Skip the internal reference cohort in the plot
  } else {
    p1 <- plotly::plot_ly(group.z[[sampleName]], x = ~z, color = I("black"), width = 700, height = 200) %>%
      ##### Add sample data
      plotly::add_markers(
        y = group.z.gene[[sampleName]]$quantile, x = group.z.gene[[sampleName]]$z,
        text = rownames(group.z.gene[[sampleName]]),
        name = "Patient",
        marker = list(size = 12, color = "black"),
        showlegend = TRUE
      ) %>%
      plotly::add_lines(
        y = group.z[[sampleName]]$quantile, x = group.z[[sampleName]]$z,
        line = list(color = "grey"),
        text = rownames(group.z[[sampleName]]),
        name = "Patient", showlegend = FALSE
      ) %>%
      ##### Add ext_cancer data
      plotly::add_markers(
        y = group.z.gene[[ext_cancer]]$quantile, x = group.z.gene[[ext_cancer]]$z,
        text = rownames(group.z.gene[[ext_cancer]]),
        name = ext_cancer,
        marker = list(size = 12, opacity = 0.5, color = "cornflowerblue"),
        showlegend = TRUE
      ) %>%
      plotly::add_lines(
        y = group.z[[ext_cancer]]$quantile, x = group.z[[ext_cancer]]$z, opacity = 0.5,
        line = list(color = "cornflowerblue", dash = "dash"),
        text = rownames(group.z[[ext_cancer]]),
        name = ext_cancer, showlegend = FALSE
      ) %>%
      ##### Add quantile lines
      plotly::add_lines(
        y = seq(0, 1, 0.1), x = rep(stats::quantile(group.z[[sampleName]]$z)[2], 11), opacity = 0.5,
        line = list(color = "gray", dash = "dash"),
        name = "Q1", showlegend = FALSE
      ) %>%
      plotly::add_lines(
        y = seq(0, 1, 0.1), x = rep(stats::quantile(group.z[[sampleName]]$z)[3], 11), opacity = 0.5,
        line = list(color = "gray", dash = "dash"),
        name = "Q2", showlegend = FALSE
      ) %>%
      plotly::add_lines(
        y = seq(0, 1, 0.1), x = rep(stats::quantile(group.z[[sampleName]]$z)[4], 11), opacity = 0.5,
        line = list(color = "gray", dash = "dash"),
        name = "Q3", showlegend = FALSE
      ) %>%
      plotly::layout(
        title = gene, xaxis = list(title = "mRNA expression (Z-score)", zeroline = FALSE, range = c(min(group.z[[sampleName]]$z) - 1.5, max(group.z[[sampleName]]$z) + 1.5)),
        yaxis = list(title = "Percentile"),
        legend = list(orientation = "v", x = 0.02, y = 1, bgcolor = "white")
      )
  }

  ##### Combine CDF plot with boxplot if this option is selected
  if (addBoxPlot) {
    p1_2 <- plotly::subplot(p1, p2, nrows = 2, shareX = TRUE, shareY = FALSE, titleY = TRUE, heights = c(0.7, 0.3)) %>%
      plotly::layout(
        xaxis = list(title = "mRNA expression (Z-score)", zeroline = FALSE, range = c(min(group.z[[sampleName]]$z) - 1.5, max(group.z[[sampleName]]$z) + 1.5)),
        yaxis = list(title = "Percentile"),
        legend = list(orientation = "v", x = 0.02, y = 1, bgcolor = "white"),
        yaxis2 = list(title = ""), xaxis2 = list(title = paste0(gene, " mRNA expression (Z-score)")), margin = list(l = 50, r = 50, b = 50, t = 50, pad = 4), autosize = FALSE,
        showlegend = TRUE, showlegend2 = FALSE
      )

    return(p1_2)
  } else {
    return(p1)
  }
  ##### Clean the space
  rm(gene, targets, data, sampleName, targets.list, group.z, group.z.gene, gene.data, gene.stats, data.z, gene.expr.df, group.colours)

  #### Clear plots to free up some memory
  if (!is.null(grDevices::dev.list())) invisible(grDevices::dev.off())
}
