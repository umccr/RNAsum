##### Generate boxplot presenting expression profiles for selected set of genes
glanceExprPlot <- function(genes, data, targets, sampleName, int_cancer, ext_cancer, comp_cancer, add_cancer = NULL, hexcode, type = "z", sort = "diff", scaling = "gene-wise", report_dir) {

  if ( comp_cancer != int_cancer ) {
    targets <- targets[ targets$Target %!in% int_cancer, ]
    data <- data[ ,rownames(targets) ]
  }

  ##### Perform Z-score transformation of the median expression values
  if ( scaling == "gene-wise" ) {

    data.z <- t(scale(t(data)))
    y_title <- "mRNA expression (Z-score)"

    if ( type == "perc" ) {
      ##### Convert a expression values into corresponding percentiles
      data.z <- t(apply(data.z, 1, perc_rank))
      y_title <- "mRNA expression (percentile)"
    }

  } else {
    data.z <- scale(data, scale = FALSE)

    if ( type == "perc" ) {
      ##### Convert a expression values into corresponding percentiles
      data.z <- t(apply(data.z, 1, perc_rank))
    }
  }

  targets$Target[ targets$Target==sampleName ] <- "Patient"

  ##### Make sure that all genes are present in the expression matrix
  genes <- genes[ genes %in% rownames(data.z) ]

  ##### Genes sorting for visualisation
  ##### Sort genes by the greatest difference between the patient and the "comp_cancer" cohort
  if ( sort == "diff" ) {
    comp_cancer.medians <- rowMedians( data.z[ genes ,targets$Target==comp_cancer ] )
    names(comp_cancer.medians) <- genes
    comp_cancer.medians.diff <- comp_cancer.medians - data.z[ genes ,targets$Target=="Patient" ]
    genes <- genes[ order(comp_cancer.medians.diff) ]

    ##### Sort genes alphabetically
  } else if (sort == "alphabetically") {
    genes <- genes[ order(genes) ]
  }

  ##### Prepare dataframe for plotly
  gene.expr.df <- NULL

  for ( gene in genes ) {
    gene.expr.df <- rbind(gene.expr.df, data.frame(gene, targets$Target, data.z[gene, ]))
  }
  colnames(gene.expr.df) <- c("Gene", "Group", "Expression")

  ##### Reorder groups
  if ( !is.null(add_cancer) ) {
    gene.expr.df$Group <- factor(gene.expr.df$Group, levels=c("Patient", int_cancer, ext_cancer, add_cancer))
    group.colours <- c(I("black"), "red", "cornflowerblue", "forestgreen")

  } else {
    gene.expr.df$Group <- factor(gene.expr.df$Group, levels=c("Patient", int_cancer, ext_cancer))
    group.colours <- c(I("black"), "red", "cornflowerblue")
  }

  p <- plot_ly( gene.expr.df, x = ~Gene, y = ~Expression, color = ~Group, type = "box", colors = group.colours, opacity=0.3, showlegend = TRUE, width = 800, height = 400 ) %>%
    add_markers(x = ~Gene[ gene.expr.df$Group %in% "Patient" ], y = ~Expression[ gene.expr.df$Group %in% "Patient" ], color = ~Group[ gene.expr.df$Group %in% "Patient" ], marker = list(size = 7), opacity=1, showlegend = FALSE) %>%

    layout(boxmode = "group", xaxis = list(title = ""), yaxis = list(title = y_title), legend = list( orientation = 'h', y = max(gene.expr.df$Expression), yancho = "top", bgcolor = "white"))

  ##### Create directory for "at glance" plots
  PlotsDir <- paste(report_dir, "glanceExprPlots", sep = "/")

  if ( !file.exists(PlotsDir) ) {
    dir.create(PlotsDir, recursive=TRUE)
  }

  ##### Save interactive plot as html file
  saveWidgetFix(p, file = paste(PlotsDir, paste0(hexcode, "_glance_expr_plot.", type, ".html"), sep = "/"))

  return( p )

  ##### Clean the space and return output
  rm(targets, data, sampleName, data.z, y_title, genes, comp_cancer.medians, comp_cancer.medians.diff, gene.expr.df, group.colours)

  #### Clear plots to free up some memory
  if(!is.null(dev.list())) invisible(dev.off())
}
