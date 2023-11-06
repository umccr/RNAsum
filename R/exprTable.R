#' Generate table with coloured cells indicating expression values for selected genes
#'
#' Generates table with coloured cells indicating expression values for selected genes.
#'
#' @param genes Selected genes.
#' @param keep_all keep all rows
#' @param data Input data.
#' @param cn_data Copy number data.
#' @param sv_data SV data.
#' @param cn_decrease Order of the CN data.
#' @param targets Target groups.
#' @param sampleName Sample name.
#' @param int_cancer Internal cancer group.
#' @param ext_cancer External cancer group.
#' @param comp_cancer Complete cancer group.
#' @param add_cancer Additonal cancer type.
#' @param genes_annot Geners annotation.
#' @param oncokb_annot OncoKb annotation.
#' @param cancer_genes Cancer genes.
#' @param mut_annot Mutation annotation.
#' @param fusion_genes Fusion genes.
#' @param type Type.
#' @param civic_clin_evid civic_clin_evid tibble from reference genes.
#' @param scaling Scaling
#'
#' @importFrom dplyr %>%
#' @return Table with coloured cells indicating expression values for selected genes
#' @export
exprTable <- function(data = NULL, genes = NULL, keep_all = FALSE, cn_data = NULL, sv_data = NULL,
                      cn_decrease = TRUE, targets, sampleName, int_cancer, ext_cancer,
                      comp_cancer, add_cancer = NULL, genes_annot = NULL,
                      oncokb_annot = NULL, cancer_genes = NULL, mut_annot = NULL,
                      fusion_genes = NULL, type = "z", scaling = "gene-wise",
                      civic_clin_evid = NULL) {
  assertthat::assert_that(
    !is.null(genes), !is.null(data), !is.null(civic_clin_evid),
    is.matrix(data), is.data.frame(targets),
    "Target" %in% colnames(targets),
    type %in% c("perc", "z"),
    is.character(genes),
    scaling %in% c("gene-wise", "group-wise")
  )
  ##### Check which of the selected genes are not present in the expression data
  genes.absent <- genes[!genes %in% rownames(data)]

  ##### Initiate dataframe for expression median values in each group
  targets.list <- unique(targets$Target)
  group.z <- as.data.frame(matrix(NA, ncol = length(targets.list), nrow = nrow(data)))
  colnames(group.z) <- targets.list
  rownames(group.z) <- rownames(data)

  ##### Perform scaling gene-wise
  if (scaling == "gene-wise") {
    ##### Calculate z-score for each group
    group.stats <- exprGroupsStats_geneWise(data, targets)[[1]]

    ##### Make sure to include only genes for which Z-scores were calculated (genes with SD = 0 across all samples will give NA)
    group.z <- group.z[rownames(group.z) %in% rownames(group.stats[[targets.list[1]]]), ]

    #### Present expression data as percentiles or z-score values (default)
    for (group in targets.list) {
      if (type == "perc") {
        group.z[, group] <- round(group.stats[[group]]$quantile, digits = 1)
      } else {
        group.z[, group] <- round(group.stats[[group]]$z, digits = 2)
      }
    }

    ##### Perform scaling group-wise
  } else {
    for (group in targets.list) {
      ##### Calculate z-score for each group
      group.stats <- exprGroupStats_groupWise(data[rownames(group.z), ], targets, group)
      group.stats <- group.stats[order(rownames(group.stats)), ]

      #### Present expression data as percentiles or z-score values (default)
      if (type == "perc") {
        group.z[, group] <- round(group.stats$quantile, digits = 1)
      } else {
        group.z[, group] <- round(group.stats$z, digits = 2)
      }
    }
  }

  ##### If additional cancer type is defined then remove it from the data
  if (!is.null(add_cancer)) {
    group.z <- group.z[, !names(group.z) %in% add_cancer]
    targets <- targets[!targets$Target %in% add_cancer, ]
    targets.list <- targets.list[!targets.list %in% add_cancer]
  }

  ##### Compute Z-scores sd for each gene across groups
  group.z <- base::cbind(
    group.z,
    SD = round(matrixStats::rowSds(as.matrix(group.z)), digits = 2)
  )

  ##### Calculate Z-score differences between investigated sample and median values in the cancer group of interest
  group.z <- base::cbind(
    group.z,
    Diff = round((group.z[, sampleName] - group.z[, comp_cancer]), digits = 2)
  )

  ##### Add NAs for genes that are absent in the expression matrix. In the "Patient vs [comp_cancer]" columns provide "0"s to facilitate interactive sorting the table. These will appear in blank cells in the table
  if (length(genes.absent) > 0) {
    NAs.df <- data.frame(matrix(NA, ncol = ncol(group.z), nrow = length(genes.absent)))
    names(NAs.df) <- names(group.z)
    rownames(NAs.df) <- genes.absent
    NAs.df[names(NAs.df) %in% "Diff"] <- 0
    group.z <- rbind(group.z, NAs.df)
  }

  ##### Change sample ID to "Patient" for better visualisation
  names(group.z)[names(group.z) == sampleName] <- "Patient"
  targets.list[targets.list == sampleName] <- "Patient"

  ##### Reorder groups
  ##### Add "Gene" column to facilitate adding annotations
  group.z <- group.z |>
    dplyr::select(dplyr::all_of(c(ext_cancer, int_cancer)), "Patient", "SD", "Diff") |>
    dplyr::mutate(Gene = rownames(group.z))

  ##### Add genes annotation
  if (!is.null(genes_annot)) {
    ##### Remove rows with duplicated gene symbols
    if ("SYMBOL" %in% colnames(genes_annot)) {
      genes_annot <- genes_annot[!duplicated(genes_annot[["SYMBOL"]]), ]
    }

    ##### Merge the dataframe with groups median expression values and gene annotations
    group.z <- merge(genes_annot, group.z, by.x = "SYMBOL", by.y = "Gene", all = TRUE, sort = FALSE)
    names(group.z) <- sub("SYMBOL", "Gene", names(group.z))
  }

  ##### Define colours for cells background for each group and the patient vs [comp_cancer] difference
  brks_clrs1 <- brks_clrs(targets.list = targets.list, group.z = group.z, step1 = 0.0005)

  ##### Subset the expression data to include only the user-defined genes
  group.z <- group.z[group.z$Gene %in% genes, ]

  #### Add variants information to the expression table - if exists. Note, "TIER" and "CONSEQUENCE" columns are required
  if (!is.null(mut_annot) && "TIER" %in% colnames(mut_annot) && length(genes) > 0) {
    mut_annot <- mut_annot |> dplyr::filter(.data$SYMBOL %in% genes)

    #### keep only varaints that has the lowest tier value. Multiple varaints detected in same gene but with higher tier will be added to additional column "CONSEQUENCE_OTHER". Applies to the ones that may have multiple mutations and hence tiers
    ##### First, create a list of genes to store multiple variants
    mut_consequence <- vector("list", length(unique(mut_annot$SYMBOL)))
    mut_consequence <- stats::setNames(mut_consequence, unique(mut_annot$SYMBOL))

    ##### Record all varaints detected in individual genes
    if (nrow(mut_annot) > 0) {
      for (i in 1:nrow(mut_annot)) {
        mut_consequence[[mut_annot$SYMBOL[i]]] <- unique(c(mut_consequence[[mut_annot$SYMBOL[i]]], mut_annot$CONSEQUENCE[i]))
      }

      mut_annot$CONSEQUENCE_OTHER <- "-"
    }

    ##### Remove the first elements since these variant consequences will be reported as the "canonical" CONSEQUENCE
    mut_consequence <- lapply(mut_consequence, function(x) x[-1])

    ##### Order variant entires based on tier info, to make sure that the varaints with the lowest tier are reported first
    mut_annot <- mut_annot[order(mut_annot$TIER), ]

    ##### Remove rows with duplicated gene symbols
    mut_annot <- mut_annot[!duplicated(mut_annot$SYMBOL), ]

    ##### Add other provided variants consequences for individual genes
    for (gene in mut_annot[["SYMBOL"]]) {
      if (length(mut_consequence[[gene]]) > 0) {
        mut_annot$CONSEQUENCE_OTHER[match(gene, mut_annot$SYMBOL)] <- mut_consequence[[gene]]
      }
    }

    #### merge the variants information with the dataframe
    group.z <- merge(group.z, mut_annot, by.x = "Gene", by.y = "SYMBOL", all = TRUE, sort = FALSE)
  }

  ##### Add CN data if provided
  if (!is.null(cn_data)) {
    ##### Get the position of "Diff" column
    col_idx <- grep("Diff", names(group.z), fixed = TRUE)

    ##### Now place the CN data after the "Diff" column
    if (length(genes) > 0) {
      group.z <- group.z |>
        dplyr::left_join(cn_data, by = "Gene") |>
        dplyr::mutate(CN = round(.data$CN, digits = 2)) |>
        dplyr::relocate("CN", .after = col_idx) |>
        dplyr::rename("Patient (CN)" = "CN")
      cn_range <- base::range(group.z[["Patient (CN)"]], na.rm = TRUE)
    } else {
      group.z <- tibble::add_column(group.z, "", .after = col_idx)
      colnames(group.z)[col_idx + 1] <- "Patient (CN)"
      cn_range <- 0
    }
  }

  ##### Add structural variants results from MANTA
  if (!is.null(sv_data) && length(genes) > 0) {
    ##### NOTE: when merging per-gene expression data with SV data from MANTA the "gene" column is used since multiple entries are possible for one gene in MANTA output
    group.z <- merge(group.z, sv_data, by.x = "Gene", by.y = "Genes", all = TRUE, sort = FALSE)
  }

  ##### Add info about known fusion genes
  if (!is.null(fusion_genes) && length(genes) > 0) {
    group.z$Fusion_gene <- NA
    group.z$Fusion_gene[group.z$Gene %in% fusion_genes] <- "Yes"
  }

  ##### Add cancer gene resources info
  if (!is.null(cancer_genes) && length(genes) > 0) {
    group.z <- merge(group.z, cancer_genes, by.x = "Gene", by.y = "Gene", all = TRUE, sort = FALSE)
  }

  ##### Include only queried genes
  group.z <- group.z[group.z$Gene %in% genes, ]
  group.z$SYMBOL <- group.z$Gene

  ##### Add links to external gene annotation resourses
  if (length(genes) > 0) {
    # create mini tbls for urls
    vicc1 <- tibble::tibble(Gene = genes) |>
      dplyr::mutate(url_vicc = glue::glue("<a href='https://search.cancervariants.org/#{.data$Gene}'>VICC</a>"))
    civic1 <- civic_clin_evid |>
      dplyr::distinct(.data$gene, .data$gene_civic_url) |>
      dplyr::mutate(gene_civic_url = glue::glue("<a href='{.data$gene_civic_url}'>CIViC</a>")) |>
      dplyr::select(Gene = "gene", url_civic = "gene_civic_url")
    onco1 <- tibble::tibble(Gene = NA, url_oncokb = NA)
    if (!is.null(oncokb_annot)) {
      onco1 <- oncokb_annot |>
        dplyr::filter(.data$OncoKB == "Yes") |>
        dplyr::mutate(url_oncokb = glue::glue("<a href='http://oncokb.org/#/gene/{.data$Gene}'>OncoKB</a>")) |>
        dplyr::select("Gene", "url_oncokb")
    }


    group.z <- group.z |>
      dplyr::left_join(vicc1, by = "Gene") |>
      dplyr::left_join(onco1, by = "Gene") |>
      dplyr::left_join(civic1, by = "Gene") |>
      tidyr::unite(col = "External resources", "url_vicc", "url_oncokb", "url_civic", sep = ", ", na.rm = TRUE, remove = TRUE) |>
      dplyr::relocate("External resources", .after = .data$Diff) |>
      dplyr::mutate(Gene = glue::glue("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene={.data$Gene}'>{.data$Gene}</a>"))
    if ("ENSEMBL" %in% colnames(group.z)) {
      group.z <- group.z |>
        dplyr::mutate(ENSEMBL = glue::glue("<a href='http://ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={.data$ENSEMBL}'>{.data$ENSEMBL}</a>"))
    }
  }

  ##### Order the data by CN values (to allow filtering based on CN information) and then by the highest absolute values for Patient vs [comp_cancer] difference (to allow filtering based on z-score differences)
  if (!is.null(cn_data) && length(genes) > 0) {
    ##### Get the position of "Patient (CN)" column
    col_idx <- grep("Patient (CN)", names(group.z), fixed = TRUE)
    group.z <- group.z[order(abs(group.z[, "Diff"]), decreasing = TRUE), ]
    group.z <- group.z[order(group.z[, col_idx], decreasing = cn_decrease), ]

    ##### Order the data by increasing TIER category (to allow filtering based on tier information) and then by the highest absolute values for "Diff" difference (to allow filtering based on z-score differences)
  } else if (!is.null(mut_annot) && length(genes) > 0) {
    group.z <- group.z[order(abs(group.z[, "Diff"]), decreasing = TRUE), ]
    group.z <- group.z[order(group.z$TIER), ]

    ##### Order the data by MANTA increasing Tier (to prioritise SVs, based on https://github.com/AstraZeneca-NGS/simple_sv_annotation/blob/master/simple_sv_annotation.py), event type and then by the highest absolute values for Patient vs [comp_cancer] difference
  } else if (!is.null(sv_data) && length(genes) > 0) {
    group.z <- group.z[order(abs(group.z[["Diff"]]), decreasing = TRUE), ]
    group.z <- group.z[order(group.z[["fusion_genes"]], decreasing = TRUE), ]
    group.z <- group.z[order(group.z[["Tier"]]), ]

    ##### Otherwise order table by the highest absolute values for Patient vs [comp_cancer] difference
  } else if (length(genes) > 0) {
    group.z <- group.z[order(abs(group.z[, "Diff"]), decreasing = TRUE), ]
  }

  ##### Remove the internal reference cohort column if the patient samples origins from other tissue. Of note, the internal reference cohort was only used to process the in-house data (including the investigated patient sample) and to correct batch-effects
  if (comp_cancer != int_cancer) {
    group.z <- group.z[, !names(group.z) %in% int_cancer]
    targets.list[match(int_cancer, targets.list)] <- "Patient"

    ##### Get the position of "Diff" column
    diff_col_idx <- grep("Diff", names(group.z), fixed = TRUE)
  } else {
    ##### Get the position of "Diff" column
    diff_col_idx <- grep("Diff", names(group.z), fixed = TRUE)
    names(group.z)[match("Diff", names(group.z))] <- paste0("Patient vs ", comp_cancer)
  }

  ##### Limit the ordered table to maximum of 2000 entries if "keep_all" is set to FALSE (default)
  if (nrow(group.z) > 2000 && !keep_all) {
    group.z <- group.z[1:2000, ]
  }

  ##### Define table height
  if (nrow(group.z) == 2) {
    table_height <- 230
    scrollY <- "67px"
  } else {
    scrollY <- "167px"
    table_height <- 318
  }

  ##### Generate a table with genes annotations and coloured expression values in each group
  dt.table <- group.z |>
    dplyr::select(-c("SYMBOL", "SD")) |>
    DT::datatable(
      filter = "none", rownames = FALSE, extensions = c("Buttons", "Scroller"),
      options = list(
        pageLength = 10, dom = "Bfrtip", buttons = c("excel", "csv", "pdf", "copy", "colvis"),
        scrollX = TRUE, scrollCollapse = TRUE, deferRender = TRUE,
        scrollY = scrollY, scroller = TRUE
      ),
      width = 800, height = table_height,
      caption = htmltools::tags$caption(style = "caption-side: top; text-align: left; color:grey; font-size:100% ;"),
      escape = FALSE
    ) |>
    DT::formatStyle(
      columns = names(group.z)[!names(group.z) %in% c("SYMBOL", "SD")],
      `font-size` = "12px", "text-align" = "center"
    ) |>
    ##### Colour cells according to the expression values quantiles in each group
    DT::formatStyle(
      columns = targets.list[1],
      backgroundColor = DT::styleInterval(
        brks_clrs1[["brks"]][[targets.list[1]]],
        brks_clrs1[["clrs"]][[targets.list[1]]]
      )
    ) |>
    DT::formatStyle(
      columns = targets.list[2],
      backgroundColor = DT::styleInterval(
        brks_clrs1[["brks"]][[targets.list[2]]],
        brks_clrs1[["clrs"]][[targets.list[2]]]
      )
    ) |>
    DT::formatStyle(
      columns = targets.list[3],
      backgroundColor = DT::styleInterval(
        brks_clrs1[["brks"]][[targets.list[3]]],
        brks_clrs1[["clrs"]][[targets.list[3]]]
      )
    ) |>
    DT::formatStyle(
      columns = names(group.z)[diff_col_idx],
      backgroundColor = DT::styleInterval(
        brks_clrs1[["brks"]][["Diff"]],
        brks_clrs1[["clrs"]][["Diff"]]
      )
    )

  if (!is.null(cn_data)) {
    dt.table <- dt.table |>
      DT::formatStyle(
        columns = "Patient (CN)",
        background = DT::styleColorBar(cn_range, "lightblue"),
        backgroundSize = "98% 88%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )

    ##### Generate a table with genes annotations and coloured expression values in each group
  }
  return(list(dt.table, group.z))
}

# Define colours for cells background for each group and the patient vs [comp_cancer] difference
brks_clrs <- function(targets.list, group.z, step1 = 0.0005) {
  ##### Initiate dataframe for expression median values in each group
  probs1 <- seq(.05, .95, step1)
  brks <- list()
  clrs <- list()
  clrs_pos.q <- round(seq(255, 150, length.out = length(probs1) / 2 + 1.5), 0) %>%
    {
      paste0("rgb(255,", ., ",", ., ")")
    }
  clrs_neg.q <- rev(round(seq(255, 150, length.out = length(probs1) / 2 - 0.5), 0)) %>%
    {
      paste0("rgb(", ., ",", ., ",", "255)")
    }
  clrs.q <- c(clrs_neg.q, clrs_pos.q)
  for (group in c(targets.list, "Diff")) {
    brks.q <- stats::quantile(group.z[, group], probs = probs1, na.rm = TRUE)
    brks.q <- c(-Inf, brks.q)
    assertthat::assert_that(length(brks.q) == length(clrs.q))
    # grab top colour for duplicated brks
    res <- tibble::tibble(brks.q = unname(brks.q), clrs.q = clrs.q) |>
      dplyr::group_by(brks.q) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup()
    brks[[group]] <- res[["brks.q"]][-1] # remove -Inf
    clrs[[group]] <- res[["clrs.q"]]
  }
  list(
    brks = brks,
    clrs = clrs
  )
}
