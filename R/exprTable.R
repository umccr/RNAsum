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
#' @param ext_links External links.
#' @param type Type.
#' @param civic_clin_evid civic_clin_evid tibble from reference genes.
#' @param scaling Scaling
#'
#' @importFrom magrittr %>%
#' @return Table with coloured cells indicating expression values for selected genes
#' @export
exprTable <- function(data = NULL, genes = NULL, keep_all = FALSE, cn_data = NULL, sv_data = NULL,
                      cn_decrease = TRUE, targets, sampleName, int_cancer, ext_cancer,
                      comp_cancer, add_cancer = NULL, genes_annot = NULL,
                      oncokb_annot = NULL, cancer_genes = NULL, mut_annot = NULL,
                      fusion_genes = NULL, ext_links = FALSE, type = "z", scaling = "gene-wise",
                      civic_clin_evid = NULL) {
  assertthat::assert_that(
    !is.null(genes), !is.null(data), !is.null(civic_clin_evid),
    is.matrix(data), is.data.frame(targets),
    "Target" %in% colnames(targets),
    type %in% c("perc", "z"),
    is.character(genes),
    scaling %in% c("gene-wise", "group-wise")
  )
  # TODO (PD): refactor this function
  ##### Check which of the selected genes are not present in the expression data
  genes.absent <- genes[genes %!in% rownames(data)]

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
    group.z <- group.z[, names(group.z) %!in% add_cancer]
    targets <- targets[targets$Target %!in% add_cancer, ]
    targets.list <- targets.list[targets.list %!in% add_cancer]
  }

  ##### Compute Z-scores sd for each gene across groups
  group.z <- cbind(group.z, round(matrixStats::rowSds(as.matrix(group.z)), digits = 2))
  names(group.z)[ncol(group.z)] <- "SD"

  ##### Calculate Z-score differences between investigated sample and median values in the cancer group of interest
  group.z <- cbind(group.z, round((group.z[, sampleName] - group.z[, comp_cancer]), digits = 2))
  names(group.z)[ncol(group.z)] <- "Diff"

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
  group.z <- cbind(group.z[, c(ext_cancer, int_cancer, "Patient")], group.z[, c("SD", "Diff")])

  ##### Add "Gene" column to facilitate adding annotations
  group.z$Gene <- rownames(group.z)

  ##### Add genes annotation
  if (!is.null(genes_annot)) {
    ##### Remove rows with duplicated gene symbols
    if ("SYMBOL" %in% names(genes_annot)) {
      genes_annot <- genes_annot[!duplicated(genes_annot$SYMBOL), ]
    }

    ##### Merge the dataframe with groups median expression values and gene annotations
    group.z <- merge(genes_annot, group.z, by.x = "SYMBOL", by.y = "Gene", all = TRUE, sort = FALSE)
    names(group.z) <- gsub("SYMBOL", "Gene", names(group.z))
  }

  ##### Define colours for cells background for each group and the patient vs [comp_cancer] difference
  ##### Initiate dataframe for expression median values in each group
  brks.q <- as.data.frame(matrix(NA, ncol = length(targets.list), nrow = length(seq(.05, .95, .0005))))
  colnames(brks.q) <- targets.list
  clrs.q <- as.data.frame(matrix(NA, ncol = length(targets.list), nrow = length(seq(.05, .95, .0005)) + 1))
  colnames(clrs.q) <- targets.list

  for (group in c(targets.list, "Diff")) {
    brks.q[[group]] <- stats::quantile(group.z[, group], probs = seq(.05, .95, .0005), na.rm = TRUE)

    clrs_pos.q <- round(seq(255, 150, length.out = length(brks.q[[group]]) / 2 + 1.5), 0) %>%
      {
        paste0("rgb(255,", ., ",", ., ")")
      }
    clrs_neg.q <- rev(round(seq(255, 150, length.out = length(brks.q[[group]]) / 2 - 0.5), 0)) %>%
      {
        paste0("rgb(", ., ",", ., ",", "255)")
      }
    clrs.q[[group]] <- c(clrs_neg.q, clrs_pos.q)
  }

  ##### Subset the expression data to include only the user-defined genes
  group.z <- group.z[group.z$Gene %in% genes, ]

  #### Add variants information to the expression table - if exists. Note, "TIER" and "CONSEQUENCE" columns are required
  if (!is.null(mut_annot) && "TIER" %in% colnames(mut_annot) && length(genes) > 0) {
    mut_annot <- mut_annot[mut_annot$SYMBOL %in% genes, ]

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
    rownames(mut_annot) <- mut_annot$SYMBOL

    ##### Add other provided variants consequences for individual genes
    for (gene in rownames(mut_annot)) {
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
      group.z <- tibble::add_column(group.z, round(cn_data[group.z$Gene, "CN"], digits = 2), .after = col_idx)
      colnames(group.z)[col_idx + 1] <- "Patient (CN)"
      cn_range <- base::range(group.z[, "Patient (CN)"], na.rm = TRUE)
    } else {
      group.z <- tibble::add_column(group.z, "", .after = col_idx)
      colnames(group.z)[col_idx + 1] <- "Patient (CN)"
      cn_range <- 0
    }
  }

  ##### Add structural variants results from MANTA
  if (!is.null(sv_data) && length(genes) > 0) {
    ##### NOTE: when merging per-gene exprssion data with SV data from MANTA the "gene" column is used since multiple entires are possible for one gene in MANTA output
    group.z <- merge(group.z, sv_data, by.x = "Gene", by.y = "Gene", all = TRUE, sort = FALSE)
  }

  ##### Add info about known fusion genes
  if (!is.null(fusion_genes) && length(genes) > 0) {
    group.z$Fusion_gene <- NA
    group.z$Fusion_gene[group.z$Gene %in% fusion_genes] <- "Yes"
  }

  ##### Add cancer gene resources info
  if (!is.null(cancer_genes) && length(genes) > 0) {
    group.z <- merge(group.z, cancer_genes, by.x = "Gene", by.y = "row.names", all = TRUE, sort = FALSE)
  }

  ##### Include only queried genes
  group.z <- group.z[group.z$Gene %in% genes, ]
  group.z$SYMBOL <- group.z$Gene

  ##### Add links to external gene annotation resourses
  if (ext_links && length(genes) > 0) {
    ##### Place the external links after the "Diff" column
    ##### Get the position of "Diff" column
    col_idx <- grep("Diff", names(group.z), fixed = TRUE)
    group.z <- tibble::add_column(group.z, NA, .after = col_idx)
    names(group.z)[col_idx + 1] <- "ext_links"

    # create mini tbls for urls
    vicc1 <- tibble::tibble(Gene = genes) |>
      dplyr::mutate(url_vicc = glue::glue("<a href='https://search.cancervariants.org/#{.data$Gene}' target='_blank'>VICC</a>"))
    civic1 <- civic_clin_evid |>
      dplyr::distinct(.data$gene, .data$gene_civic_url) |>
      dplyr::mutate(gene_civic_url = glue::glue("<a href='{.data$gene_civic_url}' target='_blank'>CIViC</a>")) |>
      dplyr::select(Gene = "gene", url_civic = "gene_civic_url")
    onco1 <- tibble::tibble(Gene = NA, url_oncokb = NA)
    if (!is.null(oncokb_annot)) {
      onco1 <- oncokb_annot |>
        dplyr::filter(.data$OncoKB == "Yes") |>
        dplyr::mutate(url_oncokb = glue::glue("<a href='http://oncokb.org/#/gene/{.data$Gene}' target='_blank'>OncoKB</a>")) |>
        dplyr::select("Gene", "url_oncokb")
    }


    x <- group.z |> tibble::as_tibble()
    xy <- x |>
      dplyr::select(-ext_links) |>
      dplyr::left_join(vicc1, by = "Gene") |>
      dplyr::left_join(onco1, by = "Gene") |>
      dplyr::left_join(civic1, by = "Gene") |>
      tidyr::unite(col = "ext_links", url_vicc, url_oncokb, url_civic, sep = ", ", na.rm = TRUE, remove = TRUE) |>
      dplyr::relocate(ext_links, .after = Diff) |>
      dplyr::mutate(gene_url = glue::glue("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene={.data$Gene}' target='_blank'>{.data$Gene}</a>"))
    if ("ENSEMBL" %in% colnames(xy)) {
      xy <- xy |>
        dplyr::mutate(ENSEMBL = glue::glue("<a href='http://ensembl.org/Homo_sapiens/Gene/Summary?db=core;g={.data$ENSEMBL}' target='_blank'>{.data$ENSEMBL}</a>"))
    }
    xy <- xy |>
      dplyr::relocate("gene_url") |>
      dplyr::select(-"Gene") |>
      dplyr::rename(
        "External resources" = "ext_links",
        "Gene" = "gene_url"
      )

    for (gene in genes) {
      # For every gene (~18K)
      ##### Provide link to VICC meta-knowledgebase ( https://search.cancervariants.org )
      group.z$ext_links[group.z$Gene == gene] <- paste0("<a href='https://search.cancervariants.org/#", gene, "' target='_blank'>VICC</a>")

      ##### Provide link to OncoKB
      if (!is.null(oncokb_annot)) {
        if ((gene %in% oncokb_annot$Gene) && ((oncokb_annot |> dplyr::filter(.data$Gene == gene) |> dplyr::pull("OncoKB")) == "Yes")) {
          group.z$ext_links[group.z$Gene == gene] <- paste(
            group.z$ext_links[group.z$Gene == gene],
            paste0(
              "<a href='http://oncokb.org/#/gene/",
              gene,
              "' target='_blank'>OncoKB</a>"
            ),
            sep = ", "
          )
        }
      }

      ##### Provide link to CIViC database druggable genes ( https://civicdb.org )
      assertthat::assert_that("gene" %in% colnames(civic_clin_evid))
      if (gene %in% civic_clin_evid[["gene"]]) {
        group.z$ext_links[group.z$Gene == gene] <- paste(
          group.z$ext_links[group.z$Gene == gene],
          paste0(
            "<a href='",
            unique(civic_clin_evid[civic_clin_evid[["gene"]] == gene, "gene_civic_url"]),
            "' target='_blank'>CIViC</a>"
          ),
          sep = ", "
        )
      }
    }

    names(group.z) <- gsub("ext_links", "External resources", names(group.z))
  }

  ##### Attach links to GeneCards and Ensembl (if provided). Here we assume that gene names are
  for (gene in genes) {
    if ("ENSEMBL" %in% names(group.z)) {
      if (!is.na(group.z$ENSEMBL[group.z$Gene == gene])) {
        group.z$ENSEMBL[group.z$Gene == gene] <- paste0(
          "<a href='http://ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",
          group.z$ENSEMBL[group.z$Gene == gene],
          "' target='_blank'>",
          group.z$ENSEMBL[group.z$Gene == gene],
          "</a>"
        )
      }
    }

    group.z$Gene[group.z$Gene == gene] <- paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "' target='_blank'>", gene, "</a>")
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
    group.z <- group.z[order(abs(group.z[, "Diff"]), decreasing = TRUE), ]
    group.z <- group.z[order(group.z$"Fusion genes", decreasing = TRUE), ]
    group.z <- group.z[order(group.z$Tier), ]

    ##### Otherwise order table by the highest absolute values for Patient vs [comp_cancer] difference
  } else if (length(genes) > 0) {
    group.z <- group.z[order(abs(group.z[, "Diff"]), decreasing = TRUE), ]
  }

  ##### Remove the internal reference cohort column if the patient samples origins from other tissue. Of note, the internal reference cohort was only used to process the in-house data (including the investigated patient sample) and to correct batch-effects
  if (comp_cancer != int_cancer) {
    group.z <- group.z[, names(group.z) %!in% int_cancer]
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
  dt.table <- DT::datatable(
    data = group.z[, names(group.z) %!in% c("SYMBOL", "SD")],
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
    DT::formatStyle(columns = names(group.z)[names(group.z) %!in% c("SYMBOL", "SD")], `font-size` = "12px", "text-align" = "center") |>
    ##### Colour cells according to the expression values quantiles in each group
    DT::formatStyle(
      columns = targets.list[1],
      backgroundColor = DT::styleInterval(brks.q[[targets.list[1]]], clrs.q[[targets.list[1]]])
    ) |>
    DT::formatStyle(
      columns = targets.list[2],
      backgroundColor = DT::styleInterval(brks.q[[targets.list[2]]], clrs.q[[targets.list[2]]])
    ) |>
    DT::formatStyle(
      columns = targets.list[3],
      backgroundColor = DT::styleInterval(brks.q[[targets.list[3]]], clrs.q[[targets.list[3]]])
    ) |>
    DT::formatStyle(
      columns = names(group.z)[diff_col_idx],
      backgroundColor = DT::styleInterval(brks.q[["Diff"]], clrs.q[["Diff"]])
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
