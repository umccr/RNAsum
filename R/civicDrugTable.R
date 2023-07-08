#' Generate table with drugs targeting selected set of genes
#'
#' Generates table with drugs targeting selected set of genes using info from
#' CIViC database (https://civicdb.org/).
#'
#' @param genes Genes of interest.
#' @param civic_var_summaries CIViC variant summaries.
#' @param civic_clin_evid CIViC clinical evidence.
#' @param evid_type Evidence type.
#' @param var_type Variant type.
#'
#' @return Table with drugs targeting selected set of genes.
#' @export
civicDrugTable <- function(genes, civic_var_summaries, civic_clin_evid, evid_type = "Predictive", var_type = NULL) {
  ##### Initialize data frame to the about drug-target info from CIViC
  drug.info <- stats::setNames(
    data.frame(matrix(ncol = 18, nrow = 0)),
    c(
      "Gene", "Variant", "variant_types", "drugs", "nct_ids", "evidence_level",
      "evidence_type", "evidence_direction", "clinical_significance", "rating",
      "civic_actionability_score", "Disease", "phenotypes", "pubmed_id",
      "variant_origin", "representative_transcript", "representative_transcript2",
      "last_review_date"
    )
  )

  evid_levels <- list(
    "A" = "A: Validated association",
    "B" = "B: Clinical evidence",
    "C" = "C: Case study",
    "D" = "D: Preclinical evidence",
    "E" = "E: Inferential association"
  )

  ##### Loop through each gene and check if they are druggable
  for (gene in genes) {
    ##### Get summary info about druggable genes
    if (gene %in% civic_clin_evid$gene) {
      ##### Extract info about all reported variants's clinical evidence for queried gene
      clin.evid.info <- civic_clin_evid[civic_clin_evid$gene == gene, ]

      ##### Use more descriptive evidence level info
      for (level in unique(clin.evid.info$evidence_level)) {
        clin.evid.info$evidence_level[clin.evid.info$evidence_level == level] <- evid_levels[[level]]
      }

      ##### Subset table to include only variants with the evidence type of interest
      clin.evid.info <- clin.evid.info |>
        dplyr::filter(.data$evidence_type == evid_type)

      if (nrow(clin.evid.info) > 0) {
        ##### Provide link to CIViC clinical evidence summary
        # TODO (PD): update the CIViC files; URL is outdated and points to insecure link.
        clin.evid.info$drugs <- paste0("<a href='", clin.evid.info$evidence_civic_url, "' target='_blank'>", clin.evid.info$drugs, "</a>")

        ##### Provide link to CIViC clinical evidence summary
        clin.evid.info$evidence_type <- paste0("<a href='", clin.evid.info$evidence_civic_url, "' target='_blank'>", clin.evid.info$evidence_type, "</a>")

        ##### Provide link to CIViC gene summary
        clin.evid.info$gene_civic_url <- paste0("<a href='", clin.evid.info$gene_civic_url, "' target='_blank'>", gene, "</a>")
        names(clin.evid.info)[names(clin.evid.info) == "gene_civic_url"] <- "Gene"

        ##### Provide link to CIViC variants summary
        clin.evid.info$variant_civic_url <- paste0("<a href='", clin.evid.info$variant_civic_url, "' target='_blank'>", clin.evid.info$variant, "</a>")
        names(clin.evid.info)[names(clin.evid.info) == "variant_civic_url"] <- "Variant"

        ##### Provide link to ClinicalTrials.gov variants summary based on NCT IDs
        for (nct_id in clin.evid.info$nct_ids) {
          if (!is.na(nct_id)) {
            ##### Deal with multiple NCT IDs (separated by comma)
            nct_id_url <- gsub(
              " '", "'",
              paste(
                gsub(
                  "/ ", "/",
                  paste(
                    "<a href='https://clinicaltrials.gov/ct2/show/",
                    unlist(strsplit(nct_id, split = ",", fixed = TRUE)),
                    "' target='_blank'>",
                    unlist(strsplit(nct_id, split = ",", fixed = TRUE)),
                    "</a>"
                  )
                ),
                collapse = ", "
              )
            )
            clin.evid.info$nct_ids[clin.evid.info$nct_ids == nct_id] <- nct_id_url
          }
        }

        ##### Provide link to PubMed variants summary
        clin.evid.info$pubmed_id <- paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed/", clin.evid.info$pubmed_id, "' target='_blank'>", clin.evid.info$pubmed_id, "</a>")

        ##### Provide link to Disease Ontology
        clin.evid.info$doid <- paste0("<a href='http://www.disease-ontology.org/?id=DOID:", clin.evid.info$doid, "' target='_blank'>", clin.evid.info$disease, "</a>")
        names(clin.evid.info)[names(clin.evid.info) == "doid"] <- "Disease"

        ##### Extract info about all variants it that gene
        var.info <- civic_var_summaries[civic_var_summaries$gene == gene, ] |>
          dplyr::select("variant", "variant_types", "civic_actionability_score") |>
          dplyr::mutate(
            variant_types = gsub("_", " ", .data$variant_types),
            variant_types = gsub(",", ", ", .data$variant_types)
          )

        ##### Merge about all variants it that gene and clinical evidence info
        clin.evid.info <- merge(clin.evid.info, var.info, by = "variant", all.x = TRUE)

        ##### Filter drug matching info depending on the variant type
        var_type.keep <- NULL

        ##### Remove entries containing "EXPRESSION", "AMPLIFICATION", "DELETION", "METHYLATION", "WILD TYPE", "FUSION", "COPY", "REARRANGEMENT", "PHOSPHORYLATION", "TRANSCRIPT", "GAIN", "LOSS"
        if (!is.null(var_type) && var_type == "mutation") {
          var_type.keep <- c(var_type.keep, grep("EXPRESSION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("AMPLIFICATION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("DELETION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("METHYLATION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("WILD TYPE", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("FUSION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("REARRANGEMENT", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("PHOSPHORYLATION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("COPY", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("TRANSCRIPT", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("GAIN", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("LOSS", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))

          clin.evid.info <- clin.evid.info[-c(unique(var_type.keep)), ]

          ##### Keep only entries containing "EXPRESSION", "FUSION", "TRANSCRIPT", "ALTERATION"
        } else if (!is.null(var_type) && var_type == "expression") {
          var_type.keep <- c(var_type.keep, grep("EXPRESSION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("FUSION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("TRANSCRIPT", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("ALTERATION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))

          clin.evid.info <- clin.evid.info[c(unique(var_type.keep)), ]

          ##### Keep only entries containing "FUSION", "ALTERATION", "[gene]-", "-[gene]"
        } else if (!is.null(var_type) && var_type == "fusion") {
          var_type.keep <- c(var_type.keep, grep("FUSION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep(paste0(gene, "-"), clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep(paste0("-", gene), clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("ALTERATION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))

          clin.evid.info <- clin.evid.info[c(unique(var_type.keep)), ]

          ##### Keep only entries containing "AMPLIFICATION", "COPY", "GAIN", "ALTERATION"
        } else if (!is.null(var_type) && var_type == "copy_gain") {
          var_type.keep <- c(var_type.keep, grep("AMPLIFICATION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("COPY", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("GAIN", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("ALTERATION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))

          clin.evid.info <- clin.evid.info[c(unique(var_type.keep)), ]

          ##### Keep only entries containing "DELETION", "COPY", "LOSS", "ALTERATION"
        } else if (!is.null(var_type) && var_type == "copy_loss") {
          var_type.keep <- c(var_type.keep, grep("DELETION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("COPY", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("LOSS", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))
          var_type.keep <- c(var_type.keep, grep("ALTERATION", clin.evid.info$variant, invert = FALSE, ignore.case = TRUE))

          clin.evid.info <- clin.evid.info[c(unique(var_type.keep)), ]
        }
      }

      if (nrow(clin.evid.info) > 0) {
        ##### Subset table to include only most important info
        clin.evid.info <- clin.evid.info[, names(drug.info)]

        ##### Add drugs info for subsequent gene
        drug.info <- rbind(drug.info, clin.evid.info)
      }
    }
  }

  ##### Use more friendly column names for the table
  names(drug.info) <- c("Gene", "Variant", "Variant type", "Drugs", "Clinical trials", "Evidence level", "Evidence type", "Evidence direction", "Clinical significance", "Trust rating", "Actionability score", "Disease", "Phenotypes", "PubMed ID", "Variant origin", "Representative transcript", "Representative transcript 2", "Review date")

  ##### Limit the info to fewer columns
  drug.info <- drug.info[, c("Gene", "Variant", "Variant type", "Drugs", "Clinical trials", "Evidence level", "Evidence direction", "Clinical significance", "Trust rating", "Actionability score", "Disease", "Phenotypes", "PubMed ID", "Representative transcript", "Representative transcript 2")]

  ##### Generate a table
  dt.table <- DT::datatable(
    data = drug.info, filter = "none",
    rownames = FALSE, extensions = c("Buttons", "Scroller"),
    options = list(
      pageLength = 10, dom = "Bfrtip",
      buttons = c("excel", "csv", "pdf", "copy", "colvis"),
      scrollX = TRUE, deferRender = TRUE, scrollY = "167px", scroller = TRUE
    ), width = 800,
    caption = htmltools::tags$caption(style = "caption-side: top; text-align: left; color:grey; font-size:100% ;"),
    escape = FALSE
  ) |>
    DT::formatStyle(columns = names(drug.info), `font-size` = "12px", "text-align" = "center") |>
    ##### Colour cells according to evidence level and trust rating
    DT::formatStyle(
      columns = "Evidence level",
      backgroundColor = DT::styleEqual(
        c("A: Validated association", "B: Clinical evidence", "C: Case study", "D: Preclinical evidence", "E: Inferential association"),
        c("mediumseagreen", "deepskyblue", "mediumpurple", "darkorange", "coral")
      )
    ) |>
    DT::formatStyle(
      columns = "Trust rating",
      backgroundColor = DT::styleEqual(c(1:5), c("coral", "azure", "lightskyblue", "palegreen", "mediumseagreen"))
    )

  return(list(dt.table, drug.info))
}
