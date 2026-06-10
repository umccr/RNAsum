#' Batch Effects Assessment Functions for RNAsum
#'
#' Functions to assess and visualize potential batch effects between clinical
#' samples and TCGA reference cohorts, particularly for protocol differences
#' like ribo-depletion vs poly-A selection.
#'
#' @name batch_assessment
NULL

#' Detect Gene ID Format
#'
#' Helper function to detect the format of gene identifiers.
#'
#' @param gene_ids Character vector of gene identifiers
#'
#' @return Character string describing the detected format
#' @keywords internal
detect_gene_id_format <- function(gene_ids) {
  if (length(gene_ids) == 0) {
    return("empty")
  }

  # Take a sample for analysis
  sample_ids <- head(gene_ids, min(1000, length(gene_ids)))

  # Check different patterns
  ensembl_count <- sum(grepl("^ENSG[0-9]{11}(\\.[0-9]+)?$", sample_ids))
  ensembl_simple <- sum(grepl("^ENSG", sample_ids))
  entrez_count <- sum(grepl("^[0-9]+$", sample_ids))
  symbol_count <- sum(grepl("^[A-Z][A-Z0-9_-]*$", sample_ids, ignore.case = FALSE))
  mixed_case_symbol <- sum(grepl("^[A-Za-z][A-Za-z0-9_-]*$", sample_ids))

  total_sample <- length(sample_ids)
  ensembl_pct <- ensembl_count / total_sample
  symbol_pct <- symbol_count / total_sample

  if (ensembl_pct > 0.8) {
    if (ensembl_count == ensembl_simple) {
      return("Ensembl gene IDs (with versions)")
    } else {
      return("Ensembl gene IDs")
    }
  } else if (symbol_pct > 0.7) {
    return("Gene symbols (HGNC-like)")
  } else if (entrez_count / total_sample > 0.8) {
    return("Entrez gene IDs")
  } else if (mixed_case_symbol / total_sample > 0.7) {
    return("Mixed-case gene symbols")
  } else {
    return(paste0("Mixed/Unknown (", ensembl_pct * 100, "% Ensembl, ", symbol_pct * 100, "% symbols)"))
  }
}

#' Built-in Cancer Gene Symbol to Ensembl ID Mapping
#'
#' Comprehensive mapping of cancer gene symbols to Ensembl gene IDs.
#' Used for automatic gene ID format conversion in batch assessment.
#'
#' @format A data frame with cancer gene symbol to Ensembl ID mappings
#' @keywords internal
.cancer_gene_mapping <- data.frame(
  symbol = c(
    "TP53", "BRCA1", "BRCA2", "KRAS", "EGFR", "PTEN", "VHL", "PIK3CA",
    "CTNNB1", "AR", "BRAF", "PTCH1", "PMS1", "MSH2", "MLH1", "MSH6",
    "IDH1", "IDH2", "KIT", "NF1", "MET", "PALB2", "ATM", "ATR",
    "ERBB2", "MYC", "SMAD4", "HOXB13", "CHEK1", "CHEK2", "CEBPA",
    "PIK3R1", "RB1", "TSC1", "TSC2", "E2F1", "SETD2", "FGFR3", "PPARG", "ZEB1",
    "APC", "CDKN2A", "NOTCH1", "GATA3", "CDH1", "RUNX1", "ARID1A", "KMT2D",
    "TET2", "DNMT3A", "ASXL1", "SF3B1", "SRSF2", "U2AF1", "ZRSR2", "EZH2",
    "KDM6A", "CREBBP", "EP300", "FBXW7", "NFE2L2", "KEAP1", "STK11", "CDKN1B",
    "FANCA", "FANCC", "FANCD2", "FANCF", "FANCG", "BRIP1", "RAD51C", "RAD51D",
    "MUTYH", "MSH3", "MSH6", "PMS2", "EPCAM", "BMPR1A", "SMAD4", "STK11",
    "PRKDC", "XRCC2", "XRCC3", "RAD50", "MRE11A", "NBN", "CHEK2", "MDM2",
    "CDKN2A", "CDKN2B", "CCND1", "CCNE1", "CDK4", "CDK6", "CCNA2", "CCNB1"
  ),
  ensembl_id = c(
    "ENSG00000141510", "ENSG00000012048", "ENSG00000139618", "ENSG00000133703",
    "ENSG00000146648", "ENSG00000171862", "ENSG00000134086", "ENSG00000198738",
    "ENSG00000184574", "ENSG00000169083", "ENSG00000157764", "ENSG00000073756",
    "ENSG00000115946", "ENSG00000095002", "ENSG00000076242", "ENSG00000119537",
    "ENSG00000138413", "ENSG00000182054", "ENSG00000134419", "ENSG00000108753",
    "ENSG00000183914", "ENSG00000165659", "ENSG00000143627", "ENSG00000149311",
    "ENSG00000141736", "ENSG00000136997", "ENSG00000141646", "ENSG00000174775",
    "ENSG00000115526", "ENSG00000183765", "ENSG00000164292", "ENSG00000121879",
    "ENSG00000139687", "ENSG00000089327", "ENSG00000103197", "ENSG00000101412",
    "ENSG00000151835", "ENSG00000068078", "ENSG00000109819", "ENSG00000169554",
    "ENSG00000134982", "ENSG00000147889", "ENSG00000148400", "ENSG00000107485",
    "ENSG00000039068", "ENSG00000159216", "ENSG00000117020", "ENSG00000167548",
    "ENSG00000168769", "ENSG00000119772", "ENSG00000163904", "ENSG00000115524",
    "ENSG00000161547", "ENSG00000159335", "ENSG00000106546", "ENSG00000181555",
    "ENSG00000147050", "ENSG00000005339", "ENSG00000188229", "ENSG00000109670",
    "ENSG00000111276", "ENSG00000079999", "ENSG00000118971", "ENSG00000111276",
    "ENSG00000163641", "ENSG00000134250", "ENSG00000169398", "ENSG00000144554",
    "ENSG00000138039", "ENSG00000103197", "ENSG00000164104", "ENSG00000164109",
    "ENSG00000132274", "ENSG00000095002", "ENSG00000119537", "ENSG00000122512",
    "ENSG00000198825", "ENSG00000165886", "ENSG00000141646", "ENSG00000118971",
    "ENSG00000166803", "ENSG00000196504", "ENSG00000196159", "ENSG00000114480",
    "ENSG00000109685", "ENSG00000060688", "ENSG00000183765", "ENSG00000135679",
    "ENSG00000147889", "ENSG00000147883", "ENSG00000123374", "ENSG00000105173",
    "ENSG00000134057", "ENSG00000065361", "ENSG00000123374", "ENSG00000105173"
  ),
  stringsAsFactors = FALSE
)

#' Convert Cancer Gene Symbols to Ensembl IDs
#'
#' Helper function to convert cancer gene symbols to Ensembl gene IDs
#' when sample data uses Ensembl format.
#'
#' @param gene_symbols Character vector of gene symbols
#' @return Character vector of Ensembl gene IDs (where mapping exists)
#' @keywords internal
convert_symbols_to_ensembl <- function(gene_symbols) {
  if (length(gene_symbols) == 0) {
    return(character(0))
  }

  # Find matches in mapping table
  mapping_subset <- .cancer_gene_mapping[.cancer_gene_mapping$symbol %in% gene_symbols, ]

  if (nrow(mapping_subset) == 0) {
    warning("No cancer gene symbols could be converted to Ensembl IDs")
    return(character(0))
  }

  # Return unique Ensembl IDs
  ensembl_ids <- unique(mapping_subset$ensembl_id)

  # Report conversion success
  n_converted <- length(ensembl_ids)
  n_total <- length(unique(gene_symbols))

  if (n_converted < n_total) {
    missing_symbols <- setdiff(gene_symbols, mapping_subset$symbol)
    if (length(missing_symbols) <= 5) {
      missing_str <- paste(missing_symbols, collapse = ", ")
    } else {
      missing_str <- paste(c(head(missing_symbols, 3), "..."), collapse = ", ")
    }
    message(sprintf("Converted %d/%d cancer gene symbols to Ensembl IDs. Missing: %s",
                   n_converted, n_total, missing_str))
  } else {
    message(sprintf("Successfully converted %d cancer gene symbols to Ensembl IDs", n_converted))
  }

  return(ensembl_ids)
}

#' Load Cancer Genes from RNAsum Databases
#'
#' Helper function to load cancer genes from available RNAsum databases.
#'
#' @param source Character, cancer gene source: "umccr", "oncokb", or "combined"
#'
#' @return Character vector of gene symbols
#' @keywords internal
load_cancer_genes <- function(source = "umccr") {

  # Validate source
  valid_sources <- c("umccr", "oncokb", "combined")
  if (!source %in% valid_sources) {
    stop("source must be one of: ", paste(valid_sources, collapse = ", "))
  }

  # First try to load from RNAsum package
  cancer_genes <- character()

  tryCatch({
    # Try to load reference gene annotations using get_refgenes
    if (requireNamespace("RNAsum", quietly = TRUE) && exists("get_refgenes")) {
      # Get reference genes (use NULL to load from default package files)
      ref_genes <- get_refgenes(p = NULL)

      if (source %in% c("umccr", "combined")) {
        # Load UMCCR cancer genes
        if ("genes_cancer" %in% names(ref_genes) && !is.null(ref_genes$genes_cancer)) {
          if ("Gene" %in% names(ref_genes$genes_cancer)) {
            umccr_genes <- ref_genes$genes_cancer$Gene
          } else if ("symbol" %in% names(ref_genes$genes_cancer)) {
            umccr_genes <- ref_genes$genes_cancer$symbol
          } else {
            # Try the first character column
            char_cols <- sapply(ref_genes$genes_cancer, is.character)
            if (any(char_cols)) {
              umccr_genes <- ref_genes$genes_cancer[[which(char_cols)[1]]]
            } else {
              umccr_genes <- character()
            }
          }
          cancer_genes <- c(cancer_genes, umccr_genes)
        }
      }

      if (source %in% c("oncokb", "combined")) {
        # Load OncoKB cancer genes
        if ("genes_oncokb" %in% names(ref_genes) && !is.null(ref_genes$genes_oncokb)) {
          if ("Hugo_Symbol" %in% names(ref_genes$genes_oncokb)) {
            oncokb_genes <- ref_genes$genes_oncokb$Hugo_Symbol
          } else if ("Gene" %in% names(ref_genes$genes_oncokb)) {
            oncokb_genes <- ref_genes$genes_oncokb$Gene
          } else {
            # Try the first character column
            char_cols <- sapply(ref_genes$genes_oncokb, is.character)
            if (any(char_cols)) {
              oncokb_genes <- ref_genes$genes_oncokb[[which(char_cols)[1]]]
            } else {
              oncokb_genes <- character()
            }
          }
          cancer_genes <- c(cancer_genes, oncokb_genes)
        }
      }
    } else {
      # Fallback to built-in cancer gene mapping when RNAsum functions not available
      message("RNAsum get_refgenes() not available, using built-in cancer gene mapping")
      cancer_genes <- .cancer_gene_mapping$symbol
    }

    # Remove duplicates and filter out empty values
    cancer_genes <- unique(cancer_genes[!is.na(cancer_genes) & cancer_genes != ""])

    if (length(cancer_genes) == 0) {
      # Final fallback using built-in mapping
      message("No cancer genes found from primary source, using built-in mapping")
      cancer_genes <- .cancer_gene_mapping$symbol
    }

    return(cancer_genes)

  }, error = function(e) {
    # Fallback to built-in cancer gene mapping on any error
    message("Error loading cancer genes from RNAsum, using built-in mapping: ", e$message)
    return(.cancer_gene_mapping$symbol)
  })
}

#' Assess Batch Effects Between Sample and Reference Cohort
#'
#' Performs comprehensive batch effect assessment including PCA analysis,
#' expression distribution comparisons, and batch severity metrics.
#'
#' @param sample_data Named vector of gene expression values for clinical sample
#' @param reference_data Matrix or data.frame of reference cohort expression
#'   (genes as rows, samples as columns)
#' @param gene_set_type Character, type of gene set to use for analysis:
#'   "top_n" (default), "cancer_genes", or "custom"
#' @param gene_subset Character vector of genes to focus analysis on. Required when
#'   gene_set_type = "custom"
#' @param n_genes Integer, number of top variable genes to use when gene_set_type = "top_n"
#' @param cancer_gene_source Character, source of cancer genes when gene_set_type = "cancer_genes":
#'   "umccr" (default), "oncokb", or "combined"
#' @param protocol_clinical Character, RNA-seq protocol for clinical sample
#'   (e.g., "ribo-depletion", "poly-A")
#' @param protocol_reference Character, RNA-seq protocol for reference cohort
#' @param output_dir Character, directory to save diagnostic plots
#' @param apply_log Logical, apply log2 transformation to data (default: TRUE)
#' @param apply_filter Logical, filter low-expressed genes (default: TRUE)
#' @param transform_method Character, transformation method: "CPM" (default), "TPM", "none"
#' @param norm_method Character, normalization method: "TMM" (default), "quantile", "none"
#' @param scaling_method Character, scaling method for Z-score transformation: "gene-wise" (default)
#'
#' @return List containing:
#'   \item{pca_results}{PCA analysis results}
#'   \item{batch_metrics}{Quantitative batch effect metrics}
#'   \item{recommendations}{Text recommendations for batch correction}
#'   \item{plots}{List of diagnostic plots}
#'   \item{gene_set_info}{Information about genes used in analysis}
#'
#' @export
assess_batch_effects <- function(sample_data,
                                  reference_data,
                                  gene_set_type = "top_n",
                                  gene_subset = NULL,
                                  n_genes = 2000,
                                  cancer_gene_source = "umccr",
                                  protocol_clinical = "unknown",
                                  protocol_reference = "TCGA_poly-A",
                                  output_dir = NULL,
                                  # Enhanced normalization parameters (User Method defaults)
                                  apply_log = TRUE,
                                  apply_filter = TRUE,
                                  transform_method = "CPM",
                                  norm_method = "TMM",
                                  scaling_method = "gene-wise") {

  # Load required libraries
  suppressMessages({
    require(ggplot2)
    require(dplyr)
    require(plotly)
    require(RColorBrewer)
  })

  # Input validation
  if (is.null(names(sample_data))) {
    stop("sample_data must be a named vector with gene identifiers")
  }

  # Validate gene_set_type
  valid_types <- c("top_n", "cancer_genes", "custom")
  if (!gene_set_type %in% valid_types) {
    stop("gene_set_type must be one of: ", paste(valid_types, collapse = ", "))
  }

  # Validate gene_subset for custom type
  if (gene_set_type == "custom" && is.null(gene_subset)) {
    stop("gene_subset must be provided when gene_set_type = 'custom'")
  }

  # Find common genes
  common_genes <- intersect(names(sample_data), rownames(reference_data))
  if (length(common_genes) < 100) {
    stop("Insufficient common genes between sample and reference (< 100)")
  }

  # Select genes for analysis based on gene_set_type
  gene_set_info <- list()

  if (gene_set_type == "top_n") {
    # Use top variable genes (original implementation)
    gene_vars <- apply(reference_data[common_genes, ], 1, var, na.rm = TRUE)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_genes, length(gene_vars))]
    gene_set_info <- list(
      type = "top_n",
      n_requested = n_genes,
      n_used = length(top_genes),
      description = paste("Top", length(top_genes), "variable genes")
    )

  } else if (gene_set_type == "cancer_genes") {
    # Load cancer genes from RNAsum databases with automatic format conversion
    tryCatch({
      cancer_genes <- load_cancer_genes(cancer_gene_source)

      # Detect gene ID format in sample data and convert cancer genes if needed
      sample_gene_format <- detect_gene_id_format(common_genes)
      cancer_gene_format <- detect_gene_id_format(cancer_genes)

      # Check if format conversion is needed
      if (grepl("Ensembl", sample_gene_format) && grepl("symbol", cancer_gene_format)) {
        message("Detected gene ID format mismatch: ")
        message("  Sample data: ", sample_gene_format)
        message("  Cancer genes: ", cancer_gene_format)
        message("  Converting cancer genes to Ensembl format...")

        # Convert cancer gene symbols to Ensembl IDs
        cancer_genes_converted <- convert_symbols_to_ensembl(cancer_genes)

        if (length(cancer_genes_converted) > 0) {
          cancer_genes <- cancer_genes_converted
          cancer_gene_format <- "Ensembl gene IDs (converted)"
        } else {
          warning("Failed to convert cancer gene symbols to Ensembl IDs")
        }
      }

      # Find intersection with converted genes
      top_genes <- intersect(cancer_genes, common_genes)

      # Handle cases where few cancer genes are found
      if (length(top_genes) < 100) {
        conversion_rate <- if(grepl("converted", cancer_gene_format)) {
          length(cancer_genes) / length(load_cancer_genes(cancer_gene_source))
        } else {
          1.0
        }

        msg <- paste0(
          "Limited cancer genes found in data (", length(top_genes), "/", length(cancer_genes), "). ",
          "Format compatibility:\n",
          "  • Sample data: ", sample_gene_format, " (", length(common_genes), " genes)\n",
          "  • Cancer genes: ", cancer_gene_format, " (", length(cancer_genes), " genes)\n"
        )

        if (conversion_rate < 0.5) {
          msg <- paste0(msg, "  • Low conversion rate suggests limited built-in mapping\n")
        }

        # Adaptive fallback strategy based on number of genes found
        if (length(top_genes) < 10) {
          warning(msg, "  Falling back to top variable genes for more robust analysis.")
          # Fallback to top variable genes
          gene_vars <- apply(reference_data[common_genes, ], 1, var, na.rm = TRUE)
          top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(1000, length(gene_vars))]
          gene_set_info$type <- "top_n_fallback"
          gene_set_info$description <- paste("Top variable genes (cancer genes insufficient:", length(intersect(cancer_genes, common_genes)), "genes)")
          gene_set_info$fallback_reason <- "Too few cancer genes found"
        } else {
          # Proceed with available cancer genes but warn user
          warning(msg, "  Proceeding with", length(top_genes), "available cancer genes.")
        }
      }

      gene_set_info <- list(
        type = "cancer_genes",
        source = cancer_gene_source,
        sample_format = sample_gene_format,
        cancer_format = cancer_gene_format,
        n_cancer_total = length(cancer_genes),
        n_used = length(top_genes),
        conversion_applied = grepl("converted", cancer_gene_format),
        description = paste("Cancer genes from", cancer_gene_source, "database:", length(top_genes), "genes")
      )

    }, error = function(e) {
      warning("Failed to load cancer genes: ", e$message, ". Falling back to top variable genes.")
      gene_vars <- apply(reference_data[common_genes, ], 1, var, na.rm = TRUE)
      top_genes <<- names(sort(gene_vars, decreasing = TRUE))[1:min(n_genes, length(gene_vars))]
      gene_set_info <<- list(
        type = "top_n_fallback",
        description = paste("Top", length(top_genes), "variable genes (cancer genes failed)")
      )
    })

  } else if (gene_set_type == "custom") {
    # Use user-provided gene subset
    top_genes <- intersect(gene_subset, common_genes)

    # Enhanced error handling for gene ID format mismatches
    if (length(top_genes) == 0) {
      # Check if this might be a gene ID format issue
      sample_gene_format <- detect_gene_id_format(common_genes)
      custom_gene_format <- detect_gene_id_format(gene_subset)

      error_msg <- paste0(
        "No genes from custom subset found in data (0/", length(gene_subset), "). ",
        "This may be due to gene ID format mismatch:\n",
        "  • Sample data format: ", sample_gene_format, "\n",
        "  • Custom genes format: ", custom_gene_format, "\n",
        "  • Suggestions:\n",
        "    - Use gene_set_type='cancer_genes' for cancer-related analysis\n",
        "    - Ensure custom genes match sample data ID format\n",
        "    - First few sample gene IDs: ", paste(head(names(sample_data), 3), collapse=", ")
      )
      stop(error_msg)
    }

    if (length(top_genes) < 50) {
      warning("Few genes in custom subset found in data (", length(top_genes), "), results may be unreliable. ",
              "Consider using gene_set_type='cancer_genes' or ensuring gene ID format compatibility.")
    }

    gene_set_info <- list(
      type = "custom",
      n_requested = length(gene_subset),
      n_used = length(top_genes),
      format_mismatch_warning = length(top_genes) < length(gene_subset) * 0.1,
      description = paste("Custom gene set:", length(top_genes), "genes")
    )
  }

  # Prepare data matrices
  ref_matrix <- reference_data[top_genes, ]
  sample_vec <- sample_data[top_genes]

  # Validate we have enough genes for analysis
  min_genes <- if(gene_set_type == "custom") 3 else 10  # More lenient for custom gene sets
  if (length(top_genes) < min_genes) {
    stop("Too few genes for batch assessment (", length(top_genes), "). Need at least ", min_genes, " genes for ", gene_set_type, " analysis.")
  }

  # Apply normalization transformations if requested
  if (apply_log || transform_method != "none" || norm_method != "none" || apply_filter) {
    cat("Applying data transformations...\n")

    # Step 1: Filtering (if requested)
    if (apply_filter) {
      # Basic filtering: remove very low expressed genes
      min_expression <- 0.1  # Less aggressive filtering for TPM data
      keep_genes <- names(sample_vec)[sample_vec > min_expression]
      common_keep <- intersect(keep_genes, rownames(ref_matrix))
      # Only filter if we retain a reasonable number of genes
      if (length(common_keep) > 1000) {
        sample_vec <- sample_vec[common_keep]
        ref_matrix <- ref_matrix[common_keep, ]
        cat("  Filtered to", length(common_keep), "genes with expression >", min_expression, "\n")
      } else {
        cat("  Skipping filtering - too few genes would remain (", length(common_keep), ")\n")
      }
    }

    # Step 2: Normalization method (TMM - before transformation)
    if (norm_method == "TMM") {
      if (requireNamespace("edgeR", quietly = TRUE) && length(sample_vec) > 0 && nrow(ref_matrix) > 0) {
        # Assume input data are TPM counts, convert to approximate counts for TMM
        # This is a simplified approach for demonstration
        combined_data <- cbind(ref_matrix, Clinical_Sample = sample_vec)
        # Convert TPM to pseudo-counts by scaling
        pseudo_counts <- round(combined_data * 1000)  # Scale factor

        # Check for valid counts
        if (all(dim(pseudo_counts) > 0) && sum(pseudo_counts, na.rm = TRUE) > 0) {
          dge <- edgeR::DGEList(counts = pseudo_counts)
          dge <- edgeR::calcNormFactors(dge, method = "TMM")

        # Step 3: Transform to CPM with TMM factors
        if (transform_method == "CPM") {
          normalized_data <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
          ref_matrix <- normalized_data[, colnames(ref_matrix)]
          sample_vec <- normalized_data[, "Clinical_Sample"]
          cat("  Applied TMM normalization + CPM transformation\n")
        } else {
          # Apply TMM factors to original data
          tmm_factors <- dge$samples$norm.factors
          ref_matrix <- sweep(ref_matrix, 2, tmm_factors[1:ncol(ref_matrix)], "/")
          sample_vec <- sample_vec / tmm_factors[ncol(combined_data)]
          cat("  Applied TMM normalization\n")
        }
        } else {
          cat("  Warning: Insufficient data for TMM normalization - skipping\n")
        }
      } else {
        cat("  Warning: edgeR package required for TMM normalization - skipping\n")
      }
    } else {
      # Step 3: Transformation method (CPM/TPM) without TMM
      if (transform_method == "CPM") {
        # For CPM transformation: counts per million
        sample_vec <- (sample_vec / sum(sample_vec, na.rm = TRUE)) * 1e6
        ref_matrix <- apply(ref_matrix, 2, function(x) (x / sum(x, na.rm = TRUE)) * 1e6)
        cat("  Applied CPM transformation\n")
      } else if (transform_method == "TPM") {
        cat("  Note: Input data assumed to be TPM already\n")
      }

      # Step 3b: Quantile normalization (alternative to TMM)
      if (norm_method == "quantile") {
        if (requireNamespace("preprocessCore", quietly = TRUE)) {
          combined_data <- cbind(ref_matrix, Clinical_Sample = sample_vec)
          normalized_data <- preprocessCore::normalize.quantiles(as.matrix(combined_data))
          dimnames(normalized_data) <- dimnames(combined_data)
          ref_matrix <- normalized_data[, colnames(ref_matrix)]
          sample_vec <- normalized_data[, "Clinical_Sample"]
          names(sample_vec) <- rownames(ref_matrix)
          cat("  Applied quantile normalization\n")
        } else {
          cat("  Warning: preprocessCore package required for quantile normalization - skipping\n")
        }
      }
    }

    # Step 4: Log transformation (after normalization and transformation)
    if (apply_log) {
      # Add small offset to avoid log of zero
      offset <- 0.25
      sample_vec <- log2(sample_vec + offset)
      ref_matrix <- log2(ref_matrix + offset)
      cat("  Applied log2 transformation (offset:", offset, ")\n")
    }
  }

  # 1. PCA Analysis
  pca_results <- perform_batch_pca(sample_vec, ref_matrix)

  # 2. Calculate batch effect metrics
  batch_metrics <- calculate_batch_metrics(sample_vec, ref_matrix)

  # 3. Generate diagnostic plots
  plots <- generate_batch_plots(sample_vec, ref_matrix, pca_results,
                                protocol_clinical, protocol_reference)

  # 4. Generate recommendations
  recommendations <- generate_batch_recommendations(batch_metrics,
                                                   protocol_clinical,
                                                   protocol_reference)

  # 5. Save plots if output directory specified
  if (!is.null(output_dir)) {
    save_batch_plots(plots, output_dir)
  }

  # Return results
  list(
    pca_results = pca_results,
    batch_metrics = batch_metrics,
    recommendations = recommendations,
    plots = plots,
    genes_analyzed = length(top_genes),
    gene_set_info = gene_set_info,
    protocols = list(clinical = protocol_clinical, reference = protocol_reference)
  )
}

#' Perform PCA Analysis for Batch Assessment
#'
#' @param sample_vec Named vector of sample expression values
#' @param ref_matrix Matrix of reference expression data
#'
#' @return List with PCA results and sample projection
#' @keywords internal
perform_batch_pca <- function(sample_vec, ref_matrix) {

  # Combine data for PCA
  combined_data <- cbind(ref_matrix, Clinical_Sample = sample_vec)

  # Remove genes with too many NAs
  complete_genes <- apply(combined_data, 1, function(x) sum(is.na(x)) < ncol(combined_data) * 0.1)
  pca_data <- combined_data[complete_genes, ]

  # Remove genes with constant or near-constant expression (causes PCA scaling issues)
  gene_vars <- apply(pca_data, 1, var, na.rm = TRUE)
  variable_genes <- gene_vars > 1e-10 & !is.na(gene_vars)  # Remove essentially constant genes
  pca_data <- pca_data[variable_genes, ]

  # Ensure we have enough genes for PCA
  if (nrow(pca_data) < 5) {
    warning("Too few variable genes for PCA (", nrow(pca_data), "). Using top variable genes from available data.")
    # Use all available genes if very few are left
    pca_data <- combined_data[complete_genes, ]
    if (nrow(pca_data) < 3) {
      stop("Insufficient genes for PCA analysis (", nrow(pca_data), " genes)")
    }
  }

  # Transpose for PCA (samples as rows)
  pca_input <- t(pca_data)

  # Handle missing values
  pca_input[is.na(pca_input)] <- rowMeans(pca_input, na.rm = TRUE)[col(pca_input)]

  # Final check for constant columns after all filtering
  col_vars <- apply(pca_input, 2, var, na.rm = TRUE)
  constant_cols <- is.na(col_vars) | col_vars < 1e-10
  if (any(constant_cols)) {
    cat("  Removing", sum(constant_cols), "constant expression genes before PCA\n")
    pca_input <- pca_input[, !constant_cols]
  }

  # Ensure we still have enough genes
  if (ncol(pca_input) < 3) {
    stop("Insufficient variable genes for PCA after filtering (", ncol(pca_input), " genes)")
  }

  # Perform PCA
  pca <- prcomp(pca_input, scale. = TRUE, center = TRUE)

  # Calculate sample distances from reference centroid
  ref_coords <- pca$x[rownames(pca$x) != "Clinical_Sample", 1:2]
  sample_coords <- pca$x["Clinical_Sample", 1:2]

  ref_centroid <- colMeans(ref_coords)
  distance_from_centroid <- sqrt(sum((sample_coords - ref_centroid)^2))

  # Calculate percentile of sample distance
  ref_distances <- sqrt(rowSums((ref_coords - rep(ref_centroid, each = nrow(ref_coords)))^2))
  distance_percentile <- mean(ref_distances < distance_from_centroid)

  list(
    pca = pca,
    sample_coords = sample_coords,
    ref_centroid = ref_centroid,
    distance_from_centroid = distance_from_centroid,
    distance_percentile = distance_percentile,
    variance_explained = summary(pca)$importance[2, 1:2]
  )
}

#' Calculate Quantitative Batch Effect Metrics
#'
#' @param sample_vec Named vector of sample expression values
#' @param ref_matrix Matrix of reference expression data
#'
#' @return List of batch effect metrics
#' @keywords internal
calculate_batch_metrics <- function(sample_vec, ref_matrix) {

  # 1. Expression distribution metrics
  ref_mean <- rowMeans(ref_matrix, na.rm = TRUE)
  ref_sd <- apply(ref_matrix, 1, sd, na.rm = TRUE)

  # Z-scores relative to reference
  z_scores <- (sample_vec - ref_mean) / ref_sd
  z_scores <- z_scores[!is.na(z_scores) & is.finite(z_scores)]

  # 2. Correlation with reference samples
  ref_cors <- apply(ref_matrix, 2, function(ref_sample) {
    common_genes <- intersect(names(sample_vec), names(ref_sample))
    cor(sample_vec[common_genes], ref_sample[common_genes], use = "complete.obs")
  })

  # 3. Rank correlation (Spearman)
  ref_rank_cors <- apply(ref_matrix, 2, function(ref_sample) {
    common_genes <- intersect(names(sample_vec), names(ref_sample))
    cor(sample_vec[common_genes], ref_sample[common_genes],
        method = "spearman", use = "complete.obs")
  })

  # 4. Expression level shift
  median_shift <- median(sample_vec - ref_mean, na.rm = TRUE)

  # 5. Variance inflation
  sample_var <- var(sample_vec, na.rm = TRUE)
  ref_var <- mean(apply(ref_matrix, 2, var, na.rm = TRUE), na.rm = TRUE)
  variance_ratio <- sample_var / ref_var

  list(
    z_score_mean = mean(abs(z_scores)),
    z_score_extreme = sum(abs(z_scores) > 2) / length(z_scores),
    correlation_median = median(ref_cors, na.rm = TRUE),
    correlation_min = min(ref_cors, na.rm = TRUE),
    rank_correlation_median = median(ref_rank_cors, na.rm = TRUE),
    median_expression_shift = median_shift,
    variance_ratio = variance_ratio,
    outlier_genes = sum(abs(z_scores) > 3)
  )
}

#' Generate Diagnostic Plots for Batch Assessment
#'
#' @param sample_vec Named vector of sample expression values
#' @param ref_matrix Matrix of reference expression data
#' @param pca_results PCA analysis results
#' @param protocol_clinical Clinical sample protocol
#' @param protocol_reference Reference cohort protocol
#'
#' @return List of ggplot objects
#' @keywords internal
generate_batch_plots <- function(sample_vec, ref_matrix, pca_results,
                                protocol_clinical, protocol_reference) {

  # 1. PCA plot
  pca_data <- data.frame(
    PC1 = pca_results$pca$x[, 1],
    PC2 = pca_results$pca$x[, 2],
    Sample_Type = c(rep("Reference", ncol(ref_matrix)), "Clinical")
  )

  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Sample_Type)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Clinical" = "red", "Reference" = "lightblue")) +
    labs(title = "PCA: Clinical Sample vs Reference Cohort",
         subtitle = paste("Clinical:", protocol_clinical, "| Reference:", protocol_reference),
         x = paste0("PC1 (", round(pca_results$variance_explained[1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(pca_results$variance_explained[2] * 100, 1), "%)")) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # 2. Expression distribution comparison
  ref_mean <- rowMeans(ref_matrix, na.rm = TRUE)
  common_genes <- intersect(names(sample_vec), names(ref_mean))

  dist_data <- data.frame(
    Reference = ref_mean[common_genes],
    Clinical = sample_vec[common_genes]
  ) %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Sample_Type", values_to = "Expression")

  dist_plot <- ggplot(dist_data, aes(x = Expression, fill = Sample_Type)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = c("Clinical" = "red", "Reference" = "lightblue")) +
    labs(title = "Expression Distribution Comparison",
         x = "Log2 Expression", y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # 3. Z-score distribution
  ref_sd <- apply(ref_matrix, 1, sd, na.rm = TRUE)
  z_scores <- (sample_vec - ref_mean) / ref_sd
  z_scores <- z_scores[is.finite(z_scores)]

  zscore_plot <- ggplot(data.frame(Z_Score = z_scores), aes(x = Z_Score)) +
    geom_histogram(bins = 50, fill = "lightcoral", alpha = 0.7) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "darkred") +
    geom_vline(xintercept = c(-3, 3), linetype = "solid", color = "red") +
    labs(title = "Z-Score Distribution (Clinical vs Reference)",
         subtitle = "Dashed lines: ±2SD | Solid lines: ±3SD",
         x = "Z-Score", y = "Count") +
    theme_minimal()

  # 4. Correlation heatmap with reference samples
  correlations <- apply(ref_matrix, 2, function(ref_sample) {
    cor(sample_vec, ref_sample, use = "complete.obs")
  })

  cor_data <- data.frame(
    Sample = 1:length(correlations),
    Correlation = correlations
  )

  cor_plot <- ggplot(cor_data, aes(x = Sample, y = Correlation)) +
    geom_point(alpha = 0.6, color = "darkblue") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    labs(title = "Correlation with Individual Reference Samples",
         subtitle = "Dashed line: 0.8 correlation threshold",
         x = "Reference Sample Index", y = "Pearson Correlation") +
    theme_minimal()

  list(
    pca = pca_plot,
    distribution = dist_plot,
    zscore = zscore_plot,
    correlation = cor_plot
  )
}

#' Generate Batch Effect Recommendations
#'
#' @param metrics List of batch effect metrics
#' @param protocol_clinical Clinical protocol
#' @param protocol_reference Reference protocol
#'
#' @return Character vector of recommendations
#' @keywords internal
generate_batch_recommendations <- function(metrics, protocol_clinical, protocol_reference) {

  recommendations <- character()

  # Protocol-based recommendations
  if (protocol_clinical != protocol_reference) {
    recommendations <- c(recommendations,
      paste("⚠️ Protocol difference detected:", protocol_clinical, "vs", protocol_reference))
  }

  # Z-score based recommendations
  if (metrics$z_score_extreme > 0.1) {
    recommendations <- c(recommendations,
      paste("⚠️ High proportion of extreme Z-scores:",
            round(metrics$z_score_extreme * 100, 1), "% > ±2SD"))
  }

  # Correlation based recommendations
  if (metrics$correlation_median < 0.7) {
    recommendations <- c(recommendations,
      "⚠️ Low correlation with reference cohort suggests batch effects")
  }

  # Expression shift recommendations
  if (abs(metrics$median_expression_shift) > 1) {
    recommendations <- c(recommendations,
      paste("⚠️ Substantial expression level shift detected:",
            round(metrics$median_expression_shift, 2), "log2 units"))
  }

  # Variance ratio recommendations
  if (metrics$variance_ratio > 2 || metrics$variance_ratio < 0.5) {
    recommendations <- c(recommendations,
      paste("⚠️ Expression variance differs substantially from reference:",
            round(metrics$variance_ratio, 2), "x"))
  }

  # General recommendations
  if (length(recommendations) == 0) {
    recommendations <- "✅ No major batch effects detected"
  } else {
    recommendations <- c(recommendations,
      "",
      "🔧 RECOMMENDATIONS:",
      "• Consider using --batch_rm flag",
      "• Validate with internal reference cohort if available",
      "• Use more stringent expression thresholds (>95th percentile)",
      "• Consider protocol-specific normalization")
  }

  recommendations
}

#' Save Batch Assessment Plots
#'
#' @param plots List of ggplot objects
#' @param output_dir Output directory path
#'
#' @keywords internal
save_batch_plots <- function(plots, output_dir) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save each plot
  for (plot_name in names(plots)) {
    filename <- file.path(output_dir, paste0("batch_assessment_", plot_name, ".png"))
    ggplot2::ggsave(filename, plots[[plot_name]], width = 10, height = 8, dpi = 300)
  }

  message("Batch assessment plots saved to: ", output_dir)
}

#' Quick Batch Effect Check
#'
#' Simplified function for rapid batch effect assessment with minimal parameters.
#'
#' @param sample_data Named vector of gene expression values
#' @param reference_data Matrix of reference expression data
#' @param gene_set_type Character, type of gene set: "top_n", "cancer_genes", "custom"
#' @param gene_subset Character vector of custom genes (when gene_set_type = "custom")
#' @param n_genes Integer, number of genes to use (when gene_set_type = "top_n")
#' @param cancer_gene_source Character, cancer gene source when gene_set_type = "cancer_genes"
#' @param save_plots Logical, whether to save diagnostic plots
#'
#' @return Character vector with assessment summary
#' @export
quick_batch_check <- function(sample_data, reference_data, gene_set_type = "top_n",
                             gene_subset = NULL, n_genes = 2000,
                             cancer_gene_source = "umccr", save_plots = FALSE) {

  output_dir <- if (save_plots) "batch_assessment_output" else NULL

  results <- assess_batch_effects(
    sample_data = sample_data,
    reference_data = reference_data,
    gene_set_type = gene_set_type,
    gene_subset = gene_subset,
    n_genes = n_genes,
    cancer_gene_source = cancer_gene_source,
    output_dir = output_dir
  )

  # Generate summary
  summary_text <- c(
    "=== BATCH EFFECT ASSESSMENT SUMMARY ===",
    "",
    paste("Gene set:", results$gene_set_info$description),
    paste("Genes analyzed:", results$genes_analyzed),
    paste("Mean |Z-score|:", round(results$batch_metrics$z_score_mean, 2)),
    paste("Extreme Z-scores:", round(results$batch_metrics$z_score_extreme * 100, 1), "%"),
    paste("Median correlation:", round(results$batch_metrics$correlation_median, 2)),
    paste("PCA distance percentile:", round(results$pca_results$distance_percentile * 100, 1), "%"),
    "",
    "RECOMMENDATIONS:",
    results$recommendations
  )

  cat(paste(summary_text, collapse = "\n"))

  return(invisible(results))
}