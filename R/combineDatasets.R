#' Combine sample expression profile with reference datasets
#'
#' Combines sample expression profile with reference datasets. Outputs a vector
#' with first element containing the merged data and second element containing
#' merged targets info
#'
#' @param sample_name Sample name.
#' @param sample_counts Sample counts.
#' @param ref_data Reference data.
#' @param report_dir Report directory.
#'
#' @return Combined sample expression profile with reference datasets.
#' @export
combineDatasets <- function(sample_name, sample_counts, ref_data, report_dir) {
  assertthat::assert_that(!is.null(sample_counts), msg = "The 'sample_counts' cannot be NULL!")
  .read_target <- function(targ, ds) {
    utils::read.table(targ, sep = "\t", as.is = TRUE, header = TRUE) |>
      tibble::as_tibble() |>
      dplyr::mutate(Dataset = ds)
  }
  .read_counts <- function(cnt, ...) {
    utils::read.table(gzfile(cnt), header = TRUE, sep = "\t", row.names = NULL) |>
      tibble::as_tibble()
  }
  target_ext <- .read_target(targ = ref_data[["ext_ref"]][["target"]], ds = ref_data[["ext_ref"]][["dataset"]])
  target_int <- .read_target(targ = ref_data[["int_ref"]][["target"]], ds = ref_data[["int_ref"]][["dataset"]])
  target_all <- dplyr::bind_rows(target_ext, target_int) |>
    dplyr::mutate(Sample_name = glue::glue("{.data$Dataset}.{.data$Sample_name}")) |>
    dplyr::bind_rows(
      tibble::tibble(Sample_name = sample_name, Target = sample_name, Dataset = sample_name)
    ) |>
    dplyr::mutate(Sample_name = make.names(.data$Sample_name))

  datasets.comb <- sample_counts
  names(datasets.comb) <- c("", sample_name)
  count_ext <- .read_counts(ref_data[["ext_ref"]][["counts"]])
  colnames(count_ext) <- glue::glue("{ref_data[['ext_ref']][['dataset']]}.{colnames(count_ext)}")
  count_int <- .read_counts(ref_data[["int_ref"]][["counts"]])
  colnames(count_int) <- glue::glue("{ref_data[['int_ref']][['dataset']]}.{colnames(count_int)}")
  datasets.comb <- base::merge(datasets.comb, count_ext, by = 1, all = FALSE, sort = TRUE)
  datasets.comb <- base::merge(datasets.comb, count_int, by = 1, all = FALSE, sort = TRUE)
  gene_list <- base::unique(c(sample_counts[["rowname"]], count_ext[[1]], count_int[[1]]))

  # TODO: refactor this a bit better
  ##### Use gene IDs as rownames
  rownames(datasets.comb) <- datasets.comb[[1]]
  datasets.comb <- datasets.comb[, -1]
  colnames(datasets.comb) <- make.names(colnames(datasets.comb))

  ##### Make sure that the target file contains info only about samples present in the data matrix
  target_all <- target_all[target_all[["Sample_name"]] %in% colnames(datasets.comb), ]

  ##### Make sure that the samples order in the data matrix is the same as in the target file
  datasets.comb <- datasets.comb[, target_all[["Sample_name"]]]

  ##### Identify genes that were not present across all per-sample files and were ommited in the merged matrix
  gene_list.missing <- gene_list[gene_list %!in% rownames(datasets.comb)]

  ##### Write list of missing genes into a file
  # TODO: not sure why we write this - does it get used?
  if (length(gene_list.missing) > 0) {
    utils::write.table(
      prepare2write(gene_list.missing),
      file = paste0(report_dir, "/", sample_name, ".RNAseq_report.missing_genes.txt"),
      sep = "\t", quote = FALSE, row.names = TRUE, append = FALSE
    )
  }

  # TODO: doing this so that output is identical with old method
  target_all <- as.data.frame(target_all)
  rownames(target_all) <- target_all[["Sample_name"]]
  target_all <- target_all |>
    dplyr::select(-.data$Sample_name)
  return(list(datasets.comb, target_all))
}
