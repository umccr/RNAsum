#' Process SV Object
#'
#' @param sv_tsv_obj SV list object read via `sv_prioritize`.
#'
#' @return List with:
#' - total SV (unmelted) variants
#' - tibble with melted variants
#' - genes involved in multi-gene events
#' @export
sv_process <- function(sv_tsv_obj) {
  total_variants <- sv_tsv_obj[["total_variants"]]
  melted <- sv_tsv_obj[["melted"]] |>
    dplyr::mutate(
      is_fusion = grepl("&", .data$Genes),
      fusion_genes = dplyr::if_else(.data$is_fusion, .data$Genes, "")
    )
  fus <- melted |>
    dplyr::filter(.data$is_fusion)
  nofus <- melted |>
    dplyr::filter(!.data$is_fusion)
  if (nrow(fus) > 0) {
    # split fusions into two rows, one per gene
    fus <- fus |>
      tidyr::separate_longer_delim(cols = "Genes", delim = "&")
  } else {
    fus <- NULL # bind_rows(NULL, x) returns x
  }

  melted <- dplyr::bind_rows(fus, nofus)
  genes_all <- melted |>
    dplyr::select("Genes") |>
    dplyr::distinct(.data$Genes)

  list(
    genes_all = genes_all,
    total_variants = total_variants,
    melted_variants = melted
  )
}

sv_prioritize <- function(sv_file){

  # Check file is not empty
  sv_all <- NULL
  if (length(readLines(con = sv_file, n = 2)) <= 1) {
    return(sv_all)
  }

  sv_all <- readr::read_tsv(sv_file, col_names = TRUE)
  total_variants <- nrow(sv_all)

  # Check if user has provided a tsv with Gene as first column
  if(colnames(sv_all)[1] == "Gene"){
    return(list(
      melted = sv_all,
      total_variants = total_variants
    ))
  }
  if(!"Gene" %in% colnames(sv_all)){
    # Assume it's an internal input. Unpack multiple annotations per region
    sv_all <- sv_all |>
      dplyr::select("annotation") |>
      dplyr::mutate(annotation = strsplit(.data$annotation, ",")) |>
      tidyr::unnest("annotation") |>
      tidyr::separate_wider_delim(
        cols = "annotation", delim = "|",
        names = c("Event", "Effect", "Genes", "Transcript", "Detail", "Tier"), too_few = "align_start"
      ) |>
      dplyr::select("Effect", "Genes") |>
      dplyr::mutate(
        # if gene_fusion, keep as-is
        gene_fusion_effect = grepl("gene_fusion", .data$Effect),
        Gene = ifelse(
          .data$gene_fusion_effect,
          .data$Genes,
          strsplit(.data$Genes, "&")
        )
      ) |>
      dplyr::select("Gene") |>
      tidyr::unnest_longer("Gene") |>
      dplyr::distinct() |>
      dplyr::arrange(.data$Gene)

    return(list(
      melted = sv_all,
      total_variants = total_variants
    ))
  }
}

