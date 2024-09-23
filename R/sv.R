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

  subset_genes <- function(genes, ind) {
    genes |>
      stringr::str_split("&") %>%
      purrr::map(
        ~ .[ind] %>%
          replace("", NA) %>%
          .[!is.na(.)]
      ) %>%
      purrr::map_chr(~ ifelse(length(.) > 0, stringr::str_c(., collapse = "&"), ""))
  }

  # Check file is not empty
  sv_all <- NULL
  if (length(readLines(con = sv_file_sash, n = 2)) <= 1) {
    return(sv_all)
  }

  sv_all <- readr::read_tsv(sv_file_sash, col_names = TRUE)
  total_variants <- nrow(sv_all)
  sv_all <- sv_all |>
    split_sv_field("AF_PURPLE", is_pct = T) |>
    dplyr::mutate(
      start = format(.data$start, big.mark = ",", trim = T),
      location = stringr::str_c(.data$chrom, ":", .data$start, sep = ""),
    ) |>
    # Unpack multiple annotations per region
    dplyr::mutate(annotation = strsplit(.data$annotation, ",")) |>
    tidyr::unnest("annotation") |>
    tidyr::separate_wider_delim(
      cols = "annotation", delim = "|",
      names = c("Event", "Effect", "Genes", "Transcript", "Detail", "Tier"), too_few = "align_start"
    ) |>
    dplyr::mutate(
      Gene = subset_genes(.data$Genes, c(1, 2)),
      Gene = ifelse((stringr::str_split(.data$Genes, "&") |> purrr::map_int(base::length)) > 2,
                    stringr::str_c(.data$Gene, "...", sep = ", "),
                    .data$Gene
      ),
      `Other affected genes` = subset_genes(.data$Genes, -c(1, 2)) |> stringr::str_replace_all("&", ", "),
      Gene = ifelse(stringr::str_detect(.data$Effect, "gene_fusion"),
                    .data$Gene,
                    .data$Gene |> stringr::str_replace_all("&", ", ")
      )
    ) |>
    dplyr::select(
      Genes = "Gene", Location = "location"
    ) |>
    # filter out empty gene rows
    dplyr::filter(.data$Genes != "") |>
    dplyr::distinct() |>
    dplyr::arrange(.data$Genes)
  total_melted <- nrow(sv_all)
  return(list(
    melted = sv_all,
    total_variants = total_variants,
    total_melted = total_melted
  ))
}
