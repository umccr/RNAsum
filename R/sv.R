#' Process Manta Object
#'
#' @param manta_tsv_obj Manta list object read via `sv_prioritize_old`.
#'
#' @return List with:
#' - total Manta (unmelted) variants
#' - tibble with melted variants
#' - genes involved in multi-gene events
#' @export
manta_process <- function(manta_tsv_obj) {
  total_variants <- manta_tsv_obj[["total_variants"]]
  melted <- manta_tsv_obj[["melted"]] |>
    dplyr::mutate(
      is_fusion = grepl("&", .data$Genes),
      fusion_genes = dplyr::if_else(is_fusion, .data$Genes, "")
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
    fus <- empty_tbl(cnames = colnames(melted))
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

sv_prioritize_old <- function(sv_file) {
  # grab the dplyr pipe (for now!)
  `%>%` <- dplyr::`%>%`
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

  format_val <- function(val, is_pct = F) {
    ifelse(!is.na(val),
      format(val, digits = 1) %>% stringr::str_c(ifelse(is_pct, "%", "")),
      NA
    )
  }

  split_sv_field <- function(.data, field, is_pct = F) {
    f_q <- rlang::enquo(field)
    f_str <- rlang::quo_name(f_q)
    f1_str <- stringr::str_c(f_str, "1")
    f2_str <- stringr::str_c(f_str, "2")
    f1_q <- rlang::sym(f1_str)
    f2_q <- rlang::sym(f2_str)
    .data %>%
      tidyr::separate(col = !!f_q, into = c(f1_str, f2_str), sep = ",", fill = "right") %>%
      dplyr::mutate(
        !!f1_q := as.double(!!f1_q) * ifelse(is_pct, 100, 1),
        !!f2_q := as.double(!!f2_q) * ifelse(is_pct, 100, 1),
        !!f_q := (!!f1_q + ifelse(is.na(!!f2_q), !!f1_q, !!f2_q)) / 2,
        !!f_q := format_val(!!f_q, is_pct),
        !!f1_q := format_val(!!f1_q, is_pct),
        !!f2_q := format_val(!!f2_q, is_pct)
      )
  }

  sv_all <- NULL
  if (length(readLines(con = sv_file, n = 2)) <= 1) {
    return(sv_all)
  }

  col_types_tab <- dplyr::tribble(
    ~Column, ~Description, ~Type,
    "caller", "Manta SV caller", "c",
    "sample", "Tumor sample name", "c",
    "chrom", "CHROM column in VCF", "c",
    "start", "POS column in VCF", "i",
    "end", "INFO/END: End position of the variant described in this record", "i",
    "svtype", "INFO/SVTYPE: Type of structural variant", "c",
    "split_read_support", "FORMAT/SR of tumor sample: Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999", "c",
    "paired_support_PE", "FORMAT/PE of tumor sample: ??", "c",
    "paired_support_PR", "FORMAT/PR of tumor sample: Spanning paired-read support for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999", "c",
    "AF_BPI", "INFO/BPI_AF: AF at each breakpoint (so AF_BPI1,AF_BPI2)", "c",
    "somaticscore", "INFO/SOMATICSCORE: Somatic variant quality score", "i",
    "tier", "INFO/SV_TOP_TIER (or 4 if missing): Highest priority tier for the effects of a variant entry", "c",
    "annotation", "INFO/SIMPLE_ANN: Simplified structural variant annotation: 'SVTYPE | EFFECT | GENE(s) | TRANSCRIPT | PRIORITY (1-4)'", "c",
    "AF_PURPLE", "INFO/PURPLE_AF: AF at each breakend (purity adjusted) (so AF_PURPLE1,AF_PURPLE2)", "c",
    "CN_PURPLE", "INFO/PURPLE_CN: CN at each breakend (purity adjusted) (so CN_PURPLE1,CN_PURPLE2)", "c",
    "CN_change_PURPLE", "INFO/PURPLE_CN_CHANGE: change in CN at each breakend (purity adjusted) (so CN_change_PURPLE1,CN_change_PURPLE2)", "c",
    "Ploidy_PURPLE", "INFO/PURPLE_PLOIDY: Ploidy of variant (purity adjusted)", "d",
    "PURPLE_status", "INFERRED if FILTER=INFERRED, or RECOVERED if has INFO/RECOVERED, else blank. INFERRED: Breakend inferred from copy number transition", "c",
    "START_BPI", "INFO/BPI_START: BPI adjusted breakend location", "i",
    "END_BPI", "INFO/BPI_END: BPI adjusted breakend location", "i",
    "ID", "ID column in VCF", "c",
    "MATEID", "INFO/MATEID: ID of mate breakend", "c",
    "ALT", "ALT column in VCF", "c"
  )
  ctypes <- paste(col_types_tab$Type, collapse = "")
  sv_all <- readr::read_tsv(sv_file, col_names = TRUE, col_types = ctypes) |>
    dplyr::select(-c("caller", "sample")) |>
    split_sv_field(AF_BPI, is_pct = T) |>
    split_sv_field(AF_PURPLE, is_pct = T) |>
    split_sv_field(CN_PURPLE) |>
    split_sv_field(CN_change_PURPLE) |>
    dplyr::mutate(
      Ploidy_PURPLE = as.double(Ploidy_PURPLE),
      Ploidy_PURPLE = format(Ploidy_PURPLE, nsmall = 2)
    ) |>
    tidyr::separate(split_read_support, c("SR (ref)", "SR (alt)"), ",", fill = "right") |>
    dplyr::mutate(SR = as.integer(`SR (alt)`)) |>
    tidyr::separate(paired_support_PR, c("PR (ref)", "PR (alt)"), ",", fill = "right") |>
    dplyr::mutate(PR = as.integer(`PR (alt)`)) |>
    tidyr::separate(paired_support_PE, c("PE (ref)", "PE (alt)"), ",", fill = "right") |>
    dplyr::mutate(PE = as.integer(`PE (alt)`)) |>
    dplyr::filter(svtype != "BND" | is.na(SR) | PR > SR) # remove BND with split read support higher than paired
  total_variants <- nrow(sv_all)
  sv_all <- sv_all |>
    # Unpack multiple annotations per region
    dplyr::mutate(annotation = strsplit(.data$annotation, ",")) |>
    tidyr::unnest("annotation") |>
    tidyr::separate_wider_delim(
      cols = "annotation", delim = "|",
      names = c("Event", "Effect", "Genes", "Transcript", "Detail", "Tier"), too_few = "align_start"
    ) |>
    dplyr::mutate(
      start = format(start, big.mark = ",", trim = T),
      end = format(end, big.mark = ",", trim = T),
      location = stringr::str_c(chrom, ":", start, sep = ""),
      location = ifelse(is.na(end), location, stringr::str_c(location))
    ) |>
    dplyr::mutate(
      Gene = subset_genes(Genes, c(1, 2)),
      Gene = ifelse((stringr::str_split(Genes, "&") |> purrr::map_int(length)) > 2,
        stringr::str_c(Gene, "...", sep = ", "),
        Gene
      ),
      `Other affected genes` = subset_genes(Genes, -c(1, 2)) |> stringr::str_replace_all("&", ", "),
      Gene = ifelse(stringr::str_detect(Effect, "gene_fusion"),
        Gene,
        Gene |> stringr::str_replace_all("&", ", ")
      )
    ) |>
    tidyr::separate_wider_delim(
      cols = "Effect", delim = "&", names = c("Effect", "Other effects"),
      too_few = "align_start", too_many = "merge"
    ) |>
    dplyr::select(
      Tier = "tier", Event = "svtype", Genes = "Gene", Effect = "Effect",
      Detail = "Detail", Location = "location", AF = "AF_PURPLE", `CN chg` = "CN_change_PURPLE",
      "SR", "PR", CN = "CN_PURPLE", Ploidy = "Ploidy_PURPLE", "PURPLE_status",
      "SR (ref)", "PR (ref)", "PE", "PE (ref)", `Somatic score` = "somaticscore",
      "Transcript", "Other effects", "Other affected genes",
      `AF at breakpoint 1` = "AF_PURPLE1", `AF at breakpoint 2` = "AF_PURPLE2",
      `CN at breakpoint 1` = "CN_PURPLE1", `CN at breakpoint 2` = "CN_PURPLE2",
      `CN change at breakpoint 1` = "CN_change_PURPLE1",
      `CN change at breakpoint 2` = "CN_change_PURPLE2",
      `AF before adjustment, bp 1` = "AF_BPI1",
      `AF before adjustment, bp 2` = "AF_BPI2"
    ) |>
    # filter out empty gene rows
    dplyr::filter(Genes != "") |>
    dplyr::distinct() |>
    dplyr::arrange(.data$Tier, .data$Effect, dplyr::desc(.data$AF), .data$Genes)
  total_melted <- nrow(sv_all)
  return(list(
    melted = sv_all,
    total_variants = total_variants,
    total_melted = total_melted
  ))
}
