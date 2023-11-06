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
  melted <- manta_tsv_obj[["melted"]]
  m1 <- melted |>
    dplyr::mutate(multigene = grepl("&", .data$Genes))
  multigenes <- m1 |>
    dplyr::filter(.data$multigene)
  nomultigenes <- m1 |>
    dplyr::filter(!.data$multigene)
  if (nrow(multigenes) > 0) {
    multigenes <- multigenes |>
      dplyr::rowwise() |>
      dplyr::mutate(g = list(unlist(strsplit(.data$Genes, split = "&", fixed = TRUE)[[1]]))) |>
      dplyr::pull("g") |>
      unlist() |>
      unique() |>
      tibble::as_tibble_col(column_name = "Genes") |>
      dplyr::arrange(.data$Genes)
  } else {
    multigenes <- empty_tbl(cnames = "Genes")
  }
  if (nrow(nomultigenes) > 0) {
    nomultigenes <- nomultigenes |>
      dplyr::select("Genes")
  } else {
    nomultigenes <- empty_tbl(cnames = "Genes")
  }

  all_genes <- dplyr::bind_rows(
    multigenes, nomultigenes
  ) |>
    dplyr::distinct(.data$Genes)




  list(
    all_genes = all_genes,
    total_variants = total_variants,
    melted_variants = melted,
    multigenes = multigenes
  )
}

sv_prioritize_old <- function(sv_file) {
  # grab the dplyr pipe
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
  sv_all <- readr::read_tsv(sv_file, col_names = TRUE) |>
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
    tidyr::unnest(annotation = strsplit(annotation, ",")) |> # Unpack multiple annotations per region
    tidyr::separate(annotation,
      c("Event", "Effect", "Genes", "Transcript", "Detail", "Tier"),
      sep = "\\|", convert = TRUE, fill = "right"
    ) %>% # Unpack annotation columns
    dplyr::mutate(
      start = format(start, big.mark = ",", trim = T),
      end = format(end, big.mark = ",", trim = T)
    ) %>%
    dplyr::mutate(
      location = stringr::str_c(chrom, ":", start, sep = ""),
      location = ifelse(is.na(end), location, stringr::str_c(location))
    ) %>%
    dplyr::arrange(Tier, Effect, desc(AF_PURPLE), Genes) %>%
    dplyr::mutate(
      Gene = subset_genes(Genes, c(1, 2)),
      Gene = ifelse((stringr::str_split(Genes, "&") %>% purrr::map_int(length)) > 2,
        stringr::str_c(Gene, "...", sep = ", "),
        Gene
      ),
      `Other affected genes` = subset_genes(Genes, -c(1, 2)) %>% stringr::str_replace_all("&", ", "),
      Gene = ifelse(stringr::str_detect(Effect, "gene_fusion"),
        Gene,
        Gene %>% stringr::str_replace_all("&", ", ")
      )
    ) %>%
    tidyr::separate(Effect, c("Effect", "Other effects"), sep = "&", fill = "right", extra = "merge") %>%
    dplyr::select(Tier = tier, Event = svtype, Genes = Gene, Effect = Effect, Detail = Detail, Location = location, AF = AF_PURPLE, `CN chg` = CN_change_PURPLE, SR, PR, CN = CN_PURPLE, Ploidy = Ploidy_PURPLE, PURPLE_status, `SR (ref)`, `PR (ref)`, PE, `PE (ref)`, `Somatic score` = somaticscore, Transcript = Transcript, `Other effects`, `Other affected genes`, `AF at breakpoint 1` = AF_PURPLE1, `AF at breakpoint 2` = AF_PURPLE2, `CN at breakpoint 1` = CN_PURPLE1, `CN at breakpoint 2` = CN_PURPLE2, `CN change at breakpoint 1` = CN_change_PURPLE1, `CN change at breakpoint 2` = CN_change_PURPLE2, `AF before adjustment, bp 1` = AF_BPI1, `AF before adjustment, bp 2` = AF_BPI2) %>%
    dplyr::distinct() |>
    # filter out empty gene rows
    dplyr::filter(Genes != "")
  total_melted <- nrow(sv_all)
  return(list(
    melted = sv_all,
    total_variants = total_variants,
    total_melted = total_melted
  ))
}
