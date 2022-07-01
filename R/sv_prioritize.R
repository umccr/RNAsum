#' Prioritize SV events
#'
#' Code from umccrise to prioritise SV events (version for "-manta.tsv" files)
#' @param sv_file Input structural variants files.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @return Prioritized variants list.
#' @export
sv_prioritize <- function(sv_file) {

  sv_all = NULL

  if (length(readLines(con = sv_file, n = 2)) > 1) {

    ##### Due to changes in PURPLE output format there are two expected column names combinations
    if ( all(c("AF_BPI", "AF_PURPLE", "CN_PURPLE", "CN_change_PURPLE", "Ploidy_PURPLE") %in% names(readr::read_tsv(sv_file, col_names = TRUE))) ) {

      sv_all <- readr::read_tsv(sv_file, col_names = TRUE) %>%
        dplyr::select(-.data$caller, -.data$sample) %>%
        split_sv_field(.data$AF_BPI, is_pct = T) %>%
        split_sv_field(.data$AF_PURPLE, is_pct = T) %>%
        split_sv_field(.data$CN_PURPLE) %>%
        split_sv_field(.data$CN_change_PURPLE) %>%
        dplyr::mutate(
          Ploidy_PURPLE = as.double(.data$Ploidy_PURPLE),
          Ploidy_PURPLE = format(.data$Ploidy_PURPLE, nsmall = 2)
        ) %>%
        tidyr::separate(.data$split_read_support, c("SR (ref)", "SR (alt)"), ",") %>%
        dplyr::mutate(SR = as.integer(.data$`SR (alt)`)) %>%
        tidyr::separate(.data$paired_support_PR, c("PR (ref)", "PR (alt)"), ",") %>%
        dplyr::mutate(PR = as.integer(.data$`PR (alt)`)) %>%
        tidyr::separate(.data$paired_support_PE, c("PE (ref)", "PE (alt)"), ",") %>%
        dplyr::mutate(PE = as.integer(.data$`PE (alt)`)) %>%
        dplyr::filter(.data$svtype != 'BND' | is.na(.data$SR) | .data$PR > .data$SR) %>%  # remove BND with split read support higher than paired
        tidyr::unnest(annotation = strsplit(.data$annotation, ',')) %>%  # Unpack multiple annotations per region
        tidyr::separate(.data$annotation,
                        c('Event', 'Effect', 'Genes', 'Transcript', 'Detail', 'Tier'),
                        sep = '\\|', convert = TRUE) %>%  # Unpack annotation columns
        dplyr::mutate(start = format(.data$start, big.mark = ',', trim = T),
                      end = format(.data$end, big.mark = ',', trim = T)) %>%
        dplyr::mutate(location = stringr::str_c(.data$chrom, ':', .data$start, sep = ''),
                      location = ifelse(is.na(.data$end), .data$location, stringr::str_c(.data$location))) %>%
        dplyr::arrange(.data$Tier, .data$Effect, dplyr::desc(.data$AF_PURPLE), .data$Genes) %>%
        dplyr::mutate(Gene = subset_genes(.data$Genes, c(1, 2)),
                      Gene = ifelse((stringr::str_split(.data$Genes, '&') %>% purrr::map_int(length)) > 2,
                                    stringr::str_c(.data$Gene, '...', sep = ', '),
                                    .data$Gene),
                      `Other affected genes` = subset_genes(.data$Genes, -c(1,2)) %>% stringr::str_replace_all('&', ', '),
                      Gene = ifelse(stringr::str_detect(.data$Effect, "gene_fusion"),
                                    .data$Gene,
                                    .data$Gene %>% stringr::str_replace_all('&', ', '))
        ) %>%
        tidyr::separate(.data$Effect, c("Effect", "Other effects"), sep = '&') %>%
        dplyr::select(Tier = .data$tier, Event = .data$svtype, .data$Gene, .data$Effect,
                      .data$Detail, Location = .data$location, AF = .data$AF_PURPLE,
                      `CN chg` = .data$CN_change_PURPLE, .data$SR, .data$PR, CN = .data$CN_PURPLE,
                      Ploidy = .data$Ploidy_PURPLE, .data$PURPLE_status, .data$`SR (ref)`, .data$`PR (ref)`,
                      .data$PE, .data$`PE (ref)`, `Somatic score` = .data$somaticscore, .data$Transcript, .data$`Other effects`,
                      .data$`Other affected genes`, `AF at breakpoint 1` = .data$AF_PURPLE1, `AF at breakpoint 2` = .data$AF_PURPLE2,
                      `CN at breakpoint 1` = .data$CN_PURPLE1, `CN at breakpoint 2` = .data$CN_PURPLE2,
                      `CN change at breakpoint 1` = .data$CN_change_PURPLE1, `CN change at breakpoint 2` = .data$CN_change_PURPLE2,
                      `AF before adjustment, bp 1` = .data$AF_BPI1, `AF before adjustment, bp 2` = .data$AF_BPI2
        ) %>%
        dplyr::distinct()
      # dplyr::mutate(chr = factor(chr, levels = c(1:22, "X", "Y", "MT"))) %>%

    } else {
      sv_all <- readr::read_tsv(sv_file, col_names = TRUE) %>%
        dplyr::select(-.data$caller, -.data$sample) %>%
        split_sv_field(.data$BPI_AF, is_pct = T) %>%
        split_sv_field(.data$AF, is_pct = T) %>%
        split_sv_field(.data$CN) %>%
        split_sv_field(.data$CN_change) %>%
        dplyr::mutate(
          Ploidy = as.double(.data$Ploidy),
          Ploidy = format(.data$Ploidy, nsmall = 2)
        ) %>%
        tidyr::separate(.data$split_read_support, c("SR (ref)", "SR (alt)"), ",") %>%
        dplyr::mutate(SR = as.integer(.data$`SR (alt)`)) %>%
        tidyr::separate(.data$paired_support_PR, c("PR (ref)", "PR (alt)"), ",") %>%
        dplyr::mutate(PR = as.integer(.data$`PR (alt)`)) %>%
        tidyr::separate(.data$paired_support_PE, c("PE (ref)", "PE (alt)"), ",") %>%
        dplyr::mutate(PE = as.integer(.data$`PE (alt)`)) %>%

        dplyr::filter(.data$svtype != 'BND' | is.na(.data$SR) | .data$PR > .data$SR) %>%  # remove BND with split read support higher than paired
        tidyr::unnest(annotation = strsplit(.data$annotation, ',')) %>%  # Unpack multiple annotations per region
        tidyr::separate(.data$annotation,
                        c('Event', 'Effect', 'Genes', 'Transcript', 'Detail', 'Tier'),
                        sep = '\\|', convert = TRUE) %>%  # Unpack annotation columns
        dplyr::mutate(start = format(.data$start, big.mark = ',', trim = T),
                      end = format(.data$end, big.mark = ',', trim = T)) %>%
        dplyr::mutate(location = stringr::str_c(.data$chrom, ':', .data$start, sep = ''),
                      location = ifelse(is.na(.data$end), .data$location, stringr::str_c(.data$location))) %>%
        dplyr::arrange(.data$Tier, .data$Effect, dplyr::desc(.data$AF), .data$Genes) %>%
        dplyr::mutate(Gene = subset_genes(.data$Genes, c(1, 2)),
                      Gene = ifelse((stringr::str_split(.data$Genes, '&') %>% purrr::map_int(length)) > 2,
                                    stringr::str_c(.data$Gene, '...', sep = ', '),
                                    .data$Gene),
                      `Other affected genes` = subset_genes(.data$Genes, -c(1,2)) %>% stringr::str_replace_all('&', ', '),
                      Gene = ifelse(stringr::str_detect(.data$Effect, "gene_fusion"),
                                    .data$Gene,
                                    .data$Gene %>% stringr::str_replace_all('&', ', '))
        ) %>%
        tidyr::separate(.data$Effect, c("Effect", "Other effects"), sep = '&') %>%
        dplyr::select(Tier = .data$tier, Event = .data$svtype, .data$Gene, .data$Effect, .data$Detail,
                      Location = .data$location, .data$AF, `CN chg` = .data$CN_change, .data$SR, .data$PR,
                      .data$CN, .data$Ploidy, .data$PURPLE_status, .data$`SR (ref)`, .data$`PR (ref)`, .data$PE,
                      .data$`PE (ref)`, `Somatic score` = .data$somaticscore, .data$Transcript,
                      .data$`Other effects`, .data$`Other affected genes`, `AF at breakpoint 1` = .data$AF1,
                      `AF at breakpoint 2` = .data$AF2, `CN at breakpoint 1` = .data$CN1,
                      `CN at breakpoint 2` = .data$CN2, `CN change at breakpoint 1` = .data$CN_change1,
                      `CN change at breakpoint 2` = .data$CN_change2, `AF before adjustment, bp 1` = .data$BPI_AF1,
                      `AF before adjustment, bp 2` = .data$BPI_AF2
        ) %>%
        dplyr::distinct()
      # dplyr::mutate(chr = factor(chr, levels = c(1:22, "X", "Y", "MT"))) %>%
    }
  } else {
    warning('No prioritized events detected')
  }
  return( sv_all )
}

subset_genes = function(genes, ind) {
  genes %>% stringr::str_split('&') %>% purrr::map(~ .[ind] %>% base::replace("", NA) %>% .[!is.na(.)]) %>% purrr::map_chr(~ ifelse(length(.) > 0, stringr::str_c(., collapse = '&'), ""))
}

format_val = function(val, is_pct = F) {
  ifelse(!is.na(val),
         format(val,  digits = 1) %>% stringr::str_c(ifelse(is_pct, "%", "")), NA)
}

split_sv_field = function(.data, field, is_pct = F) {
  f_q = rlang::enquo(field)
  f_str = rlang::quo_name(f_q)
  f1_str = stringr::str_c(f_str, '1')
  f2_str = stringr::str_c(f_str, '2')
  f1_q = rlang::sym(f1_str)
  f2_q = rlang::sym(f2_str)
  .data %>%
    tidyr::separate(!!f_q, c(f1_str, f2_str), ",") %>%
    dplyr::mutate(
      !!f1_q := as.double(!!f1_q) * ifelse(is_pct, 100, 1),
      !!f2_q := as.double(!!f2_q) * ifelse(is_pct, 100, 1),
      !!f_q  := (!!f1_q + ifelse(is.na(!!f2_q), !!f1_q, !!f2_q)) / 2,
      !!f_q  := format_val(!!f_q, is_pct),
      !!f1_q := format_val(!!f1_q, is_pct),
      !!f2_q := format_val(!!f2_q, is_pct)
    )
}
