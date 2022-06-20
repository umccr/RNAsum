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
        dplyr::select(-caller, -sample) %>%
        split_sv_field(AF_BPI, is_pct = T) %>%
        split_sv_field(AF_PURPLE, is_pct = T) %>%
        split_sv_field(CN_PURPLE) %>%
        split_sv_field(CN_change_PURPLE) %>%
        dplyr::mutate(
          Ploidy_PURPLE = as.double(Ploidy_PURPLE),
          Ploidy_PURPLE = format(Ploidy_PURPLE, nsmall = 2)
        ) %>%
        tidyr::separate(split_read_support, c("SR (ref)", "SR (alt)"), ",") %>%
        dplyr::mutate(SR = as.integer(`SR (alt)`)) %>%
        tidyr::separate(paired_support_PR, c("PR (ref)", "PR (alt)"), ",") %>%
        dplyr::mutate(PR = as.integer(`PR (alt)`)) %>%
        tidyr::separate(paired_support_PE, c("PE (ref)", "PE (alt)"), ",") %>%
        dplyr::mutate(PE = as.integer(`PE (alt)`)) %>%

        dplyr::filter(svtype != 'BND' | is.na(SR) | PR>SR) %>%  # remove BND with split read support higher than paired
        tidyr::unnest(annotation = strsplit(annotation, ',')) %>%  # Unpack multiple annotations per region
        tidyr::separate(annotation,
                        c('Event', 'Effect', 'Genes', 'Transcript', 'Detail', 'Tier'),
                        sep = '\\|', convert = TRUE) %>%  # Unpack annotation columns
        dplyr::mutate(start = format(start, big.mark = ',', trim = T),
                      end = format(end, big.mark = ',', trim = T)) %>%
        dplyr::mutate(location = stringr::str_c(chrom, ':', start, sep = ''),
                      location = ifelse(is.na(end), location, stringr::str_c(location))) %>%
        dplyr::arrange(Tier, Effect, dplyr::desc(AF_PURPLE), Genes) %>%
        dplyr::mutate(Gene = subset_genes(Genes, c(1, 2)),
                      Gene = ifelse((stringr::str_split(Genes, '&') %>% purrr::map_int(length)) > 2,
                                    stringr::str_c(Gene, '...', sep = ', '),
                                    Gene),
                      `Other affected genes` = subset_genes(Genes, -c(1,2)) %>% stringr::str_replace_all('&', ', '),
                      Gene = ifelse(stringr::str_detect(Effect, "gene_fusion"),
                                    Gene,
                                    Gene %>% stringr::str_replace_all('&', ', '))
        ) %>%
        tidyr::separate(Effect, c("Effect", "Other effects"), sep = '&') %>%
        dplyr::select(Tier = tier, Event = svtype, Gene, Effect = Effect, Detail = Detail, Location = location, AF = AF_PURPLE, `CN chg` = CN_change_PURPLE, SR, PR, CN = CN_PURPLE, Ploidy = Ploidy_PURPLE, PURPLE_status, `SR (ref)`, `PR (ref)`, PE, `PE (ref)`, `Somatic score` = somaticscore, Transcript = Transcript, `Other effects`, `Other affected genes`, `AF at breakpoint 1` = AF_PURPLE1, `AF at breakpoint 2` = AF_PURPLE2, `CN at breakpoint 1` = CN_PURPLE1, `CN at breakpoint 2` = CN_PURPLE2, `CN change at breakpoint 1` = CN_change_PURPLE1, `CN change at breakpoint 2` = CN_change_PURPLE2, `AF before adjustment, bp 1` = AF_BPI1, `AF before adjustment, bp 2` = AF_BPI2
        ) %>%
        dplyr::distinct()
      # dplyr::mutate(chr = factor(chr, levels = c(1:22, "X", "Y", "MT"))) %>%

    } else {
      sv_all <- readr::read_tsv(sv_file, col_names = TRUE) %>%
        dplyr::select(-caller, -sample) %>%
        split_sv_field(BPI_AF, is_pct = T) %>%
        split_sv_field(AF, is_pct = T) %>%
        split_sv_field(CN) %>%
        split_sv_field(CN_change) %>%
        dplyr::mutate(
          Ploidy = as.double(Ploidy),
          Ploidy = format(Ploidy, nsmall = 2)
        ) %>%
        tidyr::separate(split_read_support, c("SR (ref)", "SR (alt)"), ",") %>%
        dplyr::mutate(SR = as.integer(`SR (alt)`)) %>%
        tidyr::separate(paired_support_PR, c("PR (ref)", "PR (alt)"), ",") %>%
        dplyr::mutate(PR = as.integer(`PR (alt)`)) %>%
        tidyr::separate(paired_support_PE, c("PE (ref)", "PE (alt)"), ",") %>%
        dplyr::mutate(PE = as.integer(`PE (alt)`)) %>%

        dplyr::filter(svtype != 'BND' | is.na(SR) | PR>SR) %>%  # remove BND with split read support higher than paired
        tidyr::unnest(annotation = strsplit(annotation, ',')) %>%  # Unpack multiple annotations per region
        tidyr::separate(annotation,
                        c('Event', 'Effect', 'Genes', 'Transcript', 'Detail', 'Tier'),
                        sep = '\\|', convert = TRUE) %>%  # Unpack annotation columns
        dplyr::mutate(start = format(start, big.mark = ',', trim = T),
                      end = format(end, big.mark = ',', trim = T)) %>%
        dplyr::mutate(location = stringr::str_c(chrom, ':', start, sep = ''),
                      location = ifelse(is.na(end), location, stringr::str_c(location))) %>%
        dplyr::arrange(Tier, Effect, dplyr::desc(AF), Genes) %>%
        dplyr::mutate(Gene = subset_genes(Genes, c(1, 2)),
                      Gene = ifelse((stringr::str_split(Genes, '&') %>% purrr::map_int(length)) > 2,
                                    stringr::str_c(Gene, '...', sep = ', '),
                                    Gene),
                      `Other affected genes` = subset_genes(Genes, -c(1,2)) %>% stringr::str_replace_all('&', ', '),
                      Gene = ifelse(stringr::str_detect(Effect, "gene_fusion"),
                                    Gene,
                                    Gene %>% stringr::str_replace_all('&', ', '))
        ) %>%
        tidyr::separate(Effect, c("Effect", "Other effects"), sep = '&') %>%
        dplyr::select(Tier = tier, Event = svtype, Gene, Effect = Effect, Detail = Detail, Location = location, AF, `CN chg` = CN_change, SR, PR, CN, Ploidy, PURPLE_status, `SR (ref)`, `PR (ref)`, PE, `PE (ref)`, `Somatic score` = somaticscore, Transcript = Transcript, `Other effects`, `Other affected genes`, `AF at breakpoint 1` = AF1, `AF at breakpoint 2` = AF2, `CN at breakpoint 1` = CN1, `CN at breakpoint 2` = CN2, `CN change at breakpoint 1` = CN_change1, `CN change at breakpoint 2` = CN_change2, `AF before adjustment, bp 1` = BPI_AF1, `AF before adjustment, bp 2` = BPI_AF2
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
