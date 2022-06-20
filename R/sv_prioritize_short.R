#' Prioritize structural variants events
#'
#' Code from umccrise to prioritise SV events (version for "-sv-prioritize-manta-pass.tsv" files)
#'
#' @param sv_file Input file.
#'
#' @return Prioritized variant calls.
#' @export
sv_prioritize_short <- function(sv_file) {

  sv_all = NULL

  if (base::length(base::readLines(con = sv_file, n = 2)) > 1) {
    sv_all <- readr::read_tsv(sv_file, col_names = TRUE) %>%
      tidyr::unnest(annotation = strsplit(annotation, ',')) %>% # Unpack multiple annotations per region
      tidyr::separate(annotation,
                      c('Event', 'Annotation', 'Gene', 'Transcript', 'Priority', 'Tier'),
                      sep = '\\|', convert = TRUE) %>% # Unpack annotation columns %>%
      dplyr::mutate(start = format(start, big.mark = ',', trim = T),
                    end = format(end, big.mark = ',', trim = T)) %>%
      dplyr::mutate(Location = stringr::str_c(chrom, ':', start, sep = ''),
                    Location = ifelse(is.na(end), Location, stringr::str_c(Location))) %>%
      dplyr::mutate(SR = split_read_support, PR = paired_support_PR) %>%
      dplyr::select(Location, Gene, Priority, Tier, Annotation, Event, SR, PR) %>%
      dplyr::distinct()
    # dplyr::mutate(Chrom = factor(Chrom, levels = c(1:22, "X", "Y", "MT")))
  } else {
    warning('No prioritized events detected')
  }
  return( sv_all )
}
