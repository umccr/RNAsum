#' Read DRAGEN Mapping Metrics
#'
#' Reads the `mapping_metrics.csv` file output by DRAGEN.
#' @param x Path to the `mapping_metrics.csv` file.
#' @return A tibble with the contents of the input TSV file, or NULL if input is NULL or TSV has no data rows.
#' @examples
#' x <- system.file("rawdata/test_data/dragen",
#'   "test.mapping_metrics.csv",
#'   package = "RNAsum"
#' )
#' (d1 <- dragen_mapping_metrics_read(x))
#' @testexamples
#' expect_equal(colnames(d1)[ncol(d1)], "pct")
#' expect_equal(colnames(d1)[1], "category")
#' expect_null(dragen_mapping_metrics_read(NULL))
#' @export
dragen_mapping_metrics_read <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  abbrev_nm <- c(
    "Total input reads" = "reads_tot_input_dragen",
    "Total reads in RG" = "reads_tot_rg_dragen",
    "Number of duplicate marked reads" = "reads_num_dupmarked_dragen",
    "Number of duplicate marked and mate reads removed" = "reads_num_dupmarked_mate_reads_removed_dragen",
    "Number of unique reads (excl. duplicate marked reads)" = "reads_num_unique_dragen",
    "Reads with mate sequenced" = "reads_w_mate_seq_dragen",
    "Reads without mate sequenced" = "reads_wo_mate_seq_dragen",
    "QC-failed reads" = "reads_qcfail_dragen",
    "Mapped reads" = "reads_mapped_dragen",
    "Mapped reads adjusted for filtered mapping" = "reads_mapped_adjfilt_dragen",
    "Mapped reads R1" = "reads_mapped_r1_dragen",
    "Mapped reads R2" = "reads_mapped_r2_dragen",
    "Number of unique & mapped reads (excl. duplicate marked reads)" = "reads_num_uniq_mapped_dragen",
    "Unmapped reads" = "reads_unmapped_dragen",
    "Unmapped reads adjusted for filtered mapping" = "reads_unmapped_adjfilt_dragen",
    "Adjustment of reads matching non-reference decoys" = "reads_match_nonref_decoys_adj_dragen",
    "Singleton reads (itself mapped; mate unmapped)" = "reads_singleton_dragen",
    "Paired reads (itself & mate mapped)" = "reads_paired_dragen",
    "Properly paired reads" = "reads_paired_proper_dragen",
    "Not properly paired reads (discordant)" = "reads_discordant_dragen",
    "Paired reads mapped to different chromosomes" = "reads_paired_mapped_diff_chrom_dragen",
    "Paired reads mapped to different chromosomes (MAPQ>=10)" = "reads_paired_mapped_diff_chrom_mapq10_dragen",
    "Reads with MAPQ [40:inf)" = "reads_mapq_40_inf_dragen",
    "Reads with MAPQ [30:40)" = "reads_mapq_30_40_dragen",
    "Reads with MAPQ [20:30)" = "reads_mapq_20_30_dragen",
    "Reads with MAPQ [10:20)" = "reads_mapq_10_20_dragen",
    "Reads with MAPQ [ 0:10)" = "reads_mapq_0_10_dragen",
    "Reads with MAPQ NA (Unmapped reads)" = "reads_mapq_na_unmapped_dragen",
    "Reads with indel R1" = "reads_indel_r1_dragen",
    "Reads with indel R2" = "reads_indel_r2_dragen",
    "Total bases" = "bases_tot_dragen",
    "Total bases R1" = "bases_tot_r1_dragen",
    "Total bases R2" = "bases_tot_r2_dragen",
    "Mapped bases R1" = "bases_mapped_r1_dragen",
    "Mapped bases R2" = "bases_mapped_r2_dragen",
    "Soft-clipped bases R1" = "bases_softclip_r1_dragen",
    "Soft-clipped bases R2" = "bases_softclip_r2_dragen",
    "Mismatched bases R1" = "bases_mismatched_r1_dragen",
    "Mismatched bases R2" = "bases_mismatched_r2_dragen",
    "Mismatched bases R1 (excl. indels)" = "bases_mismatched_r1_noindels_dragen",
    "Mismatched bases R2 (excl. indels)" = "bases_mismatched_r2_noindels_dragen",
    "Q30 bases" = "bases_q30_dragen",
    "Q30 bases R1" = "bases_q30_r1_dragen",
    "Q30 bases R2" = "bases_q30_r2_dragen",
    "Q30 bases (excl. dups & clipped bases)" = "bases_q30_nodups_noclipped_dragen",
    "Total alignments" = "alignments_tot_dragen",
    "Secondary alignments" = "alignments_secondary_dragen",
    "Supplementary (chimeric) alignments" = "alignments_chimeric_dragen",
    "Estimated read length" = "read_len_dragen",
    "Bases in reference genome" = "bases_in_ref_genome_dragen",
    "Bases in target bed [% of genome]" = "bases_in_target_bed_genome_pct_dragen",
    "Average sequenced coverage over genome" = "cov_avg_seq_over_genome_dragen",
    "Insert length: mean" = "insert_len_mean_dragen",
    "Insert length: median" = "insert_len_median_dragen",
    "Insert length: standard deviation" = "insert_len_std_dev_dragen",
    "Provided sex chromosome ploidy" = "ploidy_sex_chrom_provided_dragen",
    "Estimated sample contamination" = "contamination_est_dragen",
    "DRAGEN mapping rate [mil. reads/second]" = "mapping_rate_dragen_milreads_per_sec_dragen",
    "Adjustment of reads matching filter contigs" = "reads_match_filter_contigs_adj_dragen",
    "Reads with splice junction" = "reads_splice_junction_dragen"
  )
  d <- readr::read_lines(x)
  if (length(d) == 0) { # this should never happen really
    return(NULL)
  }
  assertthat::assert_that(grepl("MAPPING/ALIGNING", d[1]))

  d |>
    tibble::enframe(name = "name", value = "value") |>
    tidyr::separate_wider_delim(cols = "value", delim = ",", names = c("category", "rg", "extra"), too_many = "merge") |>
    tidyr::separate_wider_delim(cols = "extra", delim = ",", names = c("var", "value"), too_many = "merge") |>
    tidyr::separate_wider_delim(cols = "value", delim = ",", names = c("count", "pct"), too_few = "align_start") |>
    dplyr::mutate(
      category = dplyr::case_when(
        grepl("ALIGNING SUMMARY", .data$category) ~ "summary",
        grepl("ALIGNING PER RG", .data$category) ~ "readgroup",
        TRUE ~ "unknown"
      ),
      rg = ifelse(.data$rg == "", "TOTAL", .data$rg),
      var_abbrev = dplyr::recode(.data$var, !!!abbrev_nm),
      count = sub("NA", NA, .data$count),
      count = as.numeric(.data$count),
      pct = sub("NA", NA, .data$pct),
      pct = as.numeric(.data$pct)
    ) |>
    dplyr::select(
      "category", "rg",
      "var", "var_abbrev", "count", "pct"
    )
}
