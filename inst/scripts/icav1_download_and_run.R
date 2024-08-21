# helpers for downloading GDS data for local testing
require(dracarys)
require(dplyr)
require(readr)
require(rportal, include.only = "portaldb_query_workflow")
require(glue, include.only = "glue")
require(here, include.only = "here")

# grab rnasum workflow metadata from Athena
athena_rnasum <- function(sbj) {
  q_quote <- shQuote(paste(glue("rnasum__{sbj}"), collapse = "|"))
  query1 <- glue('WHERE REGEXP_LIKE("wfr_name", {q_quote});')
  rportal::portaldb_query_workflow(query1)
}

athena_lims <- function(libid) {
  query1 <- glue("WHERE REGEXP_LIKE(\"library_id\", '{libid}');")
  rportal::portaldb_query_limsrow(query1)
}

# download gds files to a local structure reflecting the gds path starting from
# the outdir as the fs root.
rnasum_download <- function(gdsdir, outdir, token, page_size = 200, regexes) {
  dracarys::gds_files_list(gdsdir = gdsdir, token = token, page_size = page_size) |>
    dplyr::mutate(type = purrr::map_chr(.data$bname, \(x) dracarys::match_regex(x, regexes))) |>
    dplyr::select("file_id", "type", "size", "path", "bname") |>
    dplyr::filter(!is.na(.data$type)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      dname = dirname(.data$path),
      dname = sub("gds://", "", .data$dname),
      local_outdir = fs::dir_create(file.path(outdir, .data$dname)),
      outfile = file.path(local_outdir, .data$bname),
      out_dl = dracarys::gds_file_download_api(.data$file_id, .data$outfile, token)
    ) |>
    dplyr::ungroup()
}

# SBJ IDs of interest
sbj <- "SBJ04426"
lib <- "L2301428"
date1 <- "2024-08-20"
lims_raw <- athena_lims(lib)
pmeta_raw <- athena_rnasum(sbj) |>
  rportal::meta_rnasum()
lims <- lims_raw |>
  dplyr::select(library_id, sample_id, subject_id)

# generate tidy rnasum metadata from portal workflows table, and join against glims
pmeta <- pmeta_raw |>
  dplyr::left_join(lims, by = c("LibraryID" = "library_id", "SampleID" = "sample_id", "SubjectID" = "subject_id")) |>
  dplyr::select(
    gds_indir_dragen, gds_indir_umccrise, gds_indir_arriba,
    SubjectID, LibraryID, SampleID,
    rnasum_dataset,
    end_status,
    wfr_id, start, end, gds_outfile_rnasum_html,
  ) |>
  dplyr::arrange(desc(SubjectID), start) |>
  dplyr::filter(rnasum_dataset == "BRCA") |>
  dplyr::slice_head(n = 1)

# patterns of files to fish out
rnasum_file_regex <- tibble::tribble(
  ~regex, ~fun,
  "quant\\.genes\\.sf$", "DragenWtsQuantGenesSfFile",
  "quant\\.sf$", "DragenWtsQuantSfFile",
  "fusion_candidates\\.final$", "DragenWtsFusionsFinalFile",
  "fusions\\.pdf$", "ArribaFusionsPdfFile",
  "fusions\\.tsv$", "ArribaFusionsTsvFile",
  "somatic\\.pcgr\\.snvs_indels\\.tiers\\.tsv$", "PcgrTiersTsvFile",
  "purple\\.cnv\\.gene\\.tsv$", "PurpleCnvGeneTsvFile",
  "manta\\.tsv$", "MantaTsvFile",
  "mapping_metrics\\.csv$", "MapMetricsFile"
)
token <- dracarys::ica_token_validate()
outdir <- here::here("nogit/patient_data")

# melt gds_indir columns to get a single column with paths to gds directories
# of interest, and fish out files of interest from each of them, then download
meta_rnasum <- pmeta |>
  tidyr::pivot_longer(dplyr::contains("gds_indir"), names_to = "folder_type", values_to = "gds_indir") |>
  dplyr::select(SubjectID, LibraryID, rnasum_dataset, folder_type, gds_indir) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    down = list(rnasum_download(gdsdir = gds_indir, outdir = outdir, token = token, regexes = rnasum_file_regex))
  ) |>
  dplyr::ungroup()

# saveRDS(meta_rnasum, here(glue("nogit/patient_data/down_{date1}.rds")))
# meta_rnasum <- here::here(glue("nogit/patient_data/down_{date1}.rds")) |> readr::read_rds()
rnasum_params_set <- function(arriba_pdf, arriba_tsv, dataset, dragen_fusions, dragen_mapping_metrics, manta_tsv,
                              pcgr_tiers_tsv, purple_gene_tsv, report_dir, salmon,
                              sample_name, subject_id) {
  params <- list(
    arriba_pdf = arriba_pdf,
    arriba_tsv = arriba_tsv,
    batch_rm = TRUE,
    cn_gain = 95,
    cn_loss = 5,
    dataset = dataset,
    dataset_name_incl = FALSE,
    dragen_fusions = dragen_fusions,
    dragen_mapping_metrics = dragen_mapping_metrics,
    drugs = FALSE,
    filter = TRUE,
    immunogram = FALSE,
    log = TRUE,
    manta_tsv = manta_tsv,
    norm = "TMM",
    pcgr_splice_vars = TRUE,
    pcgr_tier = 4,
    pcgr_tiers_tsv = pcgr_tiers_tsv,
    project = "-",
    purple_gene_tsv = purple_gene_tsv,
    report_dir = report_dir,
    salmon = salmon,
    sample_name = sample_name,
    sample_source = "-",
    save_tables = TRUE,
    scaling = "gene-wise",
    subject_id = subject_id,
    top_genes = 5,
    transform = "CPM",
    umccrise = NULL
  )
  params
}

# now unmelt to have one row per run
d_runs <- meta_rnasum |>
  tidyr::unnest(down) |>
  dplyr::select(SubjectID, LibraryID, rnasum_dataset, type, outfile) |>
  tidyr::pivot_wider(names_from = type, values_from = outfile)

# slice to whichever run you want from d
d_runs |>
  dplyr::slice(1) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    params = list(
      rnasum_params_set(
        arriba_pdf = ArribaFusionsPdfFile,
        arriba_tsv = ArribaFusionsTsvFile,
        dataset = rnasum_dataset,
        dragen_fusions = DragenWtsFusionsFinalFile,
        dragen_mapping_metrics = MapMetricsFile,
        manta_tsv = MantaTsvFile,
        pcgr_tiers_tsv = PcgrTiersTsvFile,
        purple_gene_tsv = PurpleCnvGeneTsvFile,
        report_dir = here::here(glue::glue("nogit/patient_data/reports/{SubjectID}_{LibraryID}_{rnasum_dataset}")),
        salmon = DragenWtsQuantGenesSfFile,
        sample_name = glue::glue("{SubjectID}_{LibraryID}"),
        subject_id = SubjectID
      )
    ),
    rnasum_rmd = RNAsum::rnasum_rmd(
      out_file = here::here(glue::glue("nogit/patient_data/reports_html/{SubjectID}_{LibraryID}_{rnasum_dataset}.html")),
      quiet = FALSE, pars = params
    )
  ) |>
  dplyr::ungroup()
