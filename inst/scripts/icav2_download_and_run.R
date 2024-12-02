# helpers for downloading GDS data for local testing
require(dracarys, include.only = "s3_list_files_dir")
require(dplyr)
require(rportal, include.only = "orca_query_url")
require(glue, include.only = "glue")
require(here, include.only = "here")
require(tibble, include.only = "tibble")
require(tidyr, include.only = "pivot_longer")

# make sure you have logged into AWS and ICA
c("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_REGION") |>
  rportal::envvar_defined() |>
  stopifnot()
token_orca <- rportal::orca_jwt() |> rportal::jwt_validate()
portalRunId <- "20241201d2c4e79d" # failed run

# query orca workflow manager for payload with rnasum s3 inputs
get_rnasum_inputs <- function(prid, token_orca) {
  wfrid <- rportal::orca_prid2wfrid(prid = portalRunId, token = token)
  states <- rportal::orca_wfrid2state(wfrid, token) # view states for gived wfrid
  pldid <- states |>
    filter(.data$status == "FAILED") |>
    pull("payload")
  url <- glue("https://workflow.prod.umccr.org/api/v1/payload/{pldid}")
  pld <- rportal::orca_query_url(url, token)
  s3inputs <- pld[["data"]][["inputs"]] |>
    purrr::map(\(x) x |> stringr::str_replace("/$", ""))
  return(s3inputs)
}

rnasum_inputs <- get_rnasum_inputs(portalRunId, token)
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
outdir <- here::here("nogit/patient_data")

# melt gds_indir columns to get a single column with paths to gds directories
# of interest, and fish out files of interest from each of them, then download
meta_rnasum <- rnasum_inputs |>
  tibble::as_tibble_row() |>
  tidyr::pivot_longer(contains("Uri"), names_to = "folder_type", values_to = "s3_indir") |>
  mutate(folder_type = sub("Uri$", "", .data$folder_type)) |>
  select(subjectId, wtsTumorLibraryId, folder_type, s3_indir) |>
  rowwise() |>
  mutate(
    down = list(dracarys::dr_s3_download(s3dir = s3_indir, outdir = outdir, max_objects = 1000, regexes = rnasum_file_regex, dryrun = F))
  ) |>
  ungroup()

# saveRDS(meta_rnasum, here(glue("nogit/patient_data/down_{date1}.rds")))
# meta_rnasum <- here::here(glue("nogit/patient_data/down_{date1}.rds")) |> readr::read_rds()
rnasum_params_set <- function(arriba_pdf, arriba_tsv, dataset, dragen_fusions, dragen_mapping_metrics, sv_tsv,
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
    sv_tsv = sv_tsv,
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
  dplyr::select(subjectId, libraryId = "wtsTumorLibraryId", type, localpath) |>
  tidyr::pivot_wider(names_from = type, values_from = localpath)

# slice to whichever run you want from d
d_runs |>
  # dplyr::slice(1) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    rnasum_dataset = "PANCAN",
    params = list(
      rnasum_params_set(
        arriba_pdf = ArribaFusionsPdfFile,
        arriba_tsv = ArribaFusionsTsvFile,
        dataset = rnasum_dataset,
        dragen_fusions = DragenWtsFusionsFinalFile,
        dragen_mapping_metrics = MapMetricsFile,
        sv_tsv = MantaTsvFile,
        pcgr_tiers_tsv = PcgrTiersTsvFile,
        purple_gene_tsv = PurpleCnvGeneTsvFile,
        report_dir = here::here(glue::glue("nogit/patient_data/reports/{subjectId}_{libraryId}_{rnasum_dataset}")),
        salmon = DragenWtsQuantGenesSfFile,
        sample_name = glue::glue("{subjectId}_{libraryId}"),
        subject_id = subjectId
      )
    ),
    rnasum_rmd = RNAsum::rnasum_rmd(
      out_file = here::here(glue::glue("nogit/patient_data/reports_html/{subjectId}_{libraryId}_{rnasum_dataset}.html")),
      quiet = FALSE, pars = params
    )
  ) |>
  dplyr::ungroup()
