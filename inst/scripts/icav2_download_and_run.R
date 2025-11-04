# helpers for downloading S3 data for local testing
{
  use("dracarys", "s3_list_files_dir")
  use("dplyr")
  use("rportal", "orca_query_url")
  use("glue", "glue")
  use("here", "here")
  use("tibble", "tibble")
  use("tidyr", "pivot_longer")
}
# make sure you have logged into AWS
c("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_REGION") |>
  rportal::envvar_defined() |>
  stopifnot()
token_orca <- rportal::orca_jwt() |> rportal::jwt_validate()
portalRunId <- "2025102542c47e55" # failed run

# query orca workflow manager for payload with rnasum s3 inputs
get_rnasum_inputs <- function(prid, token) {
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

rnasum_inputs <- get_rnasum_inputs(portalRunId, token_orca)
# patterns of files to fish out
rnasum_file_regex <- tibble::tribble(
  ~regex                          , ~fun                        ,
  "quant\\.genes\\.sf$"           , "DragenWtsQuantGenesSfFile" ,
  "quant\\.sf$"                   , "DragenWtsQuantSfFile"      ,
  "fusion_candidates\\.final$"    , "DragenWtsFusionsFinalFile" ,
  "fusions\\.pdf$"                , "ArribaFusionsPdfFile"      ,
  "L[A-Z0-9]+_fusions\\.tsv$"     , "ArribaFusionsTsvFile"      , # need to account for discarded ones
  "\\.snvs_indels\\.tiers\\.tsv$" , "PcgrTiersTsvFile"          ,
  "purple\\.cnv\\.gene\\.tsv$"    , "PurpleCnvGeneTsvFile"      ,
  "manta\\.tsv$"                  , "MantaTsvFile"              ,
  "mapping_metrics\\.csv$"        , "MapMetricsFile"            ,
  "sv\\.prioritised\\.tsv$"       , "SvPrioritisedTsvFile"
)
outdir <- here::here("nogit/patient_data")

# melt s3_indir columns to get a single column with paths to S3 directories
# of interest, and fish out files of interest from each of them, then download
# 2025-11-03: now getting file paths explicitly instead of s3 dirs, so need to
# use dirnames
s3dir <- rnasum_inputs |>
  tibble::as_tibble_row() |>
  tidyr::pivot_longer(dplyr::everything()) |>
  filter(grepl("s3://", value)) |>
  mutate(s3_indir = dirname(value)) |>
  distinct(s3_indir) |>
  rowwise() |>
  mutate(
    down = list(dracarys::dr_s3_download(
      s3dir = s3_indir,
      outdir = file.path(outdir, rnasum_inputs$reportDir),
      max_objects = 1000,
      regexes = rnasum_file_regex,
      dryrun = F
    ))
  ) |>
  ungroup()

date1 <- Sys.Date()
# saveRDS(meta_rnasum, here(glue("nogit/patient_data/down_{date1}.rds")))
# meta_rnasum <- here::here(glue("nogit/patient_data/down_{date1}.rds")) |> readr::read_rds()
rnasum_params_set <- function(
  arriba_pdf,
  arriba_tsv,
  dataset,
  dragen_fusions,
  dragen_mapping_metrics,
  sv_tsv,
  pcgr_tiers_tsv,
  purple_gene_tsv,
  report_dir,
  salmon,
  sample_name,
  subject_id
) {
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

local_paths <- s3dir |>
  tidyr::unnest(down) |>
  dplyr::select(type, localpath) |>
  tibble::deframe() |>
  as.list()

params <- rnasum_params_set(
  arriba_pdf = local_paths$ArribaFusionsPdfFile,
  arriba_tsv = local_paths$ArribaFusionsTsvFile,
  dataset = rnasum_inputs$dataset,
  dragen_fusions = local_paths$DragenWtsFusionsFinalFile,
  dragen_mapping_metrics = local_paths$MapMetricsFile,
  sv_tsv = local_paths$SvPrioritisedTsvFile,
  pcgr_tiers_tsv = local_paths$PcgrTiersTsvFile,
  purple_gene_tsv = local_paths$PurpleCnvGeneTsvFile,
  report_dir = here::here(glue::glue(
    "nogit/patient_data/reports/{rnasum_inputs$reportDir}"
  )),
  salmon = local_paths$DragenWtsQuantGenesSfFile,
  sample_name = rnasum_inputs$sampleName,
  subject_id = rnasum_inputs$subjectId
)

RNAsum::rnasum_rmd(
  out_file = here::here(glue::glue(
    "nogit/patient_data/reports_html/{rnasum_inputs$reportDir}.html"
  )),
  quiet = FALSE,
  pars = params
)
