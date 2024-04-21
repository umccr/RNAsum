#--- Launch rnasum on a set of sbjs ---#
# 1. Grab metadata from prod portaldb for rnasum workflows that have succeeded
#   - this will give us paths to input gds dirs (umccrise, dragen, arriba)
# 2. Fish out the specific files required for rnasum from those input dirs, and
#    generate their presigned urls
# 3. Use the wfl.87e07ae6b46645a181e04813de535216 dev workflow to launch version
#    0.5.0--304c3ac of the rnasum workflow.
#   - this will take in a json object with the required inputs and
#     engineParameters that override the docker container used.

{
  require(dplyr)
  require(glue, include.only = "glue")
  require(here, include.only = "here")
  # devtools::install_github("umccr/rportal")
  # devtools::install_github("umccr/dracarys")
  require(rportal, include.only = "awsvault_profile")
  require(dracarys, include.only = "ica_token_validate")
  require(tibble, include.only = "tribble")
  require(tidyr, include.only = "pivot_wider")
  require(stringr, include.only = "str_replace")
}

# keep long functions in separate file
source(here("inst/scripts/cwl_launch/funcs.R"))
# connect to aws
rportal::awsvault_profile("upro")
# ica token (prod for presigned urls, dev for running)
token_prod <- dracarys::ica_token_validate(Sys.getenv("ICA_ACCESS_TOKEN_PRO"))
token_dev <- dracarys::ica_token_validate(Sys.getenv("ICA_ACCESS_TOKEN_DEV"))
# subjects from prod you want inputs from
sbjids <- c(
  "SBJ02999",
  "SBJ03038",
  "SBJ03039",
  "SBJ04660",
  "SBJ04661",
  "SBJ04662",
  "SBJ04774",
  "SBJ04775",
  "SBJ04776",
  "SBJ04777"
)
# grab rnasum workflow metadata for given sbjids
pmeta_raw <- query_workflow_rnasum(sbjids)
# grab top run of each sbj
pmeta_tidy <- pmeta_raw |>
  rportal::meta_rnasum() |>
  filter(end_status == "Succeeded", rnasum_dataset == "PANCAN") |>
  slice_head(n = 1, by = "SubjectID")

# which files to fish for in the gds dirs
rnasum_file_regex <- tibble::tribble(
  ~regex, ~fun,
  "quant\\.genes\\.sf$", "DragenWtsQuantGenesSfFile",
  "fusion_candidates\\.final$", "DragenWtsFusionsFinalFile",
  "fusions\\.pdf$", "ArribaFusionsPdfFile",
  "fusions\\.tsv$", "ArribaFusionsTsvFile",
  "somatic\\.pcgr\\.snvs_indels\\.tiers\\.tsv$", "PcgrTiersTsvFile",
  "purple\\.cnv\\.gene\\.tsv$", "PurpleCnvGeneTsvFile",
  "manta\\.tsv$", "MantaTsvFile",
  "mapping_metrics\\.csv$", "MapMetricsFile"
)

# melt gds_indir columns to get a single column with paths to gds directories
# of interest, and fish out files of interest from each of them and grab
# presigned URLs
urls <- pmeta_tidy |>
  tidyr::pivot_longer(contains("gds_indir"), names_to = "folder_type", values_to = "gds_indir") |>
  select("SubjectID", "LibraryID", "rnasum_dataset", "folder_type", "gds_indir") |>
  rowwise() |>
  mutate(
    rnasum_files = list(
      dracarys::gds_files_list_filter_relevant(
        gdsdir = gds_indir, token = token_prod, include = "PresignedUrl",
        page_size = 1000, regexes = rnasum_file_regex
      )
    )
  )

d <- urls |>
  tidyr::unnest(rnasum_files) |>
  select("SubjectID", "LibraryID", dataset = "rnasum_dataset", ftype = "type", url = "presigned_url") |>
  tidyr::pivot_wider(names_from = ftype, values_from = url)

for (i in seq_len(nrow(d))) {
  v <- "0.6.5"
  date_str <- Sys.time() |>
    stringr::str_replace_all("[- :]", "") |>
    stringr::str_replace("\\..*", "")
  report_dir <- glue("{d$SubjectID[i]}_{d$LibraryID[i]}")
  username <- "pdiakumis"
  gds_path <- glue("gds://development/{username}/temp/rnasum/{v}/{report_dir}/{date_str}")
  workdir <- glue("{gds_path}/work")
  outdir <- glue("{gds_path}/output")
  jbody <- json_obj(
    rnasum_version = v,
    workdir = workdir,
    outdir = outdir,
    arriba_pdf = d$ArribaFusionsPdfFile[i],
    arriba_tsv = d$ArribaFusionsTsvFile[i],
    dataset = d$dataset[i],
    dragen_fusions = d$DragenWtsFusionsFinalFile[i],
    dragen_mapping_metrics = d$MapMetricsFile[i],
    manta_tsv = d$MantaTsvFile[i],
    pcgr_tiers_tsv = d$PcgrTiersTsvFile[i],
    purple_gene_tsv = d$PurpleCnvGeneTsvFile[i],
    salmon = d$DragenWtsQuantGenesSfFile[i],
    sample_name = d$SubjectID[i],
    subject_id = d$SubjectID[i],
    library_id = d$LibraryID[i]
  )
  launch_rnasum(json_body = jbody, token = token_dev)
}
