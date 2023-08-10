# helpers for downloading GDS data for local testing
require(dracarys)
require(tidyverse)
require(glue, include.only = "glue")
require(here, include.only = "here")

glims_read <- function() {
  lims_key <- googledrive::drive_find("^Google LIMS$", shared_drive = "LIMS")$id
  lims <- lims_key |>
    googlesheets4::read_sheet("Sheet1", na = c(".", "", "-"), col_types = "c")
  lims |> readr::type_convert(col_types = readr::cols(.default = "c", Timestamp = "T"))
}
lims_rds <- here::here(glue("nogit/data_portal/lims.rds"))
# lims_raw <- glims_read()
# saveRDS(lims_raw, file = lims_rds)
lims_raw <- readr::read_rds(lims_rds)
pmeta_rds <- here::here(glue("nogit/data_portal/workflows.rds"))
# grab last X RNAsum runs via portal API
# pmeta_raw <- dracarys::portal_meta_read(params = "&type_name=rnasum", rows = 10)
# saveRDS(pmeta_raw, file = pmeta_rds)
pmeta_raw <- readr::read_rds(pmeta_rds)

lims <- lims_raw |>
  dplyr::select(
    Timestamp, SubjectID, SampleID, SampleName, LibraryID, ExternalSubjectID, ExternalSampleID,
    ProjectOwner, ProjectName, Type, Assay, Phenotype, Source, Quality, Topup, Workflow
  )
pmeta <- pmeta_raw

# generate tidy rnasum metadata from portal workflows table, and join against glims
meta_rnasum <- pmeta |>
  dracarys::meta_rnasum(status = c("Succeeded", "Failed")) |>
  dplyr::left_join(lims, by = c("LibraryID", "SampleID", "SubjectID")) |>
  dplyr::select(
    gds_indir_dragen, gds_indir_umccrise, gds_indir_arriba,
    SubjectID, LibraryID, SampleID, Phenotype, rnasum_dataset,
    end_status,
    # ExternalSubjectID, ProjectOwner, ProjectName, Type, Assay, Source, Quality, Workflow,
    wfr_id, start, end, gds_outfile_rnasum_html,
  ) |>
  dplyr::arrange(desc(SubjectID), start)

# melt gds_indir columns to get a single column with paths to gds directories
# of interest, so that we can fish out files of interest from each of them
meta_rnasum <- meta_rnasum |>
  tidyr::pivot_longer(dplyr::contains("gds_indir"), names_to = "folder_type", values_to = "gds_indir") |>
  dplyr::select(SubjectID, LibraryID, rnasum_dataset, folder_type, gds_indir) |>
  dplyr::mutate(
    folder_type = sub("gds_indir_", "", folder_type),
    local_outdir = here::here(glue::glue("nogit/patient_data/{SubjectID}/{LibraryID}/{rnasum_dataset}/{folder_type}")),
  )

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
)
token <- dracarys::ica_token_validate()
dryrun <- F

# download files locally
down <- meta_rnasum |>
  dplyr::rowwise() |>
  dplyr::mutate(
    down = list(dracarys::dr_gds_download(.data$gds_indir, token = token, page_size = 150, outdir = .data$local_outdir, dryrun = dryrun, regexes = rnasum_file_regex))
  ) |>
  dplyr::ungroup()

# saveRDS(down, here::here("nogit/patient_data/down_2023-08-09.rds"))
# down <- here::here("nogit/patient_data/down_2023-08-09.rds") |> readr::read_rds()

params_set <- function(arriba_pdf, arriba_tsv, dataset, dragen_fusions, manta_tsv,
                       pcgr_tiers_tsv, purple_gene_tsv, report_dir, salmon,
                       sample_name, subject_id) {
  params <- list(
    arriba_pdf = arriba_pdf,
    arriba_tsv = arriba_tsv,
    batch_rm = TRUE,
    cn_gain = 95,
    cn_loss = 5,
    dataset = dataset,
    dataset_name_incl = TRUE,
    dragen_fusions = dragen_fusions,
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
d <- down |>
  tidyr::unnest(down) |>
  dplyr::select(-c("file_id", "dname", "size", "path", "bname", "out_dl", "local_outdir", "gds_indir", "folder_type")) |>
  tidyr::pivot_wider(names_from = type, values_from = out)

# slice to whichever run you want from d
d |>
  dplyr::slice(3) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    params = list(
      params_set(
        arriba_pdf = ArribaFusionsPdfFile, arriba_tsv = ArribaFusionsTsvFile, dataset = rnasum_dataset,
        dragen_fusions = DragenWtsFusionsFinalFile, manta_tsv = MantaTsvFile,
        pcgr_tiers_tsv = PcgrTiersTsvFile, purple_gene_tsv = PurpleCnvGeneTsvFile,
        report_dir = here::here(glue::glue("nogit/patient_data/reports/{SubjectID}_{LibraryID}_{rnasum_dataset}")),
        salmon = DragenWtsQuantSfFile, sample_name = glue::glue("{SubjectID}_{LibraryID}"), subject_id = SubjectID
      )
    ),
    rnasum_rmd = rnasum_rmd(
      out_file = here::here(glue::glue("nogit/patient_data/reports_html/{SubjectID}_{LibraryID}_{rnasum_dataset}.html")),
      quiet = FALSE, pars = params
    )
  ) |>
  dplyr::ungroup()
