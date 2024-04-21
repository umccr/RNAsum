query_workflow_rnasum <- function(sbjids) {
  wfrnms <- glue("umccr__automated__rnasum__{sbjids}") |>
    paste(collapse = "|")
  q1 <- glue(
    "WHERE \"type_name\" = 'rnasum' AND REGEXP_LIKE(\"wfr_name\", '{wfrnms}') ",
    "ORDER BY \"start\" DESC;"
  )
  rportal::portaldb_query_workflow(q1)
}

json_obj <- function(rnasum_version, workdir, outdir, arriba_pdf, arriba_tsv, dataset,
                     dragen_fusions, dragen_mapping_metrics,
                     manta_tsv, pcgr_tiers_tsv, purple_gene_tsv,
                     salmon, sample_name, subject_id, library_id) {
  report_dir <- glue("{subject_id}_{library_id}")
  fc <- function(f) {
    list(
      class = "File",
      path = f
    )
  }
  list(
    "Name" = glue("rnasum_{rnasum_version}_{report_dir}"),
    "Input" = list(
      "arriba_pdf" = fc(arriba_pdf),
      "arriba_tsv" = fc(arriba_tsv),
      "dataset" = dataset,
      "dragen_fusions" = fc(dragen_fusions),
      "dragen_mapping_metrics" = fc(dragen_mapping_metrics),
      "manta_tsv" = fc(manta_tsv),
      "pcgr_tiers_tsv" = fc(pcgr_tiers_tsv),
      "purple_gene_tsv" = fc(purple_gene_tsv),
      "report_dir" = report_dir,
      "salmon" = fc(salmon),
      "sample_name" = sample_name,
      "subject_id" = subject_id
    ),
    "engineParameters" = list(
      "workDirectory" = workdir,
      "outputDirectory" = outdir,
      "overrides" = list(
        "#main" = list(
          "requirements" = list(
            "DockerRequirement" = list(
              "dockerPull" = glue("ghcr.io/umccr/rnasum:{rnasum_version}")
            )
          )
        )
      )
    )
  )
}

launch_rnasum <- function(json_body, token, wfl_id = "wfl.87e07ae6b46645a181e04813de535216", wfl_vname = "0.5.0--304c3ac") {
  token <- dracarys::ica_token_validate(token)
  base_url <- "https://aps2.platform.illumina.com/v1"
  httr::POST(
    url = glue("{base_url}/workflows/{wfl_id}/versions/{wfl_vname}:launch"),
    config = httr::add_headers(Authorization = glue::glue("Bearer {token}")),
    body = json_body,
    encode = "json"
  )
}
