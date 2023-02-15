#!/usr/bin/env Rscript

################################################################################
#   Authors: Jacek Marzec ( jacek.marzec@unimelb.edu.au ), Sehrish Kanwal, Lavinia Gordon
#   University of Melbourne Centre for Cancer Research,
#   Victorian Comprehensive Cancer Centre
#   305 Grattan St, Melbourne, VIC 3000
################################################################################

suppressMessages(library(optparse))
suppressMessages(library(RNAsum))
suppressMessages(library(glue))

option_list <- list(
  make_option("--arriba_pdf", type = "character", help = "File path of Arriba PDF output"),
  make_option("--arriba_tsv", type = "character", help = "File path of Arriba TSV output"),
  make_option("--batch_rm", default = TRUE, type = "logical", help = "Remove batch-associated effects between datasets? [def: %default]"),
  make_option("--clinical_id", type = "character", help = "ID required to match sample with the subject clinical information"),
  make_option("--clinical_info", type = "character", help = "File path to clinical information xlsx file"),
  make_option("--cn_gain", default = 95, type = "integer", help = "CN threshold value to classify genes within gained regions [def: %default]"),
  make_option("--cn_loss", default = 5, type = "integer", help = "CN threshold value to classify genes within lost regions [def: %default]"),
  make_option("--dataset", default = "PANCAN", type = "character", help = "Dataset to be used as external reference cohort [def: %default]"),
  make_option("--dataset_name_incl", default = FALSE, type = "logical", help = "Include dataset in report name? [def: %default]"),
  make_option("--dragen_fusions", type = "character", help = "File path to DRAGEN RNA-seq 'fusion_candidates.final' output"),
  make_option("--drugs", default = FALSE, type = "logical", help = "Include drug matching section in report? [def: %default]"),
  make_option("--filter", default = TRUE, type = "logical", help = "Filter out low expressed genes? [def: %default]"),
  make_option("--hide_code_btn", default = TRUE, type = "logical", help = "Hide \"Code\" button above code chunks in report? [def: %default]"),
  make_option("--immunogram", default = FALSE, type = "logical", help = "Include immunogram in report? [def: %default]"),
  make_option("--log", default = TRUE, type = "logical", help = "Log2 transform data before normalisation? [def: %default]"),
  make_option("--manta_tsv", type = "character", help = "File path to umccrise 'manta.tsv' output"),
  make_option("--norm", type = "character", help = "Normalisation method"),
  make_option("--pcgr_splice_vars", default = TRUE, type = "logical", help = "Include non-coding splice region variants reported in PCGR? [def: %default]"),
  make_option("--pcgr_tier", default = 4, type = "integer", help = "Tier threshold for reporting variants reported in PCGR [def: %default]"),
  make_option("--pcgr_tiers_tsv", type = "character", help = "File path to PCGR 'snvs_indels.tiers.tsv' output"),
  make_option("--project", type = "character", help = "Project name"),
  make_option("--purple_gene_tsv", type = "character", help = "File path to PURPLE 'purple.cnv.gene.tsv' output"),
  make_option("--report_dir", type = "character", help = "Directory path to output report"),
  make_option("--salmon", type = "character", help = "File path to salmon 'quant.sf' output"),
  make_option("--sample_name", type = "character", help = "Sample name to be presented in report"),
  make_option("--sample_source", default = "-", type = "character", help = "Type of investigated sample [def: %default]"),
  make_option("--save_tables", default = TRUE, type = "logical", help = "Save interactive summary tables as HTML? [def: %default]"),
  make_option("--scaling", default = "gene-wise", type = "character", help = "Scaling for z-score transformation (gene-wise or group-wise) [def: %default]"),
  make_option("--subject_id", type = "character", help = "Subject ID"),
  make_option("--top_genes", default = 5, type = "integer", help = "Number of top ranked genes to be presented in report"),
  make_option("--transform", default = "CPM", type = "character", help = "Transformation method to be used when converting read counts [def: %default]"),
  make_option("--umccrise", type = "character", help = "Directory path of the corresponding WGS-related data")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list, formatter = optparse::TitledHelpFormatter))

##### Check required args
stopifnot(
  "'--sample_name' and '--report_dir' are required" =
    !is.null(opt$sample_name) && !is.null(opt$report_dir)
)
stopifnot(
  "'--clinical_info' requires at least one of '--clinical_id', '--subject_id', or '--umccrise'" =
    !is.null(opt$clinical_info) && (!is.null(opt$clinical_id) || !is.null(opt$subject_id) || !is.null(opt$umccrise))
)
stopifnot(
  "Invalid '--dataset'. Please see https://umccr.github.io/RNAsum/articles/tcga_projects_summary.html" =
    toupper(opt$dataset) %in% names(RNAsum::REFERENCE_DATASETS)
)

##### Make sure that TMM, TMMwzp, RLE or upperquartile normalisation is used for
##### CPM-tansformed data and quantile normalisation is used for TPM-tansformed data
stopifnot("'--transform' must be CPM or TPM" = opt$transform %in% c("CPM", "TPM"))
if (is.null(opt$norm)) {
  if (opt$transform == "CPM") {
    opt$norm <- "TMM"
  } else if (opt$transform == "TPM") {
    opt$norm <- "quantile"
  }
}

if (opt$transform == "TPM" && !(opt$norm %in% c("quantile", "none"))) {
  warning(
    "'--norm ", opt$norm, "' normalisation is not available for TPM-transformed data ('--transform TPM')!\n",
    "Quantile normalisation ('--norm quantile') will be performed instead."
  )
  opt$norm <- "quantile"
} else if (opt$transform == "CPM" && opt$norm == "quantile") {
  stop(
    "Quantile normalisation ('--norm quantile') is available only for TPM-transformed data ('--transform TPM')!\n",
    "Use '--norm' with one of TMM, TMMwzp, RLE, upperquartile or none when using '--transform CPM'."
  )
} else if (opt$transform == "CPM" && !(opt$norm %in% c("TMM", "TMMwzp", "RLE", "upperquartile", "none"))) {
  stop(
    "Use '--norm' with one of TMM, TMMwzp, RLE, upperquartile or none when using '--transform CPM'."
  )
}

RNAsum::mkdir(opt$report_dir)

# Render Rmd report with above parameters
rmd_file <- system.file("rmd/rnasum.Rmd", package = "RNAsum")
dataset_name_incl <- ifelse(opt$dataset_name_incl, glue::glue("_{opt$dataset}"), "")
out_file_base <- glue::glue("{opt$sample_name}{toupper(dataset_name_incl)}.RNAseq_report")
rmarkdown::render(
  input = rmd_file,
  output_file = glue::glue("{out_file_base}.html"),
  output_dir = opt$report_dir,
  params = opt
)

##### Remove the associated MD file and the redundant folder with plots that are imbedded in the HTML report
unlink(file.path(opt$report_dir, glue::glue("{out_file_base}.md")))
unlink(file.path(opt$report_dir, glue::glue("{out_file_base}_files")), recursive = TRUE)
