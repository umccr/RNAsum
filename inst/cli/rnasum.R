#!/usr/bin/env Rscript

################################################################################
#   Authors: Jacek Marzec ( jacek.marzec@unimelb.edu.au ), Sehrish Kanwal, Lavinia Gordon
#   University of Melbourne Centre for Cancer Research,
#   Victorian Comprehensive Cancer Centre
#   305 Grattan St, Melbourne, VIC 3000
################################################################################

suppressMessages(library(optparse, include.only = "make_option"))
suppressMessages(library(RNAsum))
suppressMessages(library(glue, include.only = "glue"))
suppressMessages(library(fs, include.only = "dir_create"))

option_list <- list(
  make_option("--arriba_pdf", type = "character", help = "File path of Arriba PDF output."),
  make_option("--arriba_tsv", type = "character", help = "File path of Arriba TSV output."),
  make_option("--batch_rm", action = "store_true", help = "Remove batch-associated effects between datasets."),
  make_option("--cn_gain", default = 95, type = "integer", help = "CN threshold value to classify genes within gained regions. [def: %default]"),
  make_option("--cn_loss", default = 5, type = "integer", help = "CN threshold value to classify genes within lost regions. [def: %default]"),
  make_option("--dataset", default = "PANCAN", type = "character", help = "Dataset to be used as external reference cohort. [def: %default]"),
  make_option("--dataset_name_incl", action = "store_true", help = "Include dataset in report name."),
  make_option("--dragen_fusions", type = "character", help = "File path to DRAGEN RNA-seq 'fusion_candidates.final' output."),
  make_option("--dragen_mapping_metrics", type = "character", help = "File path to DRAGEN RNA-seq 'mapping_metrics.csv' output."),
  make_option("--drugs", action = "store_true", help = "Include drug matching section in report."),
  make_option("--filter", action = "store_true", help = "Filter out low expressed genes."),
  make_option("--immunogram", action = "store_true", help = "Include immunogram in report."),
  make_option("--log", action = "store_true", default = TRUE, help = "Log2 transform data before normalisation."),
  make_option("--manta_tsv", type = "character", help = "File path to umccrise 'manta.tsv' output."),
  make_option("--norm", type = "character", help = "Normalisation method."),
  make_option("--pcgr_splice_vars", action = "store_true", help = "Include non-coding splice region variants reported in PCGR."),
  make_option("--pcgr_tier", default = 4, type = "integer", help = "Tier threshold for reporting variants reported in PCGR. [def: %default]"),
  make_option("--pcgr_tiers_tsv", type = "character", help = "File path to PCGR 'snvs_indels.tiers.tsv' output."),
  make_option("--project", type = "character", help = "Project name, used for annotation purposes only."),
  make_option("--purple_gene_tsv", type = "character", help = "File path to PURPLE 'purple.cnv.gene.tsv' output."),
  make_option("--report_dir", type = "character", help = "Directory path to output report."),
  make_option("--salmon", type = "character", help = "File path to salmon 'quant.genes.sf' output."),
  make_option("--sample_name", type = "character", help = "Sample name to be presented in report."),
  make_option("--sample_source", default = "-", type = "character", help = "Type of investigated sample. [def: %default]"),
  make_option("--save_tables", action = "store_true", help = "Save interactive summary tables as HTML."),
  make_option("--scaling", default = "gene-wise", type = "character", help = "Scaling for z-score transformation (gene-wise or group-wise). [def: %default]"),
  make_option("--subject_id", type = "character", help = "Subject ID."),
  make_option("--top_genes", default = 5, type = "integer", help = "Number of top ranked genes to be presented in report."),
  make_option("--transform", default = "CPM", type = "character", help = "Transformation method to be used when converting read counts. [def: %default]"),
  make_option("--umccrise", type = "character", help = "Directory path of the corresponding WGS-related umccrise data."),
  make_option(c("--version", "-v"), action = "store_true", help = "Print RNAsum version and exit."),
  make_option("--html_dir", type = "character", default = getwd(), help = "Directory path to output final HTML report. [def: current directory (%default)]")
)

parser <- optparse::OptionParser(option_list = option_list, formatter = optparse::TitledHelpFormatter)
opt <- optparse::parse_args(parser)

if (!is.null(opt$version)) {
  cat(as.character(packageVersion("RNAsum")), "\n")
  quit("no", status = 0, runLast = FALSE)
}
# don't need these any more so NULLify to remove from params
opt$version <- NULL
opt$help <- NULL
html_dir <- opt$html_dir
opt$html_dir <- NULL

##### Convert missing flags to FALSE
flags <- c("batch_rm", "dataset_name_incl", "drugs", "filter", "immunogram", "log", "pcgr_splice_vars", "save_tables")
for (flag in flags) {
  if (is.null(opt[[flag]])) {
    opt[[flag]] <- FALSE
  }
}

##### Check required args
if (is.null(opt$sample_name) || is.null(opt$report_dir)) {
  cat("'--sample_name' and '--report_dir' are required!\n")
  quit("no", status = 1, runLast = FALSE)
}
if (!toupper(opt$dataset) %in% names(RNAsum::REFERENCE_DATASETS)) {
  cat("Invalid '--dataset'. Please see https://umccr.github.io/RNAsum/articles/tcga_projects_summary.html")
  quit("no", status = 1, runLast = FALSE)
}

##### Make sure that TMM, TMMwzp, RLE or upperquartile normalisation is used for
##### CPM-tansformed data and quantile normalisation is used for TPM-tansformed data
stopifnot("'--transform' must be CPM or TPM" = opt$transform %in% c("CPM", "TPM"))
if (is.null(opt$norm)) {
  if (opt$transform == "CPM") {
    opt$norm <- "TMM"
  } else if (opt$transform == "TPM") {
    opt$norm <- "quantile"
  } else {
    stop("YOU SHOULD NOT GET HERE!")
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

# Render Rmd report with above parameters
dataset_name_incl <- ifelse(opt$dataset_name_incl, glue::glue("_{opt$dataset}"), "")
out_file_base <- glue::glue("{opt$sample_name}{toupper(dataset_name_incl)}.RNAseq_report")
out_html <- RNAsum::rnasum_rmd(
  out_file = file.path(html_dir, glue::glue("{out_file_base}.html")),
  pars = opt
)
message(glue::glue("RNAsum HTML output at: {out_html}"))
