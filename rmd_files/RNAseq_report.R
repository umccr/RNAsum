################################################################################
#
#   File name: RNAseq_report.R
#
#   Authors: Jacek Marzec ( jacek.marzec@unimelb.edu.au ), Sehrish Kanwal, Lavinia Gordon
#
#   University of Melbourne Centre for Cancer Research,
#   Victorian Comprehensive Cancer Centre
#   305 Grattan St, Melbourne, VIC 3000
#
################################################################################

################################################################################
#
#	  Description: Script collecting user-defined parameters for the corresponding RNAseq_report.Rmd markdown script generating the "UMCCR Transcriptome Patient Summary" report. Note, only genes intersection between the sample read count file and the reference datasets expression matrices will be considered in the analyses.
#
#	  Command line use example: Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset PAAD  --bcbio_rnaseq $(pwd)/../data/test_data/final/test_sample_WTS  --report_dir (pwd)/../data/test_data/final/test_sample_WTS/RNAsum  --umccrise $(pwd)/../data/test_data/umccrised/test_sample_WGS  --clinical_info $(pwd)/../data/test_data/test_clinical_data.xlsx  --clinical_id test.subject
#
#   sample_name:  The name of the sample to be analysed and reported
#   bcbio_rnaseq: Location of the results folder from bcbio RNA-seq pipeline
#   dragen_rnaseq:Location of the results folder from dragen RNA-seq pipeline
#   dataset:      Dataset to be used as external reference cohort (default is "PANCAN")
#   report_dir:   Desired location for the report
#   ref_data_dir: Location of the reference and annotation files
#   transform:    Transformation method to be used when converting read counts. Available options are: "CPM" (default) and "TPM"
#   norm:         Normalisation method. Currently, "TMM","TMMwzp", "RLE" and "upperquartile" methods are available for CPM-transformed data and "quantile" normalisation is used for TPM-transformed data
#   batch_rm:     Remove batch-associated effects between datasets. Available options are: "TRUE" (default) and "FALSE"
#   filter:       Filtering out low expressed genes. Available options are: "TRUE" (default) and "FALSE"
#   log:          Log (base 2) transform data before normalisation. Available options are: "TRUE" (default) and "FALSE"
#   scaling:      Apply "gene-wise" (default) or "group-wise" data scaling
#   drugs:        Include drug matching section in the report. Available options are: "TRUE" and "FALSE" (default)
#   immunogram:   Include immunogram in the report. Available options are: "TRUE" and "FALSE" (default for now, as this feature is not finilised yet)
#   umccrise (optional):  Location of the corresponding umccrise output from genomic-related data (including PCGR, PURPLE and Manta output files)
#   pcgr_tier (optional): Tier threshold for reporting variants reported in PCGR (default is "4")
#   pcgr_splice_vars (optional): Include non-coding splice region variants reported in PCGR. Available options are: "TRUE" (default) and "FALSE"
#   cn_loss (optional):  CN threshold value to classify genes within lost regions (default is "5th percentile" of all CN values)
#   cn_gain (optional):  CN threshold value to classify genes within gained regions (default is "95th percentile" of all CN values)
#   clinical_info (optional): Location of xslx file with clinical information
#   clinical_id (optional):   ID required to match sample with the subject clinical information (specified in flag --clinical_info)
#   subject_id (optional):    Subject ID. If umccrise output is specified (flag --umccrise) then Subject ID is extracted from there and used to overwrite this argument
#   sample_source (optional):   Source of investigated sample (e.g. fresh frozen tissue, organoid). This information is for annotation purposes only
#   dataset_name_incl:  Include dataset in the report and sample name. Available options are: "TRUE" and "FALSE" (default)
#   project (optional):   Project name. This information is for annotation purposes only
#   top_genes:     The number of top ranked genes to be presented (default is "5")
#   save_tables:   Save interactive summary tables as HTML. Available options are: "TRUE" (default) and "FALSE"
#   hide_code_btn: Hide the "Code" button allowing to show/hide code chunks in the final HTML report. Available options are: "TRUE" (default) and "FALSE"
#   grch_version:  Human reference genome version used for genes annotation (default is "38")
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Load libraries
#===============================================================================

suppressMessages(library(optparse))
suppressMessages(library(glue))

#===============================================================================
#    Define functions
#===============================================================================

##### Create 'not in' operator
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option("--sample_name", action="store", default=NA, type='character',
              help="Desired sample name to be presented in the report"),
  make_option("--dataset", action="store", default="PANCAN", type='character',
              help="Dataset to be used as external reference cohort"),
  make_option("--bcbio_rnaseq", action="store", default=NULL, type='character',
              help="Location of the results folder from bcbio RNA-seq pipeline"),
  make_option("--dragen_rnaseq", action="store", default=NULL, type='character',
              help="Location of the results folder from Dragen RNA-seq pipeline"),
  make_option("--report_dir", action="store", default=NA, type='character',
              help="Desired location for the report"),
  make_option("--ref_data_dir", action="store", default="../data", type='character',
              help="Location of the reference and annotation files"),
  make_option("--transform", action="store", default="CPM", type='character',
              help="Transformation method to be used when converting read counts"),
  make_option("--norm", action="store", default=NA, type='character',
              help="Normalisation method"),
  make_option("--batch_rm", action="store", default=TRUE, type='logical',
              help="Remove batch-associated effects between datasets"),
  make_option("--filter", action="store", default=TRUE, type='logical',
              help="Filtering out low expressed genes"),
  make_option("--log", action="store", default=TRUE, type='logical',
              help="Log (base 2) transform data before normalisation"),
  make_option("--scaling", action="store", default="gene-wise", type='character',
              help="Scaling for z-score transformation (gene-wise or group-wise"),
  make_option("--drugs", action="store", default=FALSE, type='logical',
              help="Include drug matching section in the report"),
  make_option("--immunogram", action="store", default=FALSE, type='logical',
              help="Include immunogram in the report"),
  make_option("--umccrise", action="store", default=NULL, type='character',
              help="Location of the corresponding WGS-related data"),
  make_option("--pcgr_tier", action="store", default=4, type='integer',
              help="Tier threshold for reporting variants reported in PCGR"),
  make_option("--pcgr_splice_vars", action="store", default=TRUE, type='logical',
              help="Include non-coding splice region variants reported in PCGR"),
  make_option("--cn_loss", action="store", default=5, type='integer',
              help="CN threshold value to classify genes within lost regions"),
  make_option("--cn_gain", action="store", default=95, type='integer',
              help="CN threshold value to classify genes within gained regions"),
  make_option("--clinical_info", action="store", default=NA, type='character',
              help="Location of xslx file with clinical information"),
  make_option("--clinical_id", action="store", default=NA, type='character',
              help="ID required to match sample with the subject clinical information"),
  make_option("--subject_id", action="store", default=NA, type='character',
              help="Subject ID"),
  make_option("--sample_source", action="store", default="-", type='character',
              help="Type of investigated sample"),
  make_option("--dataset_name_incl", action="store", default=NA, type='character',
              help="Include dataset in the report name"),
  make_option("--project", action="store", default="-", type='character',
              help="Project name"),
  make_option("--top_genes", action="store", default=5, type='integer',
              help="The number of top ranked genes to be presented"),
  make_option("--save_tables", action="store", default=TRUE, type='logical',
              help="Save interactive summary tables as HTML"),
  make_option("--hide_code_btn", action="store", default=TRUE, type='logical',
              help="Hide the \"Code\" button allowing to show/hide code chunks in the final HTML report"),
  make_option("--grch_version", action="store", default=NA, type='integer',
              help="human reference genome version used for genes annotation")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in argument from command line and check if all were provide by the user
if ( (is.na(opt$sample_name) || is.null(opt$bcbio_rnaseq) || is.na(opt$report_dir)) && is.null(opt$dragen_rnaseq) ) {

  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset PAAD  --bcbio_rnaseq $(pwd)/../data/test_data/final/test_sample_WTS  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report\n\nor\n\n")
  cat("Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset PAAD  --dragen_rnaseq $(pwd)/../data/test_data/stratus/test_sample_WTS  --report_dir $(pwd)/../data/test_data/stratus/test_sample_WTS/RNAseq_report\n\n")
  q()
  
} else if ( !is.null(opt$bcbio_rnaseq) && !is.null(opt$dragen_rnaseq) ) {
  
  cat("\nOutput from only one RNA-seq pipeline, either bcbio-nextgen RNA-seq pipeline or DRAGEN RNA pipeline, is accepted at a time!\n\n")
  cat("\ncommand example:\n\nRscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset PAAD  --bcbio_rnaseq $(pwd)/../data/test_data/final/test_sample_WTS  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report\n\nor\n\n")
  cat("Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset PAAD  --dragen_rnaseq $(pwd)/../data/test_data/stratus/test_sample_WTS  --report_dir $(pwd)/../data/test_data/stratus/test_sample_WTS/RNAseq_report\n\n")
  q()
}

##### Make sure that sample ID is availabe if clincal data is provided
if ( !is.na(opt$clinical_info) && any(is.na(opt$clinical_id) && is.na(opt$subject_id) && is.na(opt$umccrise) ) ) {

  cat("ID required to match sample with the subject clinical information is missing! Please provide the ID used in the clinical data by using \"--clinical_id\" argument.\n\n")
  q()
}

##### Set default parameters
if ( is.na(opt$norm)  ) {

  if ( opt$transform == "CPM"  ) {
    opt$norm <- "TMM"

  } else if ( opt$transform == "TPM"  ) {
    opt$norm <- "quantile"
  }
}

if ( is.na(opt$dataset_name_incl)  ) {
  dataset_name_incl <- ""
} else if ( isFALSE(as.logical(opt$dataset_name_incl))  ) {
  dataset_name_incl <- ""
} else {
  dataset_name_incl <- paste0("_", opt$dataset)
}

if ( is.na(opt$grch_version)  ) {
  opt$grch_version <- 38
  ensembl_version <- 86
  ucsc_genome_assembly <- 38
  
} else if ( opt$grch_version == 38 ) {
  ensembl_version <- 86
  ucsc_genome_assembly <- 38
  
} else if ( opt$grch_version == 37 ) {
  ensembl_version <- 75
  ucsc_genome_assembly <- 19
  
} else {
  cat("\nCurrently human reference genome (GRCh) versions \"37\" and \"38\" are supported.\n\n")
  q()
}

##### Check if specified dataset type is valid
if ( toupper(opt$dataset) %!in% toupper(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "BLCA-NET", "PAAD-IPMN", "PAAD-NET", "PAAD-ACC", "LUAD-LCNEC", "PANCAN", "TEST")) ) {

  cat("\nInvalid dataset! Please use one of the following:\n\n")
  cat("[ ACC ] - this will compare the patient's data in the context of samples from TCGA Adrenocortical Carcinoma cohort\n\n")
  cat("[ BLCA ] - this will compare the patient's data in the context of samples from TCGA Bladder Urothelial Carcinoma cohort\n\n")
  cat("[ BRCA ] - this will compare the patient's data in the context of samples from TCGA Breast Invasive Carcinoma cohort\n\n")
  cat("[ CESC ] - this will compare the patient's data in the context of samples from TCGA Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma cohort\n\n")
  cat("[ CHOL ] - this will compare the patient's data in the context of samples from TCGA Cholangiocarcinoma cohort\n\n")
  cat("[ COAD ] - this will compare the patient's data in the context of samples from TCGA Colon Adenocarcinoma cohort\n\n")
  cat("[ DLBC ] - this will compare the patient's data in the context of samples from TCGA Lymphoid Neoplasm Diffuse Large B-cell Lymphoma cohort\n\n")
  cat("[ ESCA ] - this will compare the patient's data in the context of samples from TCGA Esophageal Carcinoma cohort\n\n")
  cat("[ GBM ] - this will compare the patient's data in the context of samples from TCGA Glioblastoma Multiforme cohort\n\n")
  cat("[ HNSC ] - this will compare the patient's data in the context of samples from TCGA Head and Neck Squamous Cell Carcinoma cohort\n\n")
  cat("[ KICH ] - this will compare the patient's data in the context of samples from TCGA Kidney Chromophobe cohort\n\n")
  cat("[ KIRC ] - this will compare the patient's data in the context of samples from TCGA Kidney Renal Clear Cell Carcinoma cohort\n\n")
  cat("[ KIRP ] - this will compare the patient's data in the context of samples from TCGA Kidney Renal Papillary Cell Carcinoma cohort\n\n")
  cat("[ LAML ] - this will compare the patient's data in the context of samples from TCGA Acute Myeloid Leukaemia cohort\n\n")
  cat("[ LGG ] - this will compare the patient's data in the context of samples from TCGA Brain Lower Grade Glioma cohort\n\n")
  cat("[ LIHC ] - this will compare the patient's data in the context of samples from TCGA Liver Hepatocellular Carcinoma cohort\n\n")
  cat("[ LUAD ] - this will compare the patient's data in the context of samples from TCGA Lung Adenocarcinoma cohort\n\n")
  cat("[ LUSC ] - this will compare the patient's data in the context of samples from TCGA Lung Squamous Cell Carcinoma cohort\n\n")
  cat("[ MESO ] - this will compare the patient's data in the context of samples from TCGA Mesothelioma cohort\n\n")
  cat("[ OV ] - this will compare the patient's data in the context of samples from TCGA Ovarian Serous Cystadenocarcinoma cohort\n\n")
  cat("[ PAAD ] - this will compare the patient's data in the context of samples from TCGA Pancreatic Adenocarcinoma cohort and UMCCR internal pancreatic ductal adenocarcinoma cohort\n")
  cat("[ PCPG ] - this will compare the patient's data in the context of samples from TCGA Pheochromocytoma and Paraganglioma cohort\n\n")
  cat("[ PRAD ] - this will compare the patient's data in the context of samples from TCGA Prostate Adenocarcinoma cohort\n\n")
  cat("[ READ ] - this will compare the patient's data in the context of samples from TCGA Rectum Adenocarcinoma cohort\n\n")
  cat("[ SARC ] - this will compare the patient's data in the context of samples from TCGA Sarcoma cohort\n\n")
  cat("[ SKCM ] - this will compare the patient's data in the context of samples from TCGA Skin Cutaneous Melanoma cohort\n\n")
  cat("[ STAD ] - this will compare the patient's data in the context of samples from TCGA Stomach Adenocarcinoma cohort\n\n")
  cat("[ TGCT ] - this will compare the patient's data in the context of samples from TCGA Testicular Germ Cell Tumours cohort\n\n")
  cat("[ THCA ] - this will compare the patient's data in the context of samples from TCGA Thyroid Carcinoma cohort\n\n")
  cat("[ THYM ] - this will compare the patient's data in the context of samples from TCGA Thymoma cohort\n\n")
  cat("[ UCEC ] - this will compare the patient's data in the context of samples from TCGA Uterine Corpus Endometrial Carcinoma cohort\n\n")
  cat("[ UCS ] - this will compare the patient's data in the context of samples from TCGA Uterine Carcinosarcoma cohort\n\n")
  cat("[ UVM ] - this will compare the patient's data in the context of samples from TCGA Uveal Melanoma cohort\n\n")
  cat("[ BLCA-NET ] - this will compare the patient's data in the context of samples from TCGA Bladder Urothelial Carcinoma cohort (including neuroendocrine tumour (NET) samples)\n\n")
  cat("[ PAAD-IPMN ] - this will compare the patient's data in the context of samples from TCGA Pancreatic Adenocarcinoma cohort (including intraductal papillary mucinous neoplasm (IPMN) samples) and UMCCR internal pancreatic ductal adenocarcinoma cohort\n")
  cat("[ PAAD-NET ] - this will compare the patient's data in the context of samples from TCGA Pancreatic Adenocarcinoma cohort (including neuroendocrine tumour (NET) samples) and UMCCR internal pancreatic ductal adenocarcinoma cohort\n")
  cat("[ PAAD-ACC ] - this will compare the patient's data in the context of samples from TCGA Pancreatic Adenocarcinoma cohort (including acinar cell carcinoma (ACC) samples) and UMCCR internal pancreatic ductal adenocarcinoma cohort\n")
  cat("[ LUAD-LCNEC ] - this will compare the patient's data in the context of samples from TCGA Lung Adenocarcinoma cohort (including favor large-cell neuroendocrine carcinoma (LCNEC) samples)\n")
  cat("[ PANCAN ] - this will compare the patient's data in the context of samples from all TCGA cohorts\n")
  cat("[ TEST ] - this will compare the patient's data in the context of subset of samples from TCGA Pancreatic Adenocarcinoma cohort and UMCCR internal pancreatic ductal adenocarcinoma cohort\n")
  q()
}

##### Make sure that either CPM or TPM transformation is selected
if ( opt$transform != "CPM" && opt$transform != "TPM" ) {

  cat(paste0("\nWrong transformation method was selected! \"", opt$transform, "\" transformation is not available!\n\nUse \"CPM\" or \"TPM\" transformation methods.\n\n"))
  q()
}

##### Make sure that TMM, TMMwzp, RLE or upperquartile normalisation is used for CPM-tansformed data and quantile normalisation is used for TPM-tansformed data
if ( opt$transform == "TPM" && opt$norm != "quantile" && opt$norm != "none" ) {

  cat(paste0("\nWrong normalisation method was selected! ", opt$norm, " normalisation is not available for TPM-tansformed data!\n\nQuantile normalisation will be performed for TPM-tansformed data.\n\n"))
  opt$norm <- "quantile"

} else if ( opt$transform == "CPM" && opt$norm == "quantile" ) {

  cat(paste0("\nQuantile normalisation is available only for TPM-tansformed data! \"TMM\", \"TMMwzp\", \"RLE\", \"upperquartile\" or \"none\" methods are available for CPM-tansformed data.\n\n"))
  q()

} else if ( opt$transform == "CPM" && opt$norm != "TMM" && opt$norm != "TMMwzp" && opt$norm != "RLE" && opt$norm != "upperquartile" && opt$norm != "none" ) {

  cat(paste0("\nWrong normalisation method was selected! \"TMM\", \"TMMwzp\", \"RLE\", \"upperquartile\" or \"none\" methods are available for CPM-tansformed data.\n\n"))
  q()
}

##### Create user-defined directory for the report
if ( !file.exists(opt$report_dir) ) {
  dir.create(opt$report_dir, recursive=TRUE)
}

##### Collect parameters
param_list <- list(sample_name = opt$sample_name,
               dataset = toupper(opt$dataset),
               bcbio_rnaseq = opt$bcbio_rnaseq,
               dragen_rnaseq = opt$dragen_rnaseq,
               report_dir = opt$report_dir,
               ref_data_dir = opt$ref_data_dir,
               transform = opt$transform,
               norm = opt$norm,
               batch_rm = opt$batch_rm,
               filter = opt$filter,
               log = opt$log,
               scaling = opt$scaling,
               drugs = opt$drugs,
               immunogram = opt$immunogram,
               umccrise = opt$umccrise,
               clinical_info = opt$clinical_info,
               clinical_id = opt$clinical_id,
               subject_id = opt$subject_id,
               sample_source = opt$sample_source,
               dataset_name_incl = dataset_name_incl,
               project = opt$project,
               save_tables = opt$save_tables,
               pcgr_tier = opt$pcgr_tier,
               pcgr_splice_vars = opt$pcgr_splice_vars,
               cn_loss = opt$cn_loss,
               cn_gain = opt$cn_gain,
               top_genes = opt$top_genes,
               hide_code_btn = opt$hide_code_btn,
               grch_version = as.numeric(opt$grch_version),
               ensembl_version = as.numeric(ensembl_version),
               ucsc_genome_assembly = as.numeric(ucsc_genome_assembly))

##### Pass the user-defined arguments to the RNAseq_report R markdown script and generate the report
rmarkdown::render(input = "RNAseq_report.Rmd",
                  output_file = paste0(opt$sample_name, toupper(dataset_name_incl), ".RNAseq_report.html"),
                  output_dir = opt$report_dir,
                  params = param_list )

##### Remove the assocaited MD file and the redundant folder with plots that are imbedded in the HTML report
unlink(paste0(opt$report_dir, "/", opt$sample_name, toupper(dataset_name_incl), ".RNAseq_report.md"), recursive = TRUE)
unlink(paste0(opt$report_dir, "/", opt$sample_name, toupper(dataset_name_incl), ".RNAseq_report_files"), recursive = TRUE)

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
