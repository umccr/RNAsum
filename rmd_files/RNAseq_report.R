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
#	  Command line use example: Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset PAAD  --bcbio_rnaseq $(pwd)/../data/test_data/final/test_sample_WTS  --report_dir (pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report  --umccrise $(pwd)/../data/test_data/umccrised/test_sample_WGS  --clinical_info $(pwd)/../data/test_data/test_clinical_data.xlsx  --subject_id test.subject
#
#   dataset:      Dataset to be used as external reference cohort
#   bcbio_rnaseq: Location of the results folder from bcbio RNA-seq pipeline
#   report_dir:   Desired location for the report
#   ref_data_dir: Location of the reference and annotation files
#   transform:    Transformation method to be used when converting read counts. Available options are: "CPM" (default) and "TPM"
#   norm:         Normalisation method. Currently, "TMM" is used for CPM-transformed data and "quantile" normalisation is used for TPM-transformed data
#   batch_rm:     Remove batch-associated effects between datasets. Available options are: "TRUE" (default) and "FALSE"
#   filter:       Filtering out low expressed genes. Available options are: "TRUE" (default) and "FALSE"
#   log:          Log (base 2) transform data before normalisation. Available options are: "TRUE" (default) and "FALSE"
#   scaling:      Apply "gene-wise" (default) or "group-wise" data scaling
#   umccrise (optional):  Location of the corresponding umccrise output from genomic-related data (including PCGR, PURPLE and Manta output files)
#   pcgr_tier (optional): Tier threshold for reporting variants reported in PCGR (default is "4")
#   pcgr_splice_vars (optional): Include non-coding splice region variants reported in PCGR. Available options are: "TRUE" (default) and "FALSE" 
#   cn_loss (optional):  CN threshold value to classify genes within lost regions (default is "5th percentile" of all CN values)
#   cn_gain (optional):  CN threshold value to classify genes within gained regions (default is "95th percentile" of all CN values)
#   clinical_info (optional): Location of xslx file with clinical information
#   subject_id (optional):    Subject ID required to match sample with clinical information (specified in flag --clinical_info)
#   sample_type (optional):   Type of investigated sample (e.g. fresh frozen tissue, organoid). This information is for annotation purposes only
#   project (optional):   Project name. This information is for annotation purposes only
#   top_genes:    The number of top ranked genes to be presented (default is "10")
#   dataset_name_incl:  Include dataset in the report name. Available options are: "TRUE" and "FALSE" (default)
#   save_tables:   Save interactive summary tables as HTML. Available options are: "TRUE" (default) and "FALSE"
#   hide_code_btn : Hide the "Code" button allowing to show/hide code chunks in the final HTML report. Available options are: "TRUE" (default) and "FALSE"
#   grch_version :  Human reference genome version used for genes annotation (default is "37")
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

#===============================================================================
#    Define functions
#===============================================================================

##### Create 'not in' operator
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option(c("-s", "--sample_name"), action="store", default=NA, type='character',
              help="Desired sample name to be presented in the report"),
  make_option(c("-o", "--dataset"), action="store", default=NA, type='character',
              help="Dataset to be used as external reference cohort"),
  make_option(c("-c", "--bcbio_rnaseq"), action="store", default=NA, type='character',
              help="Location of the results folder from bcbio RNA-seq pipeline"),
  make_option(c("-r", "--report_dir"), action="store", default=NA, type='character',
              help="Desired location for the report"),
  make_option(c("-v", "--ref_data_dir"), action="store", default="../data", type='character',
              help="Location of the reference and annotation files"),
  make_option(c("-t", "--transform"), action="store", default="CPM", type='character',
              help="Transformation method to be used when converting read counts"),
  make_option(c("-n", "--norm"), action="store", default=NA, type='character',
              help="Normalisation method"),
  make_option(c("-k", "--batch_rm"), action="store", default=TRUE, type='character',
              help="Remove batch-associated effects between datasets"),
  make_option(c("-f", "--filter"), action="store", default=TRUE, type='character',
              help="Filtering out low expressed genes"),
  make_option(c("-l", "--log"), action="store", default=TRUE, type='character',
              help="Log (base 2) transform data before normalisation"),
  make_option(c("-z", "--scaling"), action="store", default="gene-wise", type='character',
              help="Scaling for z-score transformation (gene-wise or group-wise"),
  make_option(c("-g", "--umccrise"), action="store", default=NULL, type='character',
              help="Location of the corresponding WGS-related data"),
  make_option(c("-j", "--pcgr_tier"), action="store", default=4, type='character',
              help="Tier threshold for reporting variants reported in PCGR"),
  make_option(c("-y", "--pcgr_splice_vars"), action="store", default=TRUE, type='character',
              help="Include non-coding splice region variants reported in PCGR"),
  make_option(c("-a", "--cn_loss"), action="store", default=5, type='character',
              help="CN threshold value to classify genes within lost regions"),
  make_option(c("-b", "--cn_gain"), action="store", default=95, type='character',
              help="CN threshold value to classify genes within gained regions"),
  make_option(c("-m", "--clinical_info"), action="store", default=NA, type='character',
              help="Location of xslx file with clinical information"),
  make_option(c("-i", "--subject_id"), action="store", default=NA, type='character',
              help="Subject ID"),
  make_option(c("-st", "--sample_type"), action="store", default=NA, type='character',
              help="Type of investigated sample"),
  make_option(c("-pr", "--project"), action="store", default=NA, type='character',
              help="Project name"),
  make_option(c("-x", "--top_genes"), action="store", default=10, type='character',
              help="The number of top ranked genes to be presented"),
  make_option(c("-w", "--dataset_name_incl"), action="store", default=NA, type='character',
              help="Include dataset in the report name"),
  make_option(c("-u", "--save_tables"), action="store", default=TRUE, type='character',
              help="Save interactive summary tables as HTML"),
  make_option(c("-d", "--hide_code_btn"), action="store", default=TRUE, type='character',
              help="Hide the \"Code\" button allowing to show/hide code chunks in the final HTML report"),
  make_option(c("-e", "--grch_version"), action="store", default=NA, type='character',
              help="human reference genome version used for genes annotation")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in argument from command line and check if all were provide by the user
if ( is.na(opt$sample_name) || is.na(opt$dataset) || is.na(opt$bcbio_rnaseq) || is.na(opt$report_dir) ) {

  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset PAAD  --bcbio_rnaseq $(pwd)/../data/test_data/final/test_sample_WTS  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report\n\n")
  q()
}

##### Make sure that sample ID is availabe if clincal data is provided
if ( !is.na(opt$clinical_info) && is.na(opt$subject_id)  ) {
  
  cat("\nSubject ID is missing! Please provide subject ID used in the clinical data by using \"--subject_id\" argument.\n\n")
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
  dataset_name_incl <- paste0(".", opt$dataset)
}

if ( is.na(opt$grch_version)  ) {
  ensembl_version <- 75
  ucsc_genome_assembly <- 19
  
} else if ( opt$grch_version == 37 ) {
  ensembl_version <- 75
  ucsc_genome_assembly <- 19
  
} else if ( opt$grch_version == 38 ) {
  ensembl_version <- 86
  ucsc_genome_assembly <- 38
  
} else {
  cat("\nCurrently human reference genome (GRCh) versions \"37\" and \"38\" are supported.\n\n")
  q()
}

##### Check if specified dataset type is valid
if ( toupper(opt$dataset) %!in% toupper(c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "BLCA-NET", "PAAD-IPMN", "PAAD-NET", "PAAD-ACC")) ) {
  
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

##### Pass the user-defined arguments to the RNAseq_report R markdown script and generate the report
rmarkdown::render(input = "RNAseq_report.Rmd", output_file = paste0(opt$sample_name, toupper(dataset_name_incl), ".RNAseq_report.html"), output_dir = opt$report_dir, params = list(sample_name = opt$sample_name, dataset = toupper(opt$dataset), bcbio_rnaseq = opt$bcbio_rnaseq, report_dir = opt$report_dir, ref_data_dir = opt$ref_data_dir, transform = opt$transform, norm = opt$norm, batch_rm = as.logical(opt$batch_rm), filter = as.logical(opt$filter), log = as.logical(opt$log), scaling = opt$scaling, umccrise = opt$umccrise, clinical_info = opt$clinical_info, subject_id = opt$subject_id, sample_type = opt$sample_type, project = opt$project, dataset_name_incl = dataset_name_incl, save_tables = as.logical(opt$save_tables), pcgr_tier = as.numeric(opt$pcgr_tier), pcgr_splice_vars = as.logical(opt$pcgr_splice_vars), cn_loss = as.numeric(opt$cn_loss), cn_gain = as.numeric(opt$cn_gain), top_genes = as.numeric(opt$top_genes), hide_code_btn = as.logical(opt$hide_code_btn), grch_version = as.numeric(opt$grch_version), ensembl_version = as.numeric(ensembl_version), ucsc_genome_assembly = as.numeric(ucsc_genome_assembly)))

##### Remove the assocaited MD file and the redundant folder with plots that are imbedded in the HTML report
unlink(paste0(opt$report_dir, "/", opt$sample_name, toupper(dataset_name_incl), ".RNAseq_report.md"), recursive = TRUE)
unlink(paste0(opt$report_dir, "/", opt$sample_name, toupper(dataset_name_incl), ".RNAseq_report_files"), recursive = TRUE)