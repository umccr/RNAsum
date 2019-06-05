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
#	  Command line use example: Rscript RNAseq_report.R  --sample_name CCR170115b_MH17T002P033_RNA  --tissue pancreas  --count_file ../data/CCR170115b_MH17T002P033_RNA-ready.counts  --report_dir ../RNAseq_report  --transform CPM  --norm TMM  --filter TRUE  --log TRUE  --sample_id 2016.249.17.MH.P033  --batch ../data//umccrised/2016_249_17_MH_P033__CCR170115b_MH17T002P033  --clinical_info ../data/clinical_data.xlsx  --plots_mode static
#
#   sample_name:   Desired sample name to be presented in the report
#   tissue:        Tissue from which the samples were derived
#   count_file:    Location and name of the read count file from bcbio RNA-seq pipeline
#   report_dir:    Desired location for the report
#   transform:    Transformation method to be used when converting read counts. Available options are: "CPM" (defualt) and "TPM"
#   norm:         Normalisation method. Currently, "TMM" is used for CPM-transformed data and "quantile" normalisation is used for TPM-transformed data
#   filter:       Filtering out low expressed genes. Available options are: "TRUE" (defualt) and "FALSE"
#   log:          Log (base 2) transform data before normalisation. Available options are: "TRUE" (defualt) and "FALSE"
#   scaling:      Apply row-wise (across samples) or column-wise (across genes in a sample) data scaling. Available options are: "sample-wise" (across samples, default) or "gene-wise" (across genes)
#   sample_id:     Sample ID
#   batch (optional):   Location of the corresponding WGS-related data (with PURPLE and Manta output files)
#   clinical_info (optional):   Location of xslx file with clinical information
#   plots_mode:    Plotting mode. Available options are: "Static" (default), "interactive" and "semi-interactive"
#   hide_code_btn :    Hide the "Code" button allowing to show/hide code chunks in the final HTML report. Available options are: "TRUE" (defualt) and "FALSE"
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
  make_option(c("-o", "--tissue"), action="store", default=NA, type='character',
              help="Tissue from which the samples were derived"),
  make_option(c("-c", "--count_file"), action="store", default=NA, type='character',
              help="Location and name of the read count file from bcbio RNA-seq pipeline"),
  make_option(c("-r", "--report_dir"), action="store", default=NA, type='character',
              help="Desired location for the report"),
  make_option(c("-t", "--transform"), action="store", default=NA, type='character',
              help="Transformation method to be used when converting read counts"),
  make_option(c("-n", "--norm"), action="store", default=NA, type='character',
              help="Normalisation method"),
  make_option(c("-f", "--filter"), action="store", default=NA, type='character',
              help="Filtering out low expressed genes"),
  make_option(c("-l", "--log"), action="store", default=NA, type='character',
              help="Log (base 2) transform data before normalisation"),
  make_option(c("-z", "--scaling"), action="store", default=NA, type='character',
              help="Scaling for z-score transformation (sample-wise or gene-wise"),
  make_option(c("-i", "--sample_id"), action="store", default=NA, type='character',
              help="Sample ID"),
  make_option(c("-b", "--batch"), action="store", default=NA, type='character',
              help="Location of the corresponding WGS-related data"),
  make_option(c("-m", "--clinical_info"), action="store", default=NA, type='character',
              help="Location of xslx file with clinical information"),
  make_option(c("-p", "--plots_mode"), action="store", default=NA, type='character',
              help="Static (default), interactive or semi-interactive mode for plots"),
  make_option(c("-d", "--hide_code_btn"), action="store", default=NA, type='character',
              help="Hide the \"Code\" button allowing to show/hide code chunks in the final HTML report")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in argument from command line and check if all were provide by the user
if ( is.na(opt$sample_name) || is.na(opt$tissue) || is.na(opt$count_file) || is.na(opt$report_dir) ) {

  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript RNAseq_report.R  --sample_name CCR170115b_MH17T002P033_RNA  --tissue pancreas  --count_file ../data/CCR170115b_MH17T002P033_RNA-ready.counts  --report_dir ../RNAseq_report\n\n")

  q()
}

##### Make sure that sample ID is availabe if clincal data is provided
if ( !is.na(opt$clinical_info) && is.na(opt$sample_id)  ) {
  
  cat("\nSample ID is missing! Please provide sample ID used in the clinical data.\n\n")
  
  q()
}

##### Set default parameters
if ( is.na(opt$transform)  ) {
  
  opt$transform <- "CPM"
}

if ( is.na(opt$norm)  ) {
  
  opt$norm <- "TMM"
}

if ( is.na(opt$filter)  ) {
  
  opt$filter <- TRUE
}

if ( is.na(opt$log)  ) {
  
  opt$log <- TRUE
}

if ( is.na(opt$scaling)  ) {
  
  opt$scaling <- "sample-wise"
}

if ( is.na(opt$plots_mode)  ) {
  
  opt$plots_mode <- "static"
}

if ( is.na(opt$hide_code_btn)  ) {
  
  opt$hide_code_btn <- TRUE
}

##### Check if specified tissue type is valid
if ( opt$tissue %!in% c("pancreas", "cervix") ) {
  
  cat("\nInvalid tissue type! Please use one of the following:\n\n")
  cat("[ pancreas ] - this will compare the patient's data in the context of samples from TCGA pancreatic ductal adenocarcinoma (PAAD) cohort and GTEx samples from healthy pancreas tissue\n")
  cat("[ cervix ] - this will compare the patient's data in the context of samples from TCGA cervical squamous cell carcinoma (CESC) cohort and GTEx samples from healthys cervix uteri tissue\n\n")
  q()
}

##### Embed static (defulat) plots (to reduce the report size), unless interactive or semi-interactive mode is specified
if ( is.na(opt$plots_mode) ) {
  
  opt$plots_mode<- "static"
  
} else if ( opt$plots_mode != "interactive" && opt$plots_mode != "semi-interactive" && opt$plots_mode != "static" ) {
  
  cat("\nThe plots mode \"", opt$plots_mode, "\" is invalid! Please select either \"static\", \"interactive\" or \"semi-interactive\" mode.\n\n")
  q()
}

##### Make sure that TMM, TMMwzp, RLE or upperquartile normalisation is used for CPM-tansformed data and quantile normalisation is used for TPM-tansformed data
if ( opt$transform == "TPM" && opt$norm == "TMM" ) {
  
  cat(paste0("\nOnly TPM normalisation is not available for TPM-tansformed data!\n\nQuantile normalisation will be performed for ", opt$transform, "-tansformed data.\n\n"))
  
  opt$norm <- "quantile"
  
} else if ( opt$transform == "CPM" && opt$norm == "quantile" ) {
  
  cat(paste0("\nQuantile normalisation is available only for TPM-tansformed data! \"TMM\", \"TMMwzp\", \"RLE\" and \"upperquartile\" methods are available for ", opt$transform, "-tansformed data.\n\n"))
  
  q()
  
} else if ( opt$transform == "CPM" &&  opt$norm != "TMM" && opt$norm != "TMMwzp" && opt$norm != "RLE" && opt$norm != "upperquartile" ) {
  
  cat(paste0("\nWrong normalisation method was selected! \"TMM\", \"TMMwzp\", \"RLE\" and \"upperquartile\" methods are available for ", opt$transform, "-tansformed data.\n\n"))
  
  q()
}

##### Pass the user-defined arguments to the RNAseq_report R markdown script and generate the report
rmarkdown::render(input = "RNAseq_report.Rmd", output_file = paste0(opt$sample_name, ".RNAseq_report.html"), output_dir = opt$report_dir, params = list(sample_name = opt$sample_name, tissue = tolower(opt$tissue), count_file = opt$count_file, report_dir = opt$report_dir, transform = opt$transform, norm = opt$norm, filter = as.logical(opt$filter), log = as.logical(opt$log), scaling = opt$scaling, sample_id = opt$sample_id, batch = opt$batch, clinical_info = opt$clinical_info, plots_mode = tolower(opt$plots_mode), hide_code_btn = as.logical(opt$hide_code_btn)))
