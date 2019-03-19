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
#	  Command line use example: Rscript RNAseq_report.R  --sample_name CCR170012_MH17T001P013  --tissue pancreas  --count_file ../data/CCR180038_SV18T002P006_RNA-ready.counts  --plots_mode static  --report_dir ../reports  --batch ../data/2016_249_18_SV_P006_1__CCR180038_SV18T002P006  --transform CPM  --norm TMM  --filter TRUE  --log TRUE
#
#   sample_name:   Desired sample name to be presented in the report
#   tissue:        Tissue from which the samples were derived
#   count_file:    Location and name of the read count file from bcbio RNA-seq pipeline
#   plots_mode:    Static (default), interactive, or semi-interactive mode for plots
#   report_dir:    Desired location for the report
#   batch (optional):   Location of the corresponding WGS-related data (with PURPLE and Manta output files)
#   transform:    Transformation method to be used when converting read counts. Available options are: "CPM" (defualt) and "TPM"
#   norm:         Normalisation method. Currently, "TMM" is used for CPM-transformed data and "quantile" normalisation is used for TPM-transformed data
#   filter:       Filtering out low expressed genes. Available options are: "TRUE" (defualt) and "FALSE"
#   log:          Log (base 2) transform data before normalisation. Available options are: "TRUE" (defualt) and "FALSE"
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
  make_option(c("-p", "--plots_mode"), action="store", default=NA, type='character',
              help="Static (default), interactive or semi-interactive mode for plots"),
  make_option(c("-r", "--report_dir"), action="store", default=NA, type='character',
              help="Desired location for the report"),
  make_option(c("-b", "--batch"), action="store", default=NA, type='character',
              help="Location of the corresponding WGS-related data"),
  make_option(c("-t", "--transform"), action="store", default=NA, type='character',
              help="Transformation method to be used when converting read counts"),
  make_option(c("-f", "--filter"), action="store", default=NA, type='character',
              help="Normalisation method"),
  make_option(c("-n", "--norm"), action="store", default=NA, type='character',
              help="Filtering out low expressed genes"),
  make_option(c("-l", "--log"), action="store", default=NA, type='character',
              help="Log (base 2) transform data before normalisation")
)

opt = parse_args(OptionParser(option_list=option_list))

opt$filter <- as.logical(opt$filter)
opt$log <- as.logical(opt$log)

##### Read in argument from command line and check if all were provide by the user
if ( is.na(opt$sample_name) || is.na(opt$tissue) || is.na(opt$count_file) || is.na(opt$report_dir) ) {

  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript RNAseq_report.R  --sample_name CCR170012_MH17T001P013  --tissue pancreas  --count_file ../data/CCR180038_SV18T002P006_RNA-ready.counts  --plots_mode static  --report_dir ../reports  --batch ../data/2016_249_18_SV_P006_1__CCR180038_SV18T002P006  --transform CPM  --norm TMM  --filter TRUE  --log TRUE\n\n")

  q()
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
rmarkdown::render(input = "RNAseq_report.Rmd", output_file = paste0(opt$sample_name, ".RNAseq_report.html"), output_dir = opt$report_dir, params = list(report_dir = opt$report_dir, sample_name = opt$sample_name, tissue = opt$tissue, plots_mode = opt$plots_mode, count_file = opt$count_file, batch = opt$batch, transform = opt$transform, filter = opt$filter, norm = opt$norm, log = opt$log))
