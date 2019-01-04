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
#	  Command line use example: Rscript RNAseq_report.R  --sample_name CCR170012_MH17T001P013  --tissue pancreas  --count_file ../data/CCR170012_MH17T001P013-ready.counts  --plots_mode static  --report_dir ../reports
#
#   sample_name:   Desired sample name to be presented in the report
#   tissue:        Tissue from which the samples were derived
#   count_file:    Location and name of the read count file from bcbio RNA-seq pipeline
#   plots_mode:    Static (default) or interactive mode for plots
#   report_dir:    Desired location for the report
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
  make_option(c("-t", "--tissue"), action="store", default=NA, type='character',
              help="Tissue from which the samples were derived"),
  make_option(c("-c", "--count_file"), action="store", default=NA, type='character',
              help="Location and name of the read count file from bcbio RNA-seq pipeline"),
  make_option(c("-p", "--plots_mode"), action="store", default=NA, type='character',
              help="Static (default) or interactive mode for plots"),
  make_option(c("-r", "--report_dir"), action="store", default=NA, type='character',
              help="Desired location for the report")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in argument from command line and check if all were provide by the user
if ( is.na(opt$sample_name) || is.na(opt$tissue) || is.na(opt$count_file) || is.na(opt$report_dir) ) {

  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript RNAseq_report.R  --sample_name CCR170012_MH17T001P013  --tissue pancreas  --count_file ../data/CCR170012_MH17T001P013-ready.counts  --plots_mode static  --report_dir ../reports\n\n")

  q()
}

##### Check if specified tissue type is valid
if ( opt$tissue %!in% c("pancreas", "cervix") ) {
  
  cat("\nInvalid tissue type! Please use one of the following:\n\n")
  cat("[ pancreas ] - this will compare the patient's data in the context of samples from TCGA pancreatic ductal adenocarcinoma (PAAD) cohort and GTEx samples from healthy pancreas tissue\n")
  cat("[ cervix ] - this will compare the patient's data in the context of samples from TCGA cervical squamous cell carcinoma (CESC) cohort and GTEx samples from healthys cervix uteri tissue\n\n")
  q()
}

##### Embed static (defulat) plots (to reduce the report size), unless interactive mode is specified
if ( is.na(opt$plots_mode) ) {
  
  opt$plots_mode<- "static"
  
} else if ( opt$plots_mode != "interactive" && opt$plots_mode != "static" ) {
  
  cat("\nThe plots mode \"", opt$plots_mode, "\" is invalid! Please select either \"static\" or \"interactive\" mode.\n\n")
  q()
}

##### Pass the user-defined arguments to the RNAseq_report R markdown script and generate the report
rmarkdown::render(input = "RNAseq_report.Rmd", output_file = paste0(opt$sample_name, ".RNAseq_report.html"), output_dir = opt$report_dir, params = list(report_dir = opt$report_dir, sample_name = opt$sample_name, tissue = opt$tissue, plots_mode = opt$plots_mode, count_file = opt$count_file))
