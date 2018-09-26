################################################################################
#
#   File name: RNAseq_ref_cohorts_report.R
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
#	  Command line use example: Rscript  RNAseq_ref_cohorts_report.R  --sample_name CCR170012_MH17T001P013  --count_file ../data/CCR170012_MH17T001P013-ready.counts  --report_dir ../reports
#
#   sample_name:   Desired sample name to be presented in the report
#   count_file:    Location and name of the read count file from bcbio RNA-seq pipeline
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
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option(c("-s", "--sample_name"), action="store", default=NA, type='character',
              help="Desired sample name to be presented in the report"),
  make_option(c("-c", "--count_file"), action="store", default=NA, type='character',
              help="Location and name of the read count file from bcbio RNA-seq pipeline"),
  make_option(c("-r", "--report_dir"), action="store", default=NA, type='character',
              help="Desired location for the report")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in argument from command line and check if all were provide by the user
if ( is.na(opt$sample_name) || is.na(opt$count_file) || is.na(opt$report_dir) ) {
  
  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\n Rscript  RNAseq_ref_cohorts_report.R  --sample_name CCR170012_MH17T001P013  --count_file ../data/CCR170012_MH17T001P013-ready.counts  --report_dir ../reports\n\n")
  
  q()
}


##### Pass the user-defined arguments to the RNAseq_report R markdown script and generate the report
rmarkdown::render(input = "RNAseq_ref_cohorts_report.Rmd", output_file = paste0(opt$sample_name, ".RNAseq_ref_cohorts_report.html"), output_dir = opt$report_dir, params = list(report_dir = opt$report_dir, sample_name = opt$sample_name, count_file = opt$count_file))
