################################################################################
#
#   File name: recount2_normals_prep.R
#
#   Authors: Jacek Marzec (jacek.marzec@unimelb.edu.au)
#
#   University of Melbourne Centre for Cancer Research,
#   Victorian Comprehensive Cancer Centre
#   305 Grattan St, Melbourne, VIC 3000
#
################################################################################

################################################################################
#
#	  Description: Script preparing expression data from normal pancreas samples from Genotype-Tissue Expression (GTEx) study using recount package (197 samples from study SRP012682, https://jhubiostatistics.shinyapps.io/recount/_w_c6ab1f3e/#tab-9377-6 and https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP012682). The script outputs read count and target files in user-defined directory.
#
#	  Command line use example: Rscript  recount2_normals_prep.R  --r_data /data/rse_gene_pancreas.Rdata  --out_dir /data
#
#   r_data:        Full path with name of the Rdata file with expression data from normal pancreas samples from Genotype-Tissue Expression (GTEx) 
#   out_dir:       Desired location for the output files
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Functions
#===============================================================================

##### Prepare object to write into a file
prepare2write <- function (x) {
  
  x2write <- cbind(rownames(x), x)
  colnames(x2write) <- c("",colnames(x))
  return(x2write)
}

##### Filter duplicate gene by selecting gene with lowest standard deviation across samples
multiGeneFilter <- function(data) {
  
  cat(paste("Processing expression data for ", dim(data)[1], " genes...\n", sep=""))
  
  #####  Compute sample variance
  data <- cbind(rownames(data), data)
  exprs <- apply(data[,-1], 2,as.numeric)
  rownames(exprs) <- rownames(data)
  nLastCol <- ncol(exprs)
  sampleVar <- apply(exprs,1,var)
  
  nGenes <- nrow(exprs)
  
  genesUnique <- as.vector(unique(data[,1]))
  
  #####  To avoid NA issues]
  genesUnique <- as.vector(na.omit(genesUnique))
  nGenesUnique <- length(genesUnique)
  
  nSamples <- ncol(exprs)
  
  #####  Create matrix for filtered data
  exprs.filtered <- array(dim=c(nGenesUnique,nLastCol))
  rownames(exprs.filtered) <- genesUnique
  colnames(exprs.filtered) <- colnames(exprs)
  
  cat("Removing duplicate genes (selecting for lowest standard deviation)...\n")
  
  ##### Create progress bar
  pb <- txtProgressBar(min = 0, max = length(genesUnique), style = 3)
  
  for (i in 1:length(genesUnique)) {
    
    #####  Index/indices of all genes matching gene (detect duplicates)
    gene <- genesUnique[i]
    idxGenes <- seq(along=1:nGenes)[data[,1]==gene]
    
    exprs.slice <- exprs[idxGenes,]
    
    #####  Find duplicates with lowest variance
    if (length(idxGenes)>1) {
      idxMinVar <- which.min(sampleVar[idxGenes])
      
      exprs.slice <- exprs.slice[idxMinVar,]
    }
    exprs.filtered[rownames(exprs.filtered)==gene,] <- as.matrix(exprs.slice)
    
    ##### Report progress
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cat(paste(nGenesUnique," unique genes names remain.\n",sep=""))
  
  rm(data, exprs, exprs.slice)
  
  return(exprs.filtered)
}

#===============================================================================
#    Load libraries
#===============================================================================

suppressMessages(library(optparse))
suppressMessages(library(recount))

#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option(c("-r", "--r_data"), action="store", default=NA, type='character',
              help="Full path with name of the Rdata file from Genotype-Tissue Expression (GTEx)"),
  make_option(c("-o", "--out_dir"), action="store", default=NA, type='character',
              help="Desired location for the output files")
)

opt = parse_args(OptionParser(option_list=option_list))

##### Read in argument from command line and check if all were provide by the user
if ( is.na(opt$r_data) || is.na(opt$out_dir) ) {
  
  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript  recount2_normals_prep.R  --r_data /data/rse_gene_pancreas.Rdata  --out_dir /data\n\n")
  
  q()
}

#===============================================================================
#    Load and prepare data
#===============================================================================

##### Load the RangedSummarizedExperiment (RSE) object with expressin data at the gene level for GTEx study (SRP012682, https://jhubiostatistics.shinyapps.io/recount/_w_c6ab1f3e/#tab-9377-6) limited to pancreas tissue
##### MOre info about RSE: https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
cat("\nLoading recount data...\n\n")
load(opt$r_data)

##### Extract phenotype data provided by the recount project
rse_GTEx.pheno <- colData(rse_gene)

##### Compute read counts
rse_GTEx <- read_counts(rse_gene, round = FALSE)

##### Extract read counts
rse_GTEx_counts <- assays(rse_GTEx)$counts

##### Change the Ensembl version to stable Ensembl gene id
rownames(rse_GTEx_counts) <- gsub('\\..*', '', rownames(rse_GTEx_counts))

##### Filter duplicated genes by selecting gene with lowest standard deviation across samples
cat("\nFiltering duplicated genes by selecting gene with lowest standard deviation across samples...\n\nNote that this step may take long time!\n\n")
rse_GTEx_counts  <- multiGeneFilter(rse_GTEx_counts)

#===============================================================================
#    Write data into files
#===============================================================================

##### Write expression data into a file
cat( paste( "\nWriting read count data from normal pancreas samples [ GTEx_normal_pancreas_Counts.exp ] and associated target file [ GTEx_normal_pancreas_Target.txt ] to [", opt$out_dir, "]...\n\n", sep=" ") )
write.table( prepare2write(rse_GTEx_counts), file = paste(opt$out_dir, "/GTEx_normal_pancreas_Counts.exp", sep=""), quote=FALSE ,sep="\t", row.names=FALSE )

##### Create target file
rse_GTEx.pheno.df <- as.data.frame( cbind(rse_GTEx.pheno$run, rep("rse_gene_pancreas.Rdata", nrow(rse_GTEx.pheno)), rep("Pancreas (normal)", nrow(rse_GTEx.pheno)), rep("", nrow(rse_GTEx.pheno))) )
names(rse_GTEx.pheno.df) <- c("Sample_name", "File_name", "Target", "Replicates")

write.table( rse_GTEx.pheno.df, file = paste(opt$out_dir, "/GTEx_normal_pancreas_Target.txt", sep=""), quote=FALSE ,sep="\t", row.names=FALSE )

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
