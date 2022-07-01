#' Combine sample expression profile with reference datasets
#'
#' Combines sample expression profile with reference datasets. Outputs a vector
#' with first element containing the merged data and second element containing 
#' merged targets info
#'
#' @param sample_name Sample name.
#' @param sample_counts Sample counts.
#' @param ref_data Reference data.
#' @param report_dir Report directory.
#' @param dataset Sample expression profile.
#'
#' @return Combined sample expression profile with reference datasets.
#' @export
combineDatasets <- function(sample_name, sample_counts, ref_data, report_dir, dataset) {

  ##### Extract info about target file for the external reference dataset
  target.ext <- utils::read.table(ref_data[["ext_ref"]][2], sep="\t", as.is=TRUE, header=TRUE)
  target.ext <- cbind(target.ext, rep(ref_data[["ext_ref"]][3], nrow(target.ext)))
  colnames(target.ext)[ncol(target.ext)] <- "Dataset"

  ##### Add prexit to sample names
  rownames(target.ext) <- paste(target.ext[,"Dataset"], target.ext[,"Sample_name"], sep = ".")
  target.ext <- target.ext[, -1]

  ##### Extract info about target file for the internal reference dataset
  target.int <- utils::read.table(ref_data[["int_ref"]][2], sep="\t", as.is=TRUE, header=TRUE)
  target.int <- cbind(target.int, rep(ref_data[["int_ref"]][3], nrow(target.int)))
  colnames(target.int)[ncol(target.int)] <- "Dataset"

  ##### Add prexit to sample names
  rownames(target.int) <- paste(target.int[,"Dataset"], target.int[,"Sample_name"], sep = ".")
  target.int <- target.int[, -1]

  target.comb <- rbind(target.ext, target.int)

  ##### Add sample info
  target.sample <- data.frame(sample_name, sample_name)
  names(target.sample) <- names(target.comb)
  rownames(target.sample) <- sample_name
  target.comb <- rbind( target.comb, target.sample )

  ##### Make syntactically valid names
  rownames(target.comb) <- make.names(rownames(target.comb))

  ##### Read sample read count file and combine it with reference datasets
  datasets.comb <- sample_counts
  names(datasets.comb) <- c("", sample_name)

  ##### list genes present in the sample read count file
  gene_list <- as.vector(datasets.comb[,1])

  ##### Loop through the expression data from different datasets and merge them into one matrix
  for ( i in 1:length(ref_data) ) {

    dataset.counts <- as.data.frame( utils::read.table(gzfile(ref_data[[i]][1]), header=TRUE, sep="\t", row.names=NULL) )

    ##### Add prexit to sample names
    colnames(dataset.counts) <- paste(unique(target.comb[,"Dataset"])[i], colnames(dataset.counts), sep = ".")

    ##### List genes present in individal files
    gene_list <- c( gene_list, as.vector(dataset.counts[,1]) )

    ##### Merge the expression datasets and make sure that the genes order is the same
    datasets.comb <- merge( datasets.comb, dataset.counts, by=1, all = FALSE, sort= TRUE)
  }

  ##### Use gene IDs as rownames
  rownames(datasets.comb) <- datasets.comb[,1]
  datasets.comb <- datasets.comb[, -1]

  ##### Make syntactically valid names
  colnames(datasets.comb) <- make.names(colnames(datasets.comb))

  ##### Make sure that the target file contains info only about samples present in the data matrix
  target.comb <- target.comb[ rownames(target.comb) %in% colnames(datasets.comb),  ]

  ##### Make sure that the samples order in the data matrix is the same as in the target file
  datasets.comb <- datasets.comb[ , rownames(target.comb) ]

  ##### Identify genes that were not present across all per-sampel files and were ommited in the merged matrix
  gene_list <- unique(gene_list)
  gene_list.missing <- gene_list[ gene_list %!in% rownames(datasets.comb) ]

  ##### Write list of missing genes into a file
  if ( length(gene_list.missing) > 0 ) {
    utils::write.table(prepare2write(gene_list.missing), file = paste0(report_dir, "/", sample_name, ".RNAseq_report.missing_genes.txt"), sep="\t", quote=FALSE, row.names=TRUE, append = FALSE )
  }

  ##### Clean the space and return output
  rm(sample_name, sample_counts, ref_data, target.ext, target.int, target.sample, dataset.counts, gene_list, gene_list.missing)
  return( list(datasets.comb, target.comb) )
}
