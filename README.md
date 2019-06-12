# RNA-seq report

This reporting tool is designed to post-process, summarise and visualise an output from *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)*. Its main application is to complement genome-based findings and to provide additional evidence for detected alterations.

NOTE, currently the pipeline is limited to report on samples from **pancreatic tissue** since only [pancreatic cancer reference cohort](https://github.com/umccr/Pancreatic-data-harmonization/tree/master/expression/in-house) have been assembled to date (see [Reference data](#reference-data) section). An exception is the *Fusion genes* report section, which relies solely on the whole-transcriptome sequencing (WTS) data from analysed sample.

## Table of contents

<!-- vim-markdown-toc GFM -->
* [Installation](#installation)
* [Workflow](#workflow)
* [Reference data](#reference-data)
  * [Internal reference cohort](#internal-reference-cohort)
  * [External reference cohort](#external-reference-cohort)
* [Testing](#testing)
* [UMCCR HPC](#umccr-hpc)
* [Usage](#usage)
  * [Arguments](#arguments)
  * [Command line use example](#command-line-use-example)


<!-- vim-markdown-toc -->


## Installation

To be done...


## Workflow

The pipeline consist of four main components illustrated and breifly described below. See the [workflow.md](workflow.md) for the complete description of the workflow.

<img src="img/RNAseq_report_workflow.png" width="100%"> 

<br/>

1. Process the RNA-seq data from *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* including **gene fusions** and per-gene **read counts**. The [gene fusion](./fusions)  candidates are re-quantified and the read counts are [normalised, transformed](img/counts_post-processing_scheme.pdf) and [converted](img/Z-score_transformation_sample_wise.pdf) into a scale that allows to present the sample's expression measurements in the context of the [reference cohorts](#reference-data).

2. Feed in **genome-based findings** from whole-genome sequencing (WGS) data to focus on genes of interest and provide additional evidence for dysregulation of mutated genes, or genes located within detected structural variants (SVs) or copy-number (CN) altered regions. The RNA-seq report pipeline is designed to be compatible with WGS patient report based on [umccrise](https://github.com/umccr/umccrise) pipeline output.


3. Collate results with knowledge derived from public **databases** and **reference cohorts** to provide additional source of evidence, e.g. about variants clinical significance or potential druggable targets, for the genes of interest and to get an idea about their expression levels in other cancer [patient cohorts](#reference-data).

4. The final product is a html-based **interactive report** with searchable tables and plots presenting expression levels of the genes of interest genes. The report consist of several sections described in [report_structure.md](report_structure.md).


## Reference data

Currently, the reference data is availbale only for **pancreatic** cancer samples.


### Internal reference cohort

In order to explore expression changes in queried sample we have built a high-quality pancreatic cancer reference cohort. This set of **40 samples** is based on WTS data generated at **UMCCR** and processed with **[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)** *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* to minimise potential batch effects between queried sample and the reference cohort and to make sure the data are comparable. The internal reference cohort assembly is summarised in [Pancreatic-data-harmonization](https://github.com/umccr/Pancreatic-data-harmonization/tree/master/expression/in-house) repository.

### External reference cohort

Additionally, pancreas adenocarcinoma (PAAD) expression data (**150 samples**) from [TCGA](https://tcga-data.nci.nih.gov/) have been processed and used as an external reference cohort (see [Pancreatic-data-harmonization](https://github.com/umccr/Pancreatic-data-harmonization/blob/master/expression/public/README.md#tcga-paad) repository for more details). This dataset, however, is expected to demonstrate prominent batch effects when compared to UMCCR WTS data due to differences in applied experimental procedures and analytical pipelines. Moreover, TCGA data may include samples from tissue material of lower quality and cellularity compared to samples processed at UMCCR.


## Testing

To be done...


## UMCCR HPC

Currently, we the pipeline is run on local machines but the aim is to set it up on [NCI Raijin](https://github.com/umccr/wiki/blob/master/computing/clusters/raijin-intro.md) and eventually on [Amazon Web Services](https://github.com/umccr/wiki/blob/master/computing/cloud/aws.md) (AWS) cloud computing.

To be done...


## Usage

To run the pipeline execure the *[RNAseq_report.R](./rmd_files/RNAseq_report.R)* script. This script catches the arguments from the command line and passes them to the *[RNAseq_report.Rmd](./rmd_files/RNAseq_report.Rmd)* script to produce the interactive HTML report.

### Arguments

Argument | Description | Required
------------ | ------------ | ------------
--sample_name | Desired sample name to be presented in the report | **Yes**
--tissue | Tissue from which the sample is derived | **Yes**
--count_file | Location and name of the read count file from *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* | **Yes**
--report_dir | Desired location for the report | **Yes**
--transform | Transformation method for converting read counts. Available options are: `CPM` (defualt) and `TPM` | No
--norm | Normalisation method. Available options are: `TMM` (for *CPM-transformed* data) and `quantile` (for *TPM-transformed* data) | No
--filter | Filtering out low expressed genes. Available options are: `TRUE` (defualt) and `FALSE` | No
--log | Log (base 2) transform data before normalisation. Available options are: `TRUE` (defualt) and `FALSE` | No
--scaling | Apply row-wise (across samples) or column-wise (across genes in a sample) data scaling. Available options are: `sample-wise` (across samples, default) or `gene-wise` (across genes) | No
--umccrise | Location of the corresponding *[umccrise](https://github.com/umccr/umccrise)* output (including [PCGR](https://github.com/sigven/pcgr), [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) and [Manta](https://github.com/Illumina/manta) output files) from genomic-related data | No
--pcgr_tier | Tier threshold for reporting variants reported in [PCGR](https://github.com/sigven/pcgr) (if [PCGR](https://github.com/sigven/pcgr) results are available, default is `3`) | No
--cn_loss | CN threshold value to classify genes within lost regions (if CN results from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) are available, default is `1.5`) | No
--cn_gain | CN threshold value to classify genes within gained regions (if CN results from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) are available, default is `3`) | No
--clinical_info | Location of *xslx* file with clinical information | No
--subject_id | Subject ID required to match sample with clinical information (if available) | No
--plots_mode | Plotting mode. Available options: `Static` (default), `interactive` and `semi-interactive` | No
--hide_code_btn | Hide the *Code* button allowing to show/hide code chunks in the final HTML report. Available options are: `TRUE` (defualt) and `FALSE` | No
--ensembl_version | Version of Ensembl database to be used for genes annotation (default is `86`) | No
--ucsc_genome_assembly | Version of UCSC Homo sapiens genome to be used for genes (default is `19`) | No

<br />

**Packages**: *[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)*, *[rapportools](https://cran.r-project.org/web/packages/rapportools/index.html)*, 
*[optparse](https://cran.r-project.org/web/packages/optparse/index.html)*, *[openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html)*, *[readr](https://cran.r-project.org/web/packages/readr/index.html)*, *[tidyverse](https://www.tidyverse.org/)*, *[dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)*, *[tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)*, *[rlang](https://cran.r-project.org/web/packages/rlang/index.html)*, *[DT](https://cran.r-project.org/web/packages/DT/index.html)*, *[kableExtra](https://cran.r-project.org/web/packages/kableExtra/index.html)*, *[matrixStats](https://cran.rstudio.com/web/packages/matrixStats/index.html)*, *[DataCombine](https://cran.r-project.org/web/packages/DataCombine/index.html)*, *[knitr](https://cran.r-project.org/web/packages/knitr/index.html)*, *[plotly](https://plot.ly/r/)*, *[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)*, *[ggforce](https://cran.r-project.org/web/packages/ggforce/index.html)*, *[magick](https://cran.r-project.org/web/packages/magick/index.html)*, *[lares](https://www.rdocumentation.org/packages/lares/versions/4.4)*, *[htmltools](https://cran.r-project.org/web/packages/htmltools/index.html)*, *[htmlwidgets](https://cran.r-project.org/web/packages/htmlwidgets/index.html)*, *[devtools](https://cran.r-project.org/web/packages/devtools/index.html)*, *[EnsDb.Hsapiens.v86](http://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v86.html)*\*, *[BSgenome.Hsapiens.UCSC.hg19](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html)*\**

\*  Ensembl based annotation version 86 is used as default but can be updated <br >
\**  UCSC Homo sapiens genome sequences version hg19 is used as default but can be updated

### Command line use example

```
Rscript RNAseq_report.R  --sample_name CCR170115b_MH17T002P033_RNA  --tissue pancreas  --count_file ../data/CCR170115b_MH17T002P033_RNA-ready.counts  --report_dir ../RNAseq_report  --transform CPM  --norm TMM  --filter TRUE  --log TRUE  --subject_id 2016.249.17.MH.P033  --umccrise ../data/umccrised/2016_249_17_MH_P033__CCR170115b_MH17T002P033  --clinical_info ../data/clinical_data.xlsx  --plots_mode semi-interactive
```

The interactive HTML report named `CCR170115b_MH17T002P033_RNA.RNAseq_report.html` will be created in `../RNAseq_report` folder.

