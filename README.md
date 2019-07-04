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
* [Usage](#usage)
  * [Arguments](#arguments)
  * [Examples](#examples)
  	 * [Required arguments only](#required-arguments-only)
  	 * [Add genome-based results](#add-genome-based-results)
  	 * [Add clinical information](#add-clinical-information)
  * [Output](#output)

<!-- vim-markdown-toc -->


## Installation

Run the [environment.yaml](envm/environment.yaml) file to create *conda* environment and install required packages. The `-p` flag should point to the *miniconda* installation path. For instance, to create `rnaseq-report` environment using *miniconda* installed in `/miniconda` directory run the following command:

```
conda env create -p /miniconda/envs/rnaseq-report --file envm/environment.yaml
```

Activate created `rnaseq-report` *conda* environment before running the pipeline

```
conda activate rnaseq-report
```


## Workflow

The pipeline consist of four main components illustrated and breifly described below. See the [workflow.md](workflow.md) for the complete description of the workflow.

<img src="img/RNAseq_report_workflow.png" width="100%"> 

<br/>

1. Process the RNA-seq data from *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* including **gene fusions** and per-gene **read counts**. The [gene fusion](./fusions)  candidates are re-quantified and the read counts are [normalised, transformed](img/counts_post-processing_scheme.pdf) and [converted](img/Z-score_transformation_gene_wise.pdf) into a scale that allows to present the sample's expression measurements in the context of the [reference cohorts](#reference-data).

2. Feed in **genome-based findings** from whole-genome sequencing (WGS) data to focus on genes of interest and provide additional evidence for dysregulation of mutated genes, or genes located within detected structural variants (SVs) or copy-number (CN) altered regions. The RNA-seq report pipeline is designed to be compatible with WGS patient report based on [umccrise](https://github.com/umccr/umccrise) pipeline output.


3. Collate results with knowledge derived from public **databases** and **reference cohorts** to provide additional source of evidence, e.g. about variants clinical significance or potential druggable targets, for the genes of interest and to get an idea about their expression levels in other cancer [patient cohorts](#reference-data).

4. The final product is a html-based **interactive report** with searchable tables and plots presenting expression levels of the genes of interest genes. The report consist of several sections described in [report_structure.md](report_structure.md).


## Reference data

Currently, the reference data is availbale only for **pancreatic** cancer samples.


### Internal reference cohort

In order to explore expression changes in queried sample we have built a high-quality pancreatic cancer reference cohort. This set of **40 samples** is based on WTS data generated at **UMCCR** and processed with **[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)** *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* to minimise potential batch effects between queried sample and the reference cohort and to make sure the data are comparable. The internal reference cohort assembly is summarised in [Pancreatic-data-harmonization](https://github.com/umccr/Pancreatic-data-harmonization/tree/master/expression/in-house) repository.

### External reference cohort

Additionally, pancreas adenocarcinoma (PAAD) expression data (**150 samples**) from [TCGA](https://tcga-data.nci.nih.gov/) have been processed and used as an external reference cohort (see [Pancreatic-data-harmonization](https://github.com/umccr/Pancreatic-data-harmonization/blob/master/expression/public/README.md#tcga-paad) repository for more details). This dataset, however, is expected to demonstrate prominent batch effects when compared to UMCCR WTS data due to differences in applied experimental procedures and analytical pipelines. Moreover, TCGA data may include samples from tissue material of lower quality and cellularity compared to samples processed at UMCCR.


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
--norm | Normalisation method. Available options are: `TMM` (for *CPM-transformed* data, defualt) and `quantile` (for *TPM-transformed* data) | No
--filter | Filtering out low expressed genes. Available options are: `TRUE` (defualt) and `FALSE` | No
--log | Log (base 2) transform data before normalisation. Available options are: `TRUE` (defualt) and `FALSE` | No
--scaling | Apply [`gene-wise`](img/Z-score_transformation_gene_wise.pdf) (default) or [`group-wise`](img/Z-score_transformation_group_wise.pdf) data scaling | No
--umccrise | Location of the corresponding *[umccrise](https://github.com/umccr/umccrise)* output (including [PCGR](https://github.com/sigven/pcgr), [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) and [Manta](https://github.com/Illumina/manta) output files) from genome-based data | No
--pcgr_tier | Tier threshold for reporting variants reported in [PCGR](https://github.com/sigven/pcgr) (if [PCGR](https://github.com/sigven/pcgr) results are available, default is `3`) | No
--cn_loss | CN threshold value to classify genes within lost regions (if CN results from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) are available, default is `1.5`) | No
--cn_gain | CN threshold value to classify genes within gained regions (if CN results from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) are available, default is `3`) | No
--clinical_info | Location of *xslx* file with clinical information | No
--subject_id | Subject ID required to match sample with clinical information (if available) | No
--plots_mode | Plotting mode. Available options: `interactive` (all possible plots will be interactive; default) and `semi-interactive` (only plots in `Input data`, `CN altered genes`, `Immune markers` and `HRD genes` sections will be interactive) or `static` (all plots will be static) | No
--hide_code_btn | Hide the *Code* button allowing to show/hide code chunks in the final HTML report. Available options are: `TRUE` (defualt) and `FALSE` | No
--grch_version | Human reference genome version used for genes annotation. Available options: `37` (default) and `38` | No

<br />

**Packages**: *[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)*, *[rapportools](https://cran.r-project.org/web/packages/rapportools/index.html)*, 
*[optparse](https://cran.r-project.org/web/packages/optparse/index.html)*, *[openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html)*, *[readr](https://cran.r-project.org/web/packages/readr/index.html)*, *[tidyverse](https://www.tidyverse.org/)*, *[dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)*, *[tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)*, *[rlang](https://cran.r-project.org/web/packages/rlang/index.html)*, *[DT](https://cran.r-project.org/web/packages/DT/index.html)*, *[kableExtra](https://cran.r-project.org/web/packages/kableExtra/index.html)*, *[matrixStats](https://cran.rstudio.com/web/packages/matrixStats/index.html)*, *[tibble](https://cran.r-project.org/web/packages/tibble/index.html)*, *[knitr](https://cran.r-project.org/web/packages/knitr/index.html)*, *[plotly](https://plot.ly/r/)*, *[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)*, *[ggforce](https://cran.r-project.org/web/packages/ggforce/index.html)*, *[pdftools](https://cran.r-project.org/web/packages/pdftools/index.html)*, *[png](https://cran.r-project.org/web/packages/png/index.html)*, *[lares](https://www.rdocumentation.org/packages/lares/versions/4.4)*, *[htmltools](https://cran.r-project.org/web/packages/htmltools/index.html)*, *[htmlwidgets](https://cran.r-project.org/web/packages/htmlwidgets/index.html)*, *[devtools](https://cran.r-project.org/web/packages/devtools/index.html)*, *[EnsDb.Hsapiens.v75](http://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v75.html)* (*[EnsDb.Hsapiens.v86](http://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v86.html)*)\*, *[BSgenome.Hsapiens.UCSC.hg19](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html)* (*[BSgenome.Hsapiens.UCSC.hg38](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html)*)\*

\*  Human reference genome ***[GRCh37](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)*** (*Ensembl* based annotation version ***75***) is used for genes annotation as default. Alternatively, human reference genome [GRCh38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39) (*Ensembl* based annotation version *86*) is used when argument `grch_version` is set to `38`.


### Examples 

Below are command line use examples for generating *Transcriptome Patient Summary* report using:

1. [required arguments only](#required-arguments-only)
2. **[genome-based results](#add-genome-based-results)**
3. [clinical information](#add-clinical-information)

**Note**:

* make sure that the created *conda* environment (see [Installation](#installation) section) is  activated

```
conda activate rnaseq-report
```

* *[RNAseq_report.R](./rmd_files/RNAseq_report.R)* script (see the beginning of [Usage](#usage) section) should be executed from [rmd_files](./rmd_files) folder

```
cd rmd_files
```

* example data is provided in [data/test_data](./data/test_data) folder


#### 1. Required arguments only

In this scenario, only expression levels of key **[`Cancer genes`](https://github.com/umccr/umccrise/blob/master/workflow.md#key-cancer-genes)**, **`Fusion genes`**, **`Immune markers`** and homologous recombination deficiency genes (**`HRD genes`**) will be reported. The genome-based findings will not be incorporated into the report, thus **no results will be provided in** ~~`Mutated genes`~~, ~~`Structural variants`~~ and ~~`CN altered genes`~~ sections. Moreover, gene fusions reported in `Fusion genes` section will not contain inforamation about evidence from genome-based data.

```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset pdac  --count_file $(pwd)/../data/test_data/final/test_sample_WTS/test_sample_WTS-ready.counts  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report
```

>The interactive HTML report named `test_sample_WTS.pdac.RNAseq_report.html` will be created in `data/test_data/final/test_sample_WTS/RNAseq_report` folder.


#### 2. Add genome-based results

This is the **preferred scenario for using** ***Transcriptome Patient Summary***, in which the genome-based will be primarily used for exploring expression levels of altered genes. The genome-based findings can be incorporated into the report by specifying location of the corresponding ***[umccrise](https://github.com/umccr/umccrise)*** output files (including results from [PCGR](https://github.com/sigven/pcgr), [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) and [Manta](https://github.com/Illumina/manta)) using `--umccrise` argument. In this scenario, **`Mutated genes`**, **`Structural variants`** and **`CN altered genes`** sections will contain information about expression levels of the mutated genes, genes located within detected structural variants (SVs) and copy-number (CN) altered regions, respectively. Genes will be ordered by increasing *variants* `TIER`, *SV* `score` and `CN` *value*, resepctively, and then by decreasing absolute values in the `Patient` vs selected `dataset` column. Moreover, gene fusions detected in WTS data and reported in **`Fusion genes`** section will be first ordered based on the evidence from genome-based data (`DNA support (gene A/B)` columns).

The *[umccrise](https://github.com/umccr/umccrise)* files are expected to be organised following the folder structure below

```
|
|____umccrised (User-defined)
  |____[SampleName]
    |____pcgr
    | |____[SampleName]-somatic.pcgr.snvs_indels.tiers.tsv
    |____purple
    | |____[SampleName].purple.gene.cnv
    |____structural
      |____[SampleName]-manta.tsv
```

```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset pdac  --count_file $(pwd)/../data/test_data/final/test_sample_WTS/test_sample_WTS-ready.counts  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report  --umccrise $(pwd)/../data/test_data/umccrised/test_sample_WGS
```

>The interactive HTML report named `test_sample_WTS.pdac.RNAseq_report.html` will be created in `data/test_data/final/test_sample_WTS/RNAseq_report` folder.

#### 3. Add clinical information

For samples derived from subjects, for which clinical information is available, a treatment regimen timeline can be added to the *Transcriptome Patient Summary* report. This can be added by specifying location of an excel spreadsheet, containing clinical information, in the `--clinical_info` argument. In this spreadsheet, at least one of the following columns is expected: `NEOADJUVANT REGIMEN`, `ADJUVANT REGIMEN`, `FIRST LINE REGIMEN`, `SECOND LINE REGIMEN` or `THIRD LINE REGIMEN`.


```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset pdac  --count_file $(pwd)/../data/test_data/final/test_sample_WTS/test_sample_WTS-ready.counts  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report  --umccrise $(pwd)/../data/test_data/umccrised/test_sample_WGS  --clinical_info $(pwd)/../data/test_data/test_clinical_data.xlsx  --subject_id test.subject
```

>The interactive HTML report named `test_sample_WTS.pdac.RNAseq_report.html` will be created in `data/test_data/final/test_sample_WTS/RNAseq_report` folder.


### Output

Will be described soon...