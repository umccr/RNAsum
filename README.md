# RNA-seq report

This reporting tool is designed to post-process, summarise and visualise an output from *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)*. Its main application is to complement genome-based findings and to provide additional evidence for detected alterations.


## Table of contents

<!-- vim-markdown-toc GFM -->
* [Installation](#installation)
* [Workflow](#workflow)
* [Reference data](#reference-data)
    * [External reference cohorts](#external-reference-cohorts)
    * [Internal reference cohort](#internal-reference-cohort)
* [Input data](#input-data)
  * [WTS](#wts)
  * [WGS](#wgs)
* [Usage](#usage)
  * [Arguments](#arguments)
  * [Examples](#examples)
  	 * [Required arguments only](#1-required-arguments-only)
  	 * [Add genome-based results](#2-add-genome-based-results)
  	 * [Add clinical information](#3-add-clinical-information)
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

1. Process the RNA-seq data from *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* including per-gene **[read counts](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)** and  **[gene fusions](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv)**. The [gene fusion](./fusions) candidates are re-quantified and the read counts are [normalised, transformed](img/counts_post-processing_scheme.png) and [converted](img/Z-score_transformation_gene_wise.png) into a scale that allows to present the sample's expression measurements in the context of the [reference cohorts](#reference-data).

2. Feed in **genome-based findings** from whole-genome sequencing (WGS) data to focus on genes of interest and provide additional evidence for dysregulation of mutated genes, or genes located within detected structural variants (SVs) or copy-number (CN) altered regions. The RNA-seq report pipeline is designed to be compatible with WGS patient report based on [umccrise](https://github.com/umccr/umccrise) pipeline output.


3. Collate results with knowledge derived from public **databases** and **reference cohorts** to provide additional source of evidence, e.g. about variants clinical significance or potential druggable targets, for the genes of interest and to get an idea about their expression levels in other cancer [patient cohorts](#reference-data).

4. The final product is a html-based **interactive report** with searchable tables and plots presenting expression levels of the genes of interest genes. The report consist of several sections described in [report_structure.md](report_structure.md).


## Reference data

The reference expression data is availbale for **33 cancer types** and were derived from [external](#external-reference-cohorts) ([TCGA](https://tcga-data.nci.nih.gov/)) and [internal](#internal-reference-cohort) (UMCCR) resources.


### External reference cohorts

In order to explore expression changes in queried sample we have built a high-quality pancreatic cancer reference cohort. 

Depending on the tissue from which the patient's sample was taken, one of **33 cancer datasets** from [TCGA](https://tcga-data.nci.nih.gov/) can be used as a reference cohort for comparing expression changes in genes of interest in investigated sample. The available cancer types are listed in [TCGA projects summary table](./TCGA_projects_summary.md). These datasets have been processed using methods described in [TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/blob/master/expression/README.md#gdc-counts-data) repository. The dataset of interest can be specified by using one of the [TCGA](https://portal.gdc.cancer.gov/) project IDs (`Project` column) for the `--dataset` argument in *[RNAseq_report.R](./rmd_files/RNAseq_report.R)* script (see [Arguments](./README.md#arguments) section). 

Of note, each dataset was **cleaned** based on the quality metrics provided in the *Merged Sample Quality Annotations* file **[merged_sample_quality_annotations.tsv](http://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf)** from [TCGA Pan-Cancer Clinical Data Resource](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018) (see [TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/tree/master/expression/README.md#data-clean-up) repository for more details, including sample inclusion criteria).


### Internal reference cohort

The [TCGA datasets](./TCGA_projects_summary.md) are expected to demonstrate prominent batch effects when compared to the input WTS data due to differences in applied experimental procedures and analytical pipelines. Moreover, TCGA data may include samples from tissue material of lower quality and cellularity compared to samples processed using local protocols. To address these issues, we have built a high-quality internal reference cohort processed using the same pipelines as the investigated sample. 

This set of **40 pancreatic cancer samples** is based on WTS data generated at **UMCCR** and processed with **[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)** *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* to minimise potential batch effects between queried sample and the reference cohorts and to make sure the data are comparable. The internal reference cohort assembly is summarised in [Pancreatic-data-harmonization](https://github.com/umccr/Pancreatic-data-harmonization/tree/master/expression/in-house) repository.

- **USED FOR BATCH-EFFECTS CORRECTION** (regardless of the input sample tissue origin) to minimise technical-related variation in the data (present due to differences in protocoles used for data generation and processing)

- **ALSO** use as a cancer group to be used to compare per-gene expression values and report in the summary tables for PDACs





## Input data

The pipeline accepts [WTS](#wts) data processed by *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)*. Additionally, the WTS data can be integrated with [WGS](#wgs)-based data processed using *[umccrise](https://github.com/umccr/umccrise)* pipeline. In the latter case, the genome-based findings from corresponding sample are incorporated into the report and are used as a primary source for expression profiles prioritisation.


### WTS 

The following output files from *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* are use in the pipeline:

Tool | Output file | Example | Required
------------ | ------------ | ------------ | ------------
[kallisto](https://pachterlab.github.io/kallisto/about) | Quantified abundances of transcripts | [abundance.tsv](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv) | **Yes**
[kallisto](https://pachterlab.github.io/kallisto/about) | Re-quantified abundances of transcripts, including fusion transcripts identified by [pizzly](https://github.com/pmelsted/pizzly) | [abundance.tsv](./data/test_data/final/test_sample_WTS/kallisto/quant_pizzly_post/abundance.tsv) | No
[pizzly](https://github.com/pmelsted/pizzly) | List of detected fusion genes | [test_sample_WTS-flat.tsv](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv) | No
[clinker](https://github.com/Oshlack/Clinker) | Plots of detected fusion genes using [pizzly](https://github.com/pmelsted/pizzly) | [EIF4A2_PTMA.pdf](./data/test_data/final/test_sample_WTS/clinker/EIF4A2_PTMA.pdf) | No

<br />

These files are expected to be organised following the folder structure below

```
|
|____final
  |____<SampleName>
    |____kallisto
    | |____abundance.tsv
    | |____quant_pizzly_post
    |   |____abundance.tsv
    |____pizzly
    | |____<SampleName>-flat.tsv
    |____clinker
      |____<GeneA_GeneB>.pdf
      |____...
      |____<GeneY_GeneZ>.pdf
```

### WGS

The following output files from *[umccrise](https://github.com/umccr/umccrise)* are accepted in the pipeline:

Tool | Output file | Example | Required
------------ | ------------ | ------------ | ------------
[PCGR](https://github.com/sigven/pcgr) | List of detected and annotated single-nucleotide variants (SNVs) and indels | [test_sample_WGS-somatic.pcgr.snvs_indels.tiers.tsv](./data/test_data/umccrised/test_sample_WGS/pcgr/test_sample_WGS-somatic.pcgr.snvs_indels.tiers.tsv) | No
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) | List of genes involved in CN altered regions | [test_sample_WGS.purple.gene.cnv](./data/test_data/umccrised/test_sample_WGS/purple/test_sample_WGS.purple.gene.cnv) | No
[Manta](https://github.com/Illumina/manta) | List of genes involved in SV regions | [test_sample_WGS-sv-prioritize-manta-pass.tsv](./data/test_data/umccrised/test_sample_WGS/structural/test_sample_WGS-sv-prioritize-manta-pass.tsv) | No

<br />


These files are expected to be organised following the folder structure below

```
|
|____umccrised
  |____<SampleName>
    |____pcgr
    | |____<SampleName>-somatic.pcgr.snvs_indels.tiers.tsv
    |____purple
    | |____<SampleName>.purple.gene.cnv
    |____structural
      |____<SampleName>-manta.tsv
```


## Usage

To run the pipeline execure the *[RNAseq_report.R](./rmd_files/RNAseq_report.R)* script. This script catches the arguments from the command line and passes them to the *[RNAseq_report.Rmd](./rmd_files/RNAseq_report.Rmd)* script to produce the interactive HTML report.

### Arguments

Argument | Description | Required
------------ | ------------ | ------------
--sample_name | Desired sample name to be presented in the report | **Yes**
--dataset | Tissue from which the sample is derived. Available options are [TCGA](https://tcga-data.nci.nih.gov/) project IDs listed in [TCGA projects summary table](./TCGA_projects_summary.md) `Project` column | **Yes**
--count_file | Location and name of the *[kallisto](https://pachterlab.github.io/kallisto/about)* read count file from *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* (see [example](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)) | **Yes**
--report_dir | Desired location for the report | **Yes**
--transform | Transformation method for converting read counts. Available options are: `CPM` (defualt) and `TPM` | No
--norm | Normalisation method. Available options are: `TMM` (for *CPM-transformed* data, defualt) and `quantile` (for *TPM-transformed* data) | No
--filter | Filtering out low expressed genes. Available options are: `TRUE` (defualt) and `FALSE` | No
--log | Log (base 2) transform data before normalisation. Available options are: `TRUE` (defualt) and `FALSE` | No
--scaling | Apply [`gene-wise`](img/Z-score_transformation_gene_wise.png) (default) or [`group-wise`](img/Z-score_transformation_group_wise.png) data scaling | No
--umccrise | Location of the corresponding *[umccrise](https://github.com/umccr/umccrise)* output (including [PCGR](https://github.com/sigven/pcgr) (see [example](./data/test_data/umccrised/test_sample_WGS/pcgr/test_sample_WGS-somatic.pcgr.snvs_indels.tiers.tsv)), [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) (see [example](./data/test_data/umccrised/test_sample_WGS/purple/test_sample_WGS.purple.gene.cnv)) and [Manta](https://github.com/Illumina/manta) (see [example](./data/test_data/umccrised/test_sample_WGS/structural/test_sample_WGS-sv-prioritize-manta-pass.tsv)) output files) from genome-based data | No
--pcgr_tier | Tier threshold for reporting variants reported in [PCGR](https://github.com/sigven/pcgr) (if [PCGR](https://github.com/sigven/pcgr) results are available, default is `3`) | No
--cn_loss | CN threshold value to classify genes within lost regions (if CN results from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) are available, default is `1.5`) | No
--cn_gain | CN threshold value to classify genes within gained regions (if CN results from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) are available, default is `3`) | No
--clinical_info | Location of *xslx* file with clinical information (see [example](./data/test_data/test_clinical_data.xlsx) ) | No
--subject_id | Subject ID required to match sample with clinical information (if available) | No
--plots_mode | Plotting mode. Available options: `interactive` (all possible plots will be interactive; default) and `semi-interactive` (only plots in `Input data`, `CN altered genes`, `Immune markers` and `HRD genes` sections will be interactive) or `static` (all plots will be static) | No
--hide_code_btn | Hide the *Code* button allowing to show/hide code chunks in the final HTML report. Available options are: `TRUE` (defualt) and `FALSE` | No
--grch_version | Human reference genome version used for genes annotation. Available options: `37` (default) and `38` | No

<br />

**Packages**: *[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)*, *[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)*, *[rapportools](https://cran.r-project.org/web/packages/rapportools/index.html)*, 
*[optparse](https://cran.r-project.org/web/packages/optparse/index.html)*, *[openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html)*, *[readr](https://cran.r-project.org/web/packages/readr/index.html)*, *[tidyverse](https://www.tidyverse.org/)*, *[dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)*, *[tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)*, *[rlang](https://cran.r-project.org/web/packages/rlang/index.html)*, *[DT](https://cran.r-project.org/web/packages/DT/index.html)*, *[kableExtra](https://cran.r-project.org/web/packages/kableExtra/index.html)*, *[matrixStats](https://cran.rstudio.com/web/packages/matrixStats/index.html)*, *[tibble](https://cran.r-project.org/web/packages/tibble/index.html)*, *[knitr](https://cran.r-project.org/web/packages/knitr/index.html)*, *[plotly](https://plot.ly/r/)*, *[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)*, *[ggforce](https://cran.r-project.org/web/packages/ggforce/index.html)*, *[pdftools](https://cran.r-project.org/web/packages/pdftools/index.html)*, *[png](https://cran.r-project.org/web/packages/png/index.html)*, *[lares](https://www.rdocumentation.org/packages/lares/versions/4.4)*, *[htmltools](https://cran.r-project.org/web/packages/htmltools/index.html)*, *[htmlwidgets](https://cran.r-project.org/web/packages/htmlwidgets/index.html)*, *[devtools](https://cran.r-project.org/web/packages/devtools/index.html)*, *[tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)*, *[EnsDb.Hsapiens.v75](http://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v75.html)* (*[EnsDb.Hsapiens.v86](http://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v86.html)*)\*, *[BSgenome.Hsapiens.UCSC.hg19](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html)* (*[BSgenome.Hsapiens.UCSC.hg38](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html)*)\*

\*  Human reference genome ***[GRCh37](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)*** (*Ensembl* based annotation version ***75***) is used for genes annotation as default. Alternatively, human reference genome [GRCh38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39) (*Ensembl* based annotation version *86*) is used when argument `grch_version` is set to `38`.


### Examples 

Below are command line use examples for generating *Transcriptome Patient Summary* report using:

1. [required arguments only](#1-required-arguments-only)
2. **[genome-based results](#2-add-genome-based-results)**
3. [clinical information](#3-add-clinical-information)

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

In this scenario, only expression levels of key **[`Cancer genes`](https://github.com/umccr/umccrise/blob/master/workflow.md#key-cancer-genes)**, **`Fusion genes`**, **`Immune markers`** and homologous recombination deficiency genes (**`HRD genes`**) will be reported. No results will be provided in ~~`Mutated genes`~~, ~~`Structural variants`~~ and ~~`CN altered genes`~~ sections. Moreover, gene fusions reported in `Fusion genes` section will not contain inforamation about evidence from genome-based data. The [TCGA](https://tcga-data.nci.nih.gov/) bladder urothelial carcinoma dataset is used as reference cohort (`--dataset BLCA`).

The input files are expected to be organised following the folder structure described in [Input data:WTS](#wts) section.

```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset BLCA  --count_file $(pwd)/../data/test_data/final/test_sample_WTS/kallisto/abundance.tsv  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report
```

>The interactive HTML report named `test_sample_WTS.blca.RNAseq_report.html` will be created in `data/test_data/final/test_sample_WTS/RNAseq_report` folder.


#### 2. Add genome-based results

This is the **preferred scenario for using** ***Transcriptome Patient Summary***, in which the genome-based findings will be used as a primary source for expression profiles prioritisation. These can be incorporated into the report by specifying location of the corresponding ***[umccrise](https://github.com/umccr/umccrise)*** output files (including results from [PCGR](https://github.com/sigven/pcgr), [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) and [Manta](https://github.com/Illumina/manta)) using `--umccrise` argument. The **`Mutated genes`**, **`Structural variants`** and **`CN altered genes`** sections will contain information about expression levels of the mutated genes, genes located within detected structural variants (SVs) and copy-number (CN) altered regions, respectively. The results in the **`Fusion genes`** section will be ordered based on the evidence from genome-based data. The [TCGA](https://tcga-data.nci.nih.gov/) cervical squamous cell carcinoma dataset is used as reference cohort (`--dataset CESC `).

The *[umccrise](https://github.com/umccr/umccrise)* output files are expected to be organised following the folder structure described in [Input data:WGS](#wgs) section.


```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset CESC  --count_file $(pwd)/../data/test_data/final/test_sample_WTS/kallisto/abundance.tsv  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report  --umccrise $(pwd)/../data/test_data/umccrised/test_sample_WGS
```

>The interactive HTML report named `test_sample_WTS.cesc.RNAseq_report.html` will be created in `data/test_data/final/test_sample_WTS/RNAseq_report` folder.

#### 3. Add clinical information

For samples derived from subjects, for which clinical information is available, a treatment regimen timeline can be added to the *Transcriptome Patient Summary* report. This can be added by specifying location of a relevant excel spreadsheet (see example [test_clinical_data.xlsx](./data/test_data/test_clinical_data.xlsx)) using the `--clinical_info` argument. In this spreadsheet, at least one of the following columns is expected: `NEOADJUVANT REGIMEN`, `ADJUVANT REGIMEN`, `FIRST LINE REGIMEN`, `SECOND LINE REGIMEN` or `THIRD LINE REGIMEN`, along with `START` and `STOP` dates of corresponding treatments. The [TCGA](https://tcga-data.nci.nih.gov/) pancreatic adenocarcinoma dataset is used as reference cohort (`--dataset PAAD `).


```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset PAAD  --count_file $(pwd)/../data/test_data/final/test_sample_WTS/kallisto/abundance.tsv  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAseq_report  --umccrise $(pwd)/../data/test_data/umccrised/test_sample_WGS  --clinical_info $(pwd)/../data/test_data/test_clinical_data.xlsx  --subject_id test.subject
```

>The interactive HTML report named `test_sample_WTS.paad.RNAseq_report.html` will be created in `data/test_data/final/test_sample_WTS/RNAseq_report` folder.


### Output

The generated html-based ***Transcriptome Patient Summary*** **report** includes searchable tables and interactive plots presenting expression levels of altered genes, as well as links to public resources describing the genes of interest. The report consist of several sections described more in detail in [report_structure.md](report_structure.md):

* [Input data](report_structure.md#input-data)
* [Clinical information\*](report_structure.md#clinical-information)
* [Mutated genes\**](report_structure.md#mutated-genes)
* [Cancer genes](report_structure.md#cancer-genes)
* [Fusion genes](report_structure.md#fusion-genes)
* [Structural variants\**](report_structure.md#structural-variants)
* [CN altered genes\**](report_structure.md#cn-altered-genes)
* [Immune markers](report_structure.md#immune-markers)
* [HRD genes](report_structure.md#hrd-genes)
* [Drug matching](report_structure.md#drug-matching)
* [Addendum](report_structure.md#addendum)

\* if clinical information is available; see `--clinical_info` [argument](#arguments) <br />
\** if genome-based results are available; see `--umccrise ` [argument](#arguments)
