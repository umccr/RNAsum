# RNAsum

RNA-seq reporting workflow designed to post-process, summarise and visualise an output from *[bcbio-nextgen RNA-seq](https://bcbio-nextgen.readthedocs.io/en/latest/contents/bulk_rnaseq.html)* or *[Dragen RNA](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html)* pipelines. Its main application is to complement genome-based findings from [umccrise](https://github.com/umccr/umccrise) pipeline and to provide additional evidence for detected alterations.


## Table of contents

<!-- vim-markdown-toc GFM -->
* [Installation](#installation)
* [Workflow](#workflow)
* [Reference data](#reference-data)
    * [External reference cohorts](#external-reference-cohorts)
    * [Internal reference cohort](#internal-reference-cohort)
* [Input data](#input-data)
  * [WTS](#wts)
     * [bcbio-nextgen](#bcbio-nextgen)
     * [Dragen RNA](#dragen-rna)
  * [WGS](#wgs)
* [Usage](#usage)
  * [Arguments](#arguments)
  * [Examples](#examples)
  	 * [WTS data only](#1-wts-data-only)
  	 * [WTS and WGS data](#2-wts-and-wgs-data)
  	 * [WTS WGS and clinical data](#3-wts-wgs-and-clinical-data)
  * [Output](#output)
* [Docker](#docker)

<!-- vim-markdown-toc -->


## Installation

Run the following to create a directory "rnasum" and install into it

```
mkdir rnasum
cd rnasum
source <(curl -s https://raw.githubusercontent.com/umccr/RNAseq-Analysis-Report/master/install.sh)
```

It will generate `load_rnasum.sh` script that can be sourced to load the `rnasum` environment:

```
source load_rnasum.sh
```


## Workflow

The pipeline consist of five main components illustrated and briefly described below. See the [workflow.md](workflow.md) for the complete description of the **[data processing workflow](workflow.md)**.

<img src="img/RNAsum_workflow.png" width="100%"> 

<br/>

1. Collect patient sample WTS data from *[bcbio-nextgen RNA-seq](https://bcbio-nextgen.readthedocs.io/en/latest/contents/bulk_rnaseq.html)* or *[Dragen RNA](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html)* pipeline including per-gene **[read counts](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)** and **[gene fusions](./data/test_data/final/test_sample_WTS/arriba/fusions.tsv)**.

2. Add expression data from [reference cohorts](#reference-data) to get an idea about expression levels of genes of interest in other cancer [patient cohorts](#reference-data). The read counts are [normalised, transformed](img/counts_post-processing_scheme.png) and [converted](img/Z-score_transformation_gene_wise.png) into a scale that allows to present the sample's expression measurements in the context of the [reference cohorts](#reference-data).

3. Feed in **genome-based findings** from whole-genome sequencing (WGS) data to focus on genes of interest and provide additional evidence for dysregulation of mutated genes, or genes located within detected structural variants (SVs) or copy-number (CN) altered regions. The ***RNAsum*** pipeline is designed to be compatible with WGS patient report based on ***[umccrise](https://github.com/umccr/umccrise)*** pipeline output.


4. Collate results with knowledge derived from in-house resources and public **databases** to provide additional source of evidence for clinical significance of altered genes, e.g. to flag variants with clinical significance or potential druggable targets.

5. The final product is a html-based **interactive report** with searchable tables and plots presenting expression levels of the genes of interest genes. The report consist of several sections described in [report_structure.md](report_structure.md).


## Reference data

The reference expression data is availbale for **33 cancer types** and were derived from [external](#external-reference-cohorts) ([TCGA](https://tcga-data.nci.nih.gov/)) and [internal](#internal-reference-cohort) ([UMCCR](https://research.unimelb.edu.au/centre-for-cancer-research/our-research/precision-oncology-research-group)) resources.


### External reference cohorts

In order to explore expression changes in queried sample we have built a high-quality pancreatic cancer reference cohort. 

Depending on the tissue from which the patient's sample was taken, one of **33 cancer datasets** from [TCGA](https://tcga-data.nci.nih.gov/) can be used as a reference cohort for comparing expression changes in genes of interest in investigated sample. Additionally, 10 samples from each of the 33 datasets were combined to create **[Pan-Cancer dataset](./TCGA_projects_summary.md#pan-cancer-dataset)**, and for some cohorts **[extended sets](./TCGA_projects_summary.md#extended-datasets)** are also available. All available datasets are listed in **[TCGA projects summary table](./TCGA_projects_summary.md)**. These datasets have been processed using methods described in [TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/blob/master/expression/README.md#gdc-counts-data) repository. The dataset of interest can be specified by using one of the [TCGA](https://portal.gdc.cancer.gov/) project IDs (`Project` column) for the `--dataset` argument in *[RNAseq_report.R](./rmd_files/RNAseq_report.R)* script (see [Arguments](./README.md#arguments) section). 

###### Note

Each dataset was **cleaned** based on the quality metrics provided in the *Merged Sample Quality Annotations* file **[merged_sample_quality_annotations.tsv](http://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf)** from [TCGA PanCanAtlas initiative webpage](https://gdc.cancer.gov/about-data/publications/pancanatlas) (see [TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/tree/master/expression/README.md#data-clean-up) repository for more details, including sample inclusion criteria).


### Internal reference cohort

The publically available TCGA datasets are expected to demonstrate prominent [batch effects](https://www.ncbi.nlm.nih.gov/pubmed/20838408) when compared to the in-house WTS data due to differences in applied experimental procedures and analytical pipelines. Moreover, TCGA data may include samples from tissue material of lower quality and cellularity compared to samples processed using local protocols. To address these issues, we have built a high-quality internal reference cohort processed using the same pipelines as input data (see [Data pre-processing](./workflow.md#data-pre-processing) section on the [workflow](./workflow.md) page). 

This internal reference set of **40 pancreatic cancer samples** is based on WTS data generated at **[UMCCR](https://research.unimelb.edu.au/centre-for-cancer-research/our-research/precision-oncology-research-group)** and processed with **[bcbio-nextgen RNA-seq](https://bcbio-nextgen.readthedocs.io/en/latest/contents/bulk_rnaseq.html)** pipeline to minimise potential batch effects between investigated samples and the reference cohort and to make sure the data are comparable. The internal reference cohort assembly is summarised in [Pancreatic-data-harmonization](https://github.com/umccr/Pancreatic-data-harmonization/tree/master/expression/in-house) repository.

###### Note

The are two rationales for using the internal reference cohort:

1. In case of **pancreatic cancer samples** this cohort is used (I) in ***batch effects correction***, as well as (II) as a reference point for ***comparing per-gene expression levels*** observed in the investigated single-subject data and data from other pancreatic cancer patients.

2. In case of samples from **any cancer type** the data from the internal reference cohort is used in ***batch effects correction*** procedure performed to minimise technical-related variation in the data.

## Input data

The pipeline accepts [WTS](#wts) data processed by *[bcbio-nextgen RNA-seq](https://bcbio-nextgen.readthedocs.io/en/latest/contents/bulk_rnaseq.html)* or *[Dragen RNA](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html)*  pipeline. Additionally, the WTS data can be integrated with [WGS](#wgs)-based data processed using *[umccrise](https://github.com/umccr/umccrise)* pipeline. In the latter case, the genome-based findings from corresponding sample are incorporated into the report and are used as a primary source for expression profiles prioritisation.


### WTS 

The only required WTS input data are **read counts** provided in quantification file from either *[bcbio-nextgen RNA-seq](https://bcbio-nextgen.readthedocs.io/en/latest/contents/bulk_rnaseq.html)* or *[Dragen RNA](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html)* pipeline.

#### bcbio-nextgen

The **read counts** are provided within quantification file from [kallisto](https://pachterlab.github.io/kallisto/about) (see example *[abundance.tsv](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)* file and its [description](https://pachterlab.github.io/kallisto/starting#results)). The per-transcript abundances are reported in *estimated counts* (`est_counts`) and in *[Transcripts Per Million](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)* (`tpm`), which are then converted to per-gene estimates. Additionally, a list of **[fusion genes](./fusions)** detected by [arriba](https://arriba.readthedocs.io/en/latest/) and [pizzly](https://github.com/pmelsted/pizzly) can be provided (see example *[fusions.tsv](./data/test_data/final/test_sample_WTS/arriba/fusions.tsv)* and *[test_sample_WTS-flat.tsv](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv)*). 

Table below lists all input data accepted in the pipeline:

Input file | Tool | Example | Required
------------ | ------------ | ------------ | ------------
Quantified **abundances** of transcripts | [kallisto](https://pachterlab.github.io/kallisto/about) | [abundance.tsv](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv) | **Yes**
List of detected **fusion genes** | [arriba](https://arriba.readthedocs.io/en/latest/) </br> [pizzly](https://github.com/pmelsted/pizzly) | [fusions.tsv](./data/test_data/final/test_sample_WTS/arriba/fusions.tsv) </br> [test_sample_WTS-flat.tsv](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv) | No
Plots of detected **fusion genes** using [arriba](https://arriba.readthedocs.io/en/latest/) | [arriba](https://arriba.readthedocs.io/en/latest/) | [fusions.pdf](./data/test_data/final/test_sample_WTS/arriba/fusions.pdf) | No

<br />

These files are expected to be organised following the folder structure below

```
|
|____<SampleName>
  |____kallisto
  | |____abundance.tsv
  |____pizzly
  | |____<SampleName>-flat.tsv
  |____arriba
    |____fusions.pdf
    |____fusions.tsv
```

###### Note

[Fusion genes](./fusions) detected by [pizzly](https://github.com/pmelsted/pizzly) are expected to be listed in the [flat table](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv). By default two output tables are provided: (1) *\<sample_name\>-flat.tsv* listing all gene fusion candidates and (2) *\<sample_name\>-flat-filtered.tsv* listing only gene fusions remaining after filtering step. However, this workflow makes use of gene fusions listed in the **unfiltered** [pizzly](https://github.com/pmelsted/pizzly) output file (see example [test_sample_WTS-flat.tsv](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv)) since it was noted that some genuine fusions (based on WGS data and curation efforts) are excluded in the filtered [pizzly](https://github.com/pmelsted/pizzly) output file.


#### Dragen RNA

The **read counts** are provided within quantification file from [salmon](https://salmon.readthedocs.io/en/latest/salmon.html) (see example *[TEST.quant.sf](./data/test_data/stratus/test_sample_WTS/TEST.quant.sf)* file and its [description](https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats)). The per-transcript abundances are reported in *estimated counts* (`NumReads `) and in *[Transcripts Per Million](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)* (`TPM `), which are then converted to per-gene estimates. Additionally, a list of **[fusion genes](./fusions)** can be provided (see example *[TEST.fusion_candidates.final](./data/test_data/stratus/test_sample_WTS/TEST.fusion_candidates.final)*). 

Table below lists all input data accepted in the pipeline:

Input file | Tool | Example | Required
------------ | ------------ | ------------ | ------------
Quantified **abundances** of transcripts | [salmon](https://salmon.readthedocs.io/en/latest/salmon.html) | [TEST.quant.sf](./data/test_data/stratus/test_sample_WTS/TEST.quant.sf) | **Yes**
List of detected **fusion genes** | [Dragen RNA](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html) | [TEST.fusion_candidates.final](./data/test_data/stratus/test_sample_WTS/TEST.fusion_candidates.final) | No

<br />

These files are expected to be organised following the folder structure below

```
|
|____<SampleName>
  |____<SampleName>quant.sf
  |____<SampleName>.fusion_candidates.final
```


<br />

### WGS

The following *[umccrise](https://github.com/umccr/umccrise)* output files are accepted as input data in the pipeline:

Input file | Tool | Example | Required
------------ | ------------ | ------------ | ------------
List of detected and annotated single-nucleotide variants (**SNVs**) and indels | [PCGR](https://github.com/sigven/pcgr) | [test_subject__test_sample_WGS-somatic.pcgr.snvs_indels.tiers.tsv](./data/test_data/umccrised/test_subject__test_sample_WGS/pcgr/test_subject__test_sample_WGS-somatic.pcgr.snvs_indels.tiers.tsv) | No
List of genes involved in **CN** altered regions | [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) | [test_subject__test_sample_WGS.purple.gene.cnv](./data/test_data/umccrised/test_subject__test_sample_WGS/purple/test_subject__test_sample_WGS.purple.gene.cnv) | No
List of genes involved in **SV** regions | [Manta](https://github.com/Illumina/manta) | [test_subject__test_sample_WGS-sv-prioritize-manta-pass.tsv](./data/test_data/umccrised/test_subject__test_sample_WGS/structural/test_subject__test_sample_WGS-sv-prioritize-manta-pass.tsv) | No

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
--sample_name | The name of the sample to be analysed and reported | **Yes**
--bcbio_rnaseq | Location of the results folder from *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/bulk_rnaseq.html)* | **Yes***
--dragen_rnaseq | Location of the results folder from *[Dragen RNA pipeline](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html)* | **Yes***
--report_dir | Desired location for the report | **Yes**
--dataset | Dataset to be used as external reference cohort. Available options are [TCGA](https://tcga-data.nci.nih.gov/) project IDs listed in [TCGA projects summary table](./TCGA_projects_summary.md) `Project` column (default is `PANCAN`) | No
--transform | Transformation method for converting read counts. Available options are: `CPM` (default) and `TPM` | No
--norm | Normalisation method. Available options are: `TMM` (default), `TMMwzp`, `RLE`, `upperquartile` or `none` for *CPM-transformed* data, and `quantile` (default) or `none` for *TPM-transformed* data | No
--batch_rm | Remove batch-associated effects between datasets. Available options are: `TRUE` (default) and `FALSE`  | No
--filter | Filtering out low expressed genes. Available options are: `TRUE` (default) and `FALSE` | No
--log | Log (base 2) transform data before normalisation. Available options are: `TRUE` (default) and `FALSE` | No
--scaling | Apply [`gene-wise`](img/Z-score_transformation_gene_wise.png) (default) or [`group-wise`](img/Z-score_transformation_group_wise.png) data scaling. Available options are: `TRUE` (default) and `FALSE` | No
--drugs | Include drug matching section in the report. Available options are: `TRUE` and `FALSE` (default) | No
--immunogram | Include immunogram in the report. Available options are: `TRUE` and `FALSE` (default) | No
--umccrise | Location of the corresponding *[umccrise](https://github.com/umccr/umccrise)* output (including [PCGR](https://github.com/sigven/pcgr) (see [example](./data/test_data/umccrised/test_subject__test_sample_WGS/pcgr/test_subject__test_sample_WGS-somatic.pcgr.snvs_indels.tiers.tsv)), [Manta](https://github.com/Illumina/manta) (see [example](./data/test_data/umccrised/test_subject__test_sample_WGS/structural/test_subject__test_sample_WGS-sv-prioritize-manta-pass.tsv)) and [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) (see [example](./data/test_data/umccrised/test_subject__test_sample_WGS/purple/test_subject__test_sample_WGS.purple.gene.cnv)) output files) from genome-based data | No
--pcgr_tier | [Tier](https://pcgr.readthedocs.io/en/latest/tier_systems.html#tier-model-2-pcgr-acmg) threshold for reporting variants reported in [PCGR](https://github.com/sigven/pcgr) (if [PCGR](https://github.com/sigven/pcgr) results are available, default is `4`) | No
--pcgr_splice_vars | Include non-coding `splice_region_variant`s reported in [PCGR](https://github.com/sigven/pcgr) (if [PCGR](https://github.com/sigven/pcgr) results are available). Available options are: `TRUE` (default) and `FALSE` | No
--cn_loss | CN threshold value to classify genes within lost regions (if CN results from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) are available, default is `5th percentile` of all CN values) | No
--cn_gain | CN threshold value to classify genes within gained regions (if CN results from [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) are available, default is `95th percentile` of all CN values) | No
--clinical_info | Location of *xslx* file with clinical information (see [example](./data/test_data/test_clinical_data.xlsx) ) | No
--clinical_id | ID required to match sample with the subject clinical information (if available) | No
--subject_id | Subject ID. Note, if `umccrise` is specified (flag `--umccrise`) then subject ID is extracted from `umccrise` output files and used to overwrite this argument | No
--sample_source | Source of investigated sample (e.g. fresh frozen tissue, organoid; for annotation purposes only) | No
--sample_name_mysql | Desired sample name for MySQL insert command. By default value in `--sample_name` is used | No
--project | Project name (for annotation purposes only) | No
--top_genes | The number of top ranked genes to be presented (default is `5`) | No
--dataset_name_incl | Include dataset in the report name. Available options are: `TRUE` and `FALSE` (default) | No
--save_tables | Save interactive summary tables as HTML files. Available options are: `TRUE` (default) and `FALSE` | No
--hide_code_btn | Hide the *Code* button allowing to show/hide code chunks in the final HTML report. Available options are: `TRUE` (default) and `FALSE` | No
--grch_version | Human reference genome version used for genes annotation. Available options: `37` and `38` (default) | No

<br />

\* Location of the results folder from either *[bcbio-nextgen RNA-seq](https://bcbio-nextgen.readthedocs.io/en/latest/contents/bulk_rnaseq.html)* or *[Dragen RNA](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html)* pipeline is required.

**Packages**: required packages are listed in [environment.yaml](envm/environment.yaml) file.

###### Note

Human reference genome ***[GRCh38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39)*** (*Ensembl* based annotation version ***86***) is used for genes annotation as default. Alternatively, human reference genome [GRCh37](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/) (*Ensembl* based annotation version *75*) is used when argument `grch_version` is set to `37`.


### Examples

Below are command line use examples for generating *Patient Transcriptome Summary* report using:

1. [WTS data only](#1-wts-data-only)
2. **[WTS and WGS data](#2-wts-and-wgs-data)**
3. [WTS WGS and clinical data](#3-wts-wgs-and-clinical-data)

###### Note

* make sure that the created *conda* environment (see [Installation](#installation) section) is  activated

```
conda activate rnasum
```

* *[RNAseq_report.R](./rmd_files/RNAseq_report.R)* script (see the beginning of [Usage](#usage) section) should be executed from [rmd_files](./rmd_files) folder

```
cd rmd_files
```

* example data is provided in [data/test_data](./data/test_data) folder
* Usually the data processing and report generation would take less than **20 minutes** using **16GB RAM** memory and **1 CPU**


#### 1. WTS data only

In this scenario, only [WTS](#wts) data will be used and only expression levels of key **[`Cancer genes`](https://github.com/umccr/umccrise/blob/master/workflow.md#key-cancer-genes)**, **`Fusion genes`**, **`Immune markers`** and homologous recombination deficiency genes (**`HRD genes`**) will be reported. Moreover, gene fusions reported in `Fusion genes` section will not contain inforamation about evidence from genome-based data. A subset of the [TCGA](https://tcga-data.nci.nih.gov/) pancreatic adenocarcinoma dataset is used as reference cohort (`--dataset TEST `).

The input files are expected to be organised following the folder structure described in [Input data:WTS](#wts) section.


##### bcbio-nextgen

```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset TEST  --bcbio_rnaseq $(pwd)/../data/test_data/final/test_sample_WTS  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAsum  --save_tables FALSE
```

>The interactive HTML report named `test_sample_WTS.RNAsum.html` will be created in `data/test_data/final/test_sample_WTS/RNAsum` folder.

##### Dragen RNA

```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset TEST  --dragen_rnaseq $(pwd)/../data/test_data/stratus/test_sample_WTS  --report_dir $(pwd)/../data/test_data/stratus/test_sample_WTS/RNAsum  --save_tables FALSE
```

>The interactive HTML report named `test_sample_WTS.RNAsum.html` will be created in `data/test_data/stratus/test_sample_WTS/RNAsum` folder.


#### 2. WTS and WGS data

This is the **most frequent and preferred case**, in which the [WGS](#wgs)-based findings will be used as a primary source for expression profiles prioritisation. The genome-based results can be incorporated into the report by specifying location of the corresponding ***[umccrise](https://github.com/umccr/umccrise)*** output files (including results from [PCGR](https://github.com/sigven/pcgr), [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) and [Manta](https://github.com/Illumina/manta)) using `--umccrise` argument. The **`Mutated genes`**, **`Structural variants`** and **`CN altered genes`** sections will contain information about expression levels of the mutated genes, genes located within detected structural variants (SVs) and copy-number (CN) altered regions, respectively. The results in the **`Fusion genes`** section will be ordered based on the evidence from genome-based data. A subset of the [TCGA](https://tcga-data.nci.nih.gov/) pancreatic adenocarcinoma dataset is used as reference cohort (`--dataset TEST `).

The *[umccrise](https://github.com/umccr/umccrise)* output files are expected to be organised following the folder structure described in [Input data:WGS](#wgs) section.

##### bcbio-nextgen

```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset TEST  --bcbio_rnaseq $(pwd)/../data/test_data/final/test_sample_WTS  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAsum  --umccrise $(pwd)/../data/test_data/umccrised/test_subject__test_sample_WGS  --save_tables FALSE
```

>The interactive HTML report named `test_sample_WTS.RNAsum.html` will be created in `data/test_data/final/test_sample_WTS/RNAsum` folder.

##### Dragen RNA

```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset TEST  --dragen_rnaseq $(pwd)/../data/test_data/stratus/test_sample_WTS  --report_dir $(pwd)/../data/test_data/stratus/test_sample_WTS/RNAsum  --umccrise $(pwd)/../data/test_data/umccrised/test_subject__test_sample_WGS  --save_tables FALSE
```

>The interactive HTML report named `test_sample_WTS.RNAsum.html` will be created in `data/test_data/stratus/test_sample_WTS/RNAsum` folder.


#### 3. WTS WGS and clinical data

For samples derived from subjects, for which clinical information is available, a treatment regimen timeline can be added to the *Patient Transcriptome Summary* report. This can be added by specifying location of a relevant excel spreadsheet (see example [test_clinical_data.xlsx](./data/test_data/test_clinical_data.xlsx)) using the `--clinical_info` argument. In this spreadsheet, at least one of the following columns is expected: `NEOADJUVANT REGIMEN`, `ADJUVANT REGIMEN`, `FIRST LINE REGIMEN`, `SECOND LINE REGIMEN` or `THIRD LINE REGIMEN`, along with `START` and `STOP` dates of corresponding treatments. A subset of the [TCGA](https://tcga-data.nci.nih.gov/) pancreatic adenocarcinoma dataset is used as reference cohort (`--dataset TEST `).

##### bcbio-nextgen

```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset TEST  --bcbio_rnaseq $(pwd)/../data/test_data/final/test_sample_WTS  --report_dir $(pwd)/../data/test_data/final/test_sample_WTS/RNAsum  --umccrise $(pwd)/../data/test_data/umccrised/test_subject__test_sample_WGS  --clinical_info $(pwd)/../data/test_data/test_clinical_data.xlsx --save_tables FALSE
```

>The interactive HTML report named `test_sample_WTS.RNAsum.html` will be created in `data/test_data/final/test_sample_WTS/RNAsum` folder.

##### Dragen RNA

```
Rscript RNAseq_report.R  --sample_name test_sample_WTS  --dataset TEST  --dragen_rnaseq $(pwd)/../data/test_data/stratus/test_sample_WTS  --report_dir $(pwd)/../data/test_data/stratus/test_sample_WTS/RNAsum  --umccrise $(pwd)/../data/test_data/umccrised/test_subject__test_sample_WGS  --clinical_info $(pwd)/../data/test_data/test_clinical_data.xlsx  --save_tables FALSE
```

>The interactive HTML report named `test_sample_WTS.RNAsum.html` will be created in `data/test_data/stratus/test_sample_WTS/RNAsum` folder.

### Output

The pipeline generates html-based ***Patient Transcriptome Summary*** **[report](#report)** and [results](#results) folder within user-defined `output` folder:

```
|
|____<output>
  |____<SampleName>.<output>.html
  |____results
    |____exprTables
    |____glanceExprPlots
    |____...
```

#### Report

The generated html-based ***Patient Transcriptome Summary*** **report** includes searchable tables and interactive plots presenting expression levels of altered genes, as well as links to public resources describing the genes of interest. The report consist of several sections, including:

* [Input data](report_structure.md#input-data)
* [Clinical information\*](report_structure.md#clinical-information)
* [Findings summary](report_structure.md#findings-summary)
* [Mutated genes\**](report_structure.md#mutated-genes)
* [Fusion genes](report_structure.md#fusion-genes)
* [Structural variants\**](report_structure.md#structural-variants)
* [CN altered genes\**](report_structure.md#cn-altered-genes)
* [Immune markers](report_structure.md#immune-markers)
* [HRD genes](report_structure.md#hrd-genes)
* [Cancer genes](report_structure.md#cancer-genes)
* [Drug matching](report_structure.md#drug-matching)
* [Addendum](report_structure.md#addendum)

\* if clinical information is available; see `--clinical_info` [argument](#arguments) <br />
\** if genome-based results are available; see `--umccrise` [argument](#arguments)

Detailed description of the **[report structure](report_structure.md)**, including **[results prioritisation](report_structure.md)** and **[visualisation](report_structure.md)** is available in [report_structure.md](report_structure.md).
 
#### Results

The `results` folder contains intermediate files, including plots and tables that are presented in the [report](#report).

## Docker

 - Pull ready to run docker image from DockerHub

 `docker pull umccr/rnasum:0.3.2`
 
 - An example command to use this pulled docker container is:

 ```
 docker run --rm -v /path/to/RNAseq-report/RNAseq-Analysis-Report/envm/wts-report-wrapper.sh:/work/test.sh -v /path/to/RNAseq-report/RNAseq-Analysis-Report/data:/work c18db89d3093 /work/test.sh
 ```
 
 - Assumptions

 	- You are running the RNAsum container against the [RNAsum code](https://github.com/umccr/RNAsum/) and  [test/reference data](https://github.com/umccr/RNAsum/tree/master/data/)
