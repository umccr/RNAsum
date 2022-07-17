
-   <a href="#rnasum" id="toc-rnasum">RNAsum</a>
    -   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#workflow" id="toc-workflow">Workflow</a>
    -   <a href="#reference-data" id="toc-reference-data">Reference data</a>
        -   <a href="#external-reference-cohorts"
            id="toc-external-reference-cohorts">External reference cohorts</a>
        -   <a href="#internal-reference-cohort"
            id="toc-internal-reference-cohort">Internal reference cohort</a>
    -   <a href="#input-data" id="toc-input-data">Input data</a>
        -   <a href="#wts" id="toc-wts">WTS</a>
        -   <a href="#wgs" id="toc-wgs">WGS</a>
    -   <a href="#usage" id="toc-usage">Usage</a>
        -   <a href="#run-via-aws-batch" id="toc-run-via-aws-batch">Run via
            AWS-Batch</a>
        -   <a href="#examples" id="toc-examples">Examples</a>
        -   <a href="#output" id="toc-output">Output</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# RNAsum

`RNAsum` is an R package that can post-process, summarise and visualise
outputs from [DRAGEN
RNA](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html)
or [bcbio-nextgen
RNA-seq](https://bcbio-nextgen.readthedocs.io/en/latest/contents/bulk_rnaseq.html)
pipelines. Its main application is to complement genome-based findings
from the [umccrise](https://github.com/umccr/umccrise) pipeline and to
provide additional evidence for detected alterations.

## Installation

-   **R** package can be installed directly from the [GitHub
    source](https://github.com/umccr/RNAsum):

``` r
remotes::install_github("umccr/RNAsum") # latest master commit
remotes::install_github("umccr/RNAsum@v0.0.X") # version 0.0.X
remotes::install_github("umccr/RNAsum@abcde") # commit abcde
remotes::install_github("umccr/RNAsum#123") # PR 123
```

-   **Conda** package is available from the Anaconda [umccr
    channel](https://anaconda.org/umccr/r-rnasum):

``` bash
conda install r-rnasum==0.0.X -c umccr -c conda-forge -c bioconda
```

-   **Docker** image is available at the [GitHub Container
    Registy](https://github.com/umccr/RNAsum/pkgs/container/rnasum):

``` bash
docker pull ghcr.io/umccr/rnasum:latest
```

## Workflow

The pipeline consists of five main components illustrated and briefly
described below. For more details, see [workflow.md](/workflow.md).

<img src="img/RNAsum_workflow.png" width="100%">

1.  Collect patient **WTS data** from the `DRAGEN RNA` or
    `bcbio-nextgen RNA-seq` pipeline including per-gene [read
    counts](./inst/rawdata/test_data/final/test_sample_WTS/kallisto/abundance.tsv)
    and [gene
    fusions](./inst/rawdata/test_data/final/test_sample_WTS/arriba/fusions.tsv).

2.  Add expression data from **[reference cohorts](#reference-data)** to
    get an idea about the expression levels of genes of interest in
    other cancer patient cohorts. The read counts are [normalised,
    transformed](img/counts_post-processing_scheme.png) and
    [converted](img/Z-score_transformation_gene_wise.png) into a scale
    that allows to present the patient’s expression measurements in the
    context of the reference cohorts.

3.  Supply **genome-based findings** from whole-genome sequencing (WGS)
    data to focus on genes of interest and to provide additional
    evidence for dysregulation of mutated genes, or genes located within
    detected structural variants (SVs) or copy-number (CN) altered
    regions. `RNAsum` is designed to be compatible with WGS patient
    outputs generated from `umccrise`.

4.  Collate results with knowledge derived from in-house resources and
    public databases to provide additional sources of evidence for
    clinical significance of altered genes e.g. to flag variants with
    clinical significance or potential druggable targets.

5.  The final product is an interactive HTML report with searchable
    tables and plots presenting expression levels of the genes of
    interest. The report consists of several sections described
    [here](./report_structure.md).

## Reference data

The reference expression data are available for **33 cancer types** and
were derived from [external](#external-reference-cohorts)
([TCGA](https://tcga-data.nci.nih.gov/)) and
[internal](#internal-reference-cohort)
([UMCCR](https://research.unimelb.edu.au/centre-for-cancer-research/our-research/precision-oncology-research-group))
resources.

### External reference cohorts

In order to explore expression changes in the patient, we have built a
high-quality pancreatic cancer reference cohort.

Depending on the tissue from which the patient’s sample was taken, one
of **33 cancer datasets** from TCGA can be used as a reference cohort
for comparing expression changes in genes of interest of the patient.
Additionally, 10 samples from each of the 33 TCGA datasets were combined
to create the **[Pan-Cancer
dataset](./TCGA_projects_summary.md#pan-cancer-dataset)**, and for some
cohorts **[extended
sets](./TCGA_projects_summary.md#extended-datasets)** are also
available. All available datasets are listed in **[TCGA projects summary
table](./TCGA_projects_summary.md)**. These datasets have been processed
using methods described in the
[TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/blob/master/expression/README.md#gdc-counts-data)
repository. The dataset of interest can be specified by using one of the
TCGA project IDs for the `RNAsum` `--dataset` argument (see
[Arguments](./README.md#arguments)).

**Note**

Each dataset was **cleaned** based on the quality metrics provided in
the *Merged Sample Quality Annotations* file (download the TSV from
[here](http://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf)),
linked to from the [TCGA PanCanAtlas initiative
webpage](https://gdc.cancer.gov/about-data/publications/pancanatlas)
(see the
[TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/tree/master/expression/README.md#data-clean-up)
repository for more details, including sample inclusion criteria).

### Internal reference cohort

The publicly available TCGA datasets are expected to demonstrate
prominent [batch effects](https://www.ncbi.nlm.nih.gov/pubmed/20838408)
when compared to the in-house WTS data due to differences in applied
experimental procedures and analytical pipelines. Moreover, TCGA data
may include samples from tissue material of lower quality and
cellularity compared to samples processed using local protocols. To
address these issues, we have built a high-quality internal reference
cohort processed using the same pipelines as input data (see [data
pre-processing](./workflow.md#data-pre-processing)).

This internal reference set of **40 pancreatic cancer samples** is based
on WTS data generated at
**[UMCCR](https://research.unimelb.edu.au/centre-for-cancer-research/our-research/precision-oncology-research-group)**
and processed with the **bcbio-nextgen RNA-seq** pipeline to minimise
potential batch effects between investigated samples and the reference
cohort and to make sure the data are comparable. The internal reference
cohort assembly is summarised in the
[Pancreatic-data-harmonization](https://github.com/umccr/Pancreatic-data-harmonization/tree/master/expression/in-house)
repository.

**Note**

The are two rationales for using the internal reference cohort:

1.  In case of **pancreatic cancer samples** this cohort is used:
    -   in ***batch effects correction***
    -   as a reference point for ***comparing per-gene expression
        levels*** observed in the data of the patient of interest and
        data from other pancreatic cancer patients.
2.  In case of samples from **any cancer type** the data from the
    internal reference cohort is used in the ***batch effects
    correction*** procedure performed to minimise technical-related
    variation in the data.

## Input data

`RNAsum` accepts [WTS](#wts) data processed by the `DRAGEN RNA` or
`bcbio-nextgen RNA-seq` pipeline. Additionally, the WTS data can be
integrated with [WGS](#wgs)-based data processed using the `umccrise`
pipeline. In the latter case, the genome-based findings from the
corresponding patient sample are incorporated into the report and are
used as a primary source for expression profile prioritisation.

### WTS

The only required WTS input data are **read counts** provided in a
quantification file from either the `DRAGEN RNA` or
`bcbio-nextgen RNA-seq` pipeline.

#### DRAGEN RNA

The table below lists all input data accepted in `RNAsum`:

| Input file                           | Tool                                                                                                                                                 | Example                                                                                                | Required |
|--------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------|----------|
| Quantified transcript **abundances** | [salmon](https://salmon.readthedocs.io/en/latest/salmon.html) ([description](https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats)) | [TEST.quant.sf](/inst/rawdata/test_data/dragen/TEST.quant.sf)                                          | **Yes**  |
| **Fusion gene** list                 | [DRAGEN RNA](https://sapac.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html) | [TEST.fusion_candidates.final](/inst/rawdata/test_data/dragen/test_sample_WTS.fusion_candidates.final) | No       |

These files are expected to be organised in the following structure:

``` text
|
|____<SampleName>
  |____<SampleName>quant.sf
  |____<SampleName>.fusion_candidates.final
```

#### bcbio-nextgen (legacy)

The table below lists all input data accepted in `RNAsum`:

| Input file                           | Tool                                                                                                                            | Example                                                                                                                                                                                    | Required |
|--------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------|
| Quantified transcript **abundances** | [kallisto](https://pachterlab.github.io/kallisto/about) ([description](https://pachterlab.github.io/kallisto/starting#results)) | [abundance.tsv](/inst/rawdata/test_data/final/test_sample_WTS/kallisto/abundance.tsv), [Transcripts Per Million](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)         | **Yes**  |
| **Fusion gene** list                 | [arriba](https://arriba.readthedocs.io/en/latest), [pizzly](https://github.com/pmelsted/pizzly)                                 | [fusions.tsv](/inst/rawdata/test_data/final/test_sample_WTS/arriba/fusions.tsv), [test_sample_WTS-flat.tsv](/inst/rawdata/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv) | No       |
| **Fusion gene** plots                | [arriba plot](https://arriba.readthedocs.io/en/latest/visualization/)                                                           | [fusions.pdf](/inst/rawdata/test_data/final/test_sample_WTS/arriba/fusions.pdf)                                                                                                            | No       |

These files are expected to be organised in the following structure:

``` text
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

**Note**

`pizzly` outputs two [fusion gene](./fusions) files by default: -
*\<sample_name\>-flat.tsv* lists all (**unfiltered**) the gene fusions -
*\<sample_name\>-flat-filtered.tsv* lists **filtered** gene fusions

`RNAsum` makes use of gene fusions listed in the **unfiltered** file
since it was noted that some genuine fusions (based on WGS data and
curation efforts) get filtered out.

### WGS

`RNAsum` is designed to be compatible with WGS outputs generated from
`umccrise`.

The table below lists all input data accepted in `RNAsum`:

| Input file      | Tool                                                                    | Example                                                                                                                   | Required |
|-----------------|-------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------|----------|
| **SNVs/Indels** | [PCGR](https://github.com/sigven/pcgr)                                  | [pcgr.snvs_indels.tiers.tsv](/inst/rawdata/test_data/umccrised/test_sample_WGS/small_variants/pcgr.snvs_indels.tiers.tsv) | No       |
| **CNVs**        | [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) | [purple.cnv.gene.tsv](/inst/rawdata/test_data/umccrised/test_sample_WGS/purple/purple.gene.cnv)                           | No       |
| **SVs**         | [Manta](https://github.com/Illumina/manta)                              | [sv-prioritize-manta.tsv](/inst/rawdata/test_data/umccrised/test_sample_WGS/structural/sv-prioritize-manta.tsv)           | No       |

These files are expected to be organised in the following structure:

``` text
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

``` bash
rnasum_cli=$(Rscript -e 'x = system.file("cli", package = "RNAsum"); cat(x, "\n")' | xargs)
export PATH="${rnasum_cli}:${PATH}"
```

    $ rnasum.R --version
    rnasum.R x.x.x

    $ rnasum.R --help
    bash: line 8: rnasum.R: command not found

**Note**

Human reference genome
***[GRCh38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39)***
(*Ensembl* based annotation version ***86***) is used for gene
annotation by default. Alternatively, human reference genome
[GRCh37](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/)
(*Ensembl* based annotation version *75*) is used when argument
`grch_version` is set to `37`.

### Run via AWS-Batch

At UMCCR we run `RNAsum` in production via
[AWS-batch](https://aws.amazon.com/batch/). Before running, collate the
following information:

-   S3 path to umccrise results.
-   S3 path to WTS data results (bcbio).
-   [Reference
    dataset](https://github.com/umccr/RNAsum/blob/master/TCGA_projects_summary.md)
    to use for the sample.

Next, follow these steps:

1.  Configure the AWS CLI to access data on S3 (for the UMCCR account
    see our
    [wiki](https://github.com/umccr/wiki/blob/master/computing/cloud/amazon/aws_cli.md#configuration)).
2.  Run the `aws lambda invoke` command, based on the example below:

``` bash
aws lambda invoke \
  --region ap-southeast-2 \
  --function-name wts_report_trigger_lambda_prod \
  --cli-binary-format raw-in-base64-out \
  --payload '{"dataDirWGS":"Scott-SFRC/SBJ0000/WGS/2021-11-24/umccrised/SBJ0000__SB0000_something", "dataDirWTS":"Scott-SFRC/SBJ0000/WTS/2021-11-07/final/SBJ0000_something", "refDataset":"UCEC"}' \
  /tmp/lambda.output
```

The `RNAsum` HTML report for this sample will be created in the
`umccr-primary-data-prod` S3 bucket, under
`Scott-SFRC/SBJ0000/WTS/2021-11-07/RNAsum`. To use primary data from a
different S3 bucket, use `"dataBucket":"a-different-bucket"` in the
`aws lambda invoke` command. Also set `resultBucket` accordingly.

### Examples

Below are `RNAsum` CLI commands for generating HTML reports under
different data availability scenarios:

1.  [WTS data only](#1-wts-data-only)
2.  [WTS and WGS data](#2-wts-and-wgs-data)
3.  [WTS WGS and clinical data](#3-wts-wgs-and-clinical-data)

**Note**

-   Example data is provided in the </inst/rawdata/test_data> folder.
-   The `RNAsum` runtime should be less than **20 minutes** using **16GB
    RAM** memory and **1 CPU**.

#### 1. WTS data only

In this scenario, only [WTS](#wts) data will be used and only expression
levels of key
**[`Cancer genes`](https://github.com/umccr/umccrise/blob/master/workflow.md#key-cancer-genes)**,
**`Fusion genes`**, **`Immune markers`** and homologous recombination
deficiency genes (**`HRD genes`**) will be reported. Moreover, gene
fusions reported in the `Fusion genes` report section will not contain
information about evidence from genome-based data. A subset of the TCGA
pancreatic adenocarcinoma dataset is used as the reference cohort
(`--dataset TEST`).

##### DRAGEN RNA

``` bash
rnasum.R \
  --sample_name test_sample_WTS \
  --dataset TEST \
  --dragen_rnaseq $(pwd)/../rawdata/test_data/dragen \
  --report_dir $(pwd)/../rawdata/test_data/dragen/RNAsum \
  --save_tables FALSE
```

The HTML report `test_sample_WTS.RNAsum.html` will be created in the
`../rawdata/test_data/dragen/RNAsum` folder.

##### bcbio-nextgen (legacy)

``` bash
rnasum.R \
  --sample_name test_sample_WTS \
  --dataset TEST \
  --bcbio_rnaseq $(pwd)/../rawdata/test_data/final/test_sample_WTS \
  --report_dir $(pwd)/../rawdata/test_data/final/test_sample_WTS/RNAsum \
  --save_tables FALSE
```

The HTML report `test_sample_WTS.RNAsum.html` will be created in the
`../rawdata/test_data/final/test_sample_WTS/RNAsum` folder.

#### 2. WTS and WGS data

This is the **most frequent and preferred case**, in which the
[WGS](#wgs)-based findings will be used as a primary source for
expression profile prioritisation. The genome-based results can be
incorporated into the report by specifying the location of the
corresponding `umccrise` output files (including results from `PCGR`,
`PURPLE`, and `Manta`) using the `--umccrise` argument. The
**`Mutated genes`**, **`Structural variants`** and
**`CN altered genes`** report sections will contain information about
expression levels of the mutated genes, genes located within detected
SVs and CN altered regions, respectively. The results in the
**`Fusion genes`** section will be ordered based on the evidence from
genome-based data. A subset of the TCGA pancreatic adenocarcinoma
dataset is used as reference cohort (`--dataset TEST`).

##### DRAGEN RNA

``` bash
rnasum.R \
  --sample_name test_sample_WTS \
  --dataset TEST \
  --dragen_rnaseq $(pwd)/../rawdata/test_data/dragen \
  --report_dir $(pwd)/../rawdata/test_data/dragen/RNAsum \
  --umccrise $(pwd)/../rawdata/test_data/umccrised/test_sample_WGS \
  --save_tables FALSE
```

The HTML report `test_sample_WTS.RNAsum.html` will be created in the
`../rawdata/test_data/dragen/RNAsum` folder.

##### bcbio-nextgen (legacy)

``` bash
rnasum.R \
  --sample_name test_sample_WTS \
  --dataset TEST \
  --bcbio_rnaseq $(pwd)/../rawdata/test_data/final/test_sample_WTS \
  --report_dir $(pwd)/../rawdata/test_data/final/test_sample_WTS/RNAsum \
  --umccrise $(pwd)/../rawdata/test_data/umccrised/test_sample_WGS \
  --save_tables FALSE
```

The HTML report `test_sample_WTS.RNAsum.html` will be created in the
`../rawdata/test_data/final/test_sample_WTS/RNAsum` folder.

#### 3. WTS WGS and clinical data

For samples derived from subjects, for which clinical information is
available, a treatment regimen timeline can be added to the HTML report.
This can be added by specifying the location of a relevant excel
spreadsheet (see example
[test_clinical_data.xlsx](./data/test_data/test_clinical_data.xlsx))
using the `--clinical_info` argument. In this spreadsheet, at least one
of the following columns is expected: `NEOADJUVANT REGIMEN`,
`ADJUVANT REGIMEN`, `FIRST LINE REGIMEN`, `SECOND LINE REGIMEN` or
`THIRD LINE REGIMEN`, along with `START` and `STOP` dates of
corresponding treatments. A subset of the TCGA pancreatic adenocarcinoma
dataset is used as the reference cohort (`--dataset TEST`).

##### DRAGEN RNA

``` bash
rnasum.R \
  --sample_name test_sample_WTS \
  --dataset TEST \
  --dragen_rnaseq $(pwd)/../rawdata/test_data/dragen \
  --report_dir $(pwd)/../rawdata/test_data/dragen/RNAsum \
  --umccrise $(pwd)/../rawdata/test_data/umccrised/test_sample_WGS \
  --save_tables FALSE \
  --clinical_info $(pwd)/../rawdata/test_data/test_clinical_data.xlsx \
  --save_tables FALSE
```

The HTML report `test_sample_WTS.RNAsum.html` will be created in the
`../rawdata/test_data/stratus/test_sample_WTS_dragen_v3.9.3/RNAsum`
folder.

##### bcbio-nextgen (legacy)

``` bash
rnasum.R \
  --sample_name test_sample_WTS \
  --dataset TEST \
  --bcbio_rnaseq $(pwd)/../rawdata/test_data/final/test_sample_WTS \
  --report_dir $(pwd)/../rawdata/test_data/final/test_sample_WTS/RNAsum \
  --umccrise $(pwd)/../rawdata/test_data/umccrised/test_sample_WGS \
  --clinical_info $(pwd)/../rawdata/test_data/test_clinical_data.xlsx \
  --save_tables FALSE
```

The HTML report `test_sample_WTS.RNAsum.html` will be created in the
`../rawdata/test_data/final/test_sample_WTS/RNAsum` folder.

### Output

The pipeline generates a HTML ***Patient Transcriptome Summary***
**[report](#report)** and a [results](#results) folder:

    |
    |____<output>
      |____<SampleName>.<output>.html
      |____results
        |____exprTables
        |____glanceExprPlots
        |____...

#### Report

The generated HTML report includes searchable tables and interactive
plots presenting expression levels of altered genes, as well as links to
public resources describing the genes of interest. The report consists
of several sections, including:

-   Input data
-   Clinical information\*
-   Findings summary
-   Mutated genes\*\*
-   Fusion genes
-   Structural variants\*\*
-   CN altered genes\*\*
-   Immune markers
-   HRD genes
-   Cancer genes
-   Drug matching

\* if clinical information is available; see `--clinical_info` argument
<br /> \*\* if genome-based results are available; see `--umccrise`
argument

Detailed description of the report structure, including result
prioritisation and visualisation is available
[here](report_structure.md).

#### Results

The `results` folder contains intermediate files, including plots and
tables that are presented in the HTML report.
