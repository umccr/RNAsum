## RNA-seq report workflow

The description of the main components involved in the **[read counts](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)** data processing, **[gene fusions](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv)** detection, integration with **WGS**-based data (processed using *[umccrise](https://github.com/umccr/umccrise)* pipeline) and results presentation in the *Transcriptome Patient Summary* report. 

## Table of contents

<!-- vim-markdown-toc GFM -->
* [Data pre-processing](#data-pre-processing)
    * [Read counts](#read-counts)
    * [Fusion genes detection](#fusion-genes-detection)
* [Counts post-processing](#counts-post-processing])
    * [Data collection](#data-collection)
    * [Transformation](#transformation)
    * [Filtering](#filtering)
    * [Normalisation](#normalisation)
    * [Combination](#combination)
    * [Batch-effects correction](#batch-effects-correction])
* [Data scaling](#data-scaling)
  * [Gene-wise Z-score transformation](#gene-wise-z-score-transformation)
  * [Group-wise centering](#group-wise-centering)
* [Integration with WGS-based results](#data-scaling)
* [Report generation](#report-generation)

<!-- vim-markdown-toc -->

## Data pre-processing

The RNA-seq data is processed using *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* including per-gene **[read counts](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)** calculation and **[gene fusions](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv)** detection. 


### Read counts


### Fusion genes detection

All [fusions](./fusions) listed in the unfiltered [pizzly](https://github.com/pmelsted/pizzly) output file (see example [test_sample_WTS-flat.tsv](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv)) are reported since it was noted that some genuine fusions (based on WGS data and curation efforts) are excluded in the filtered [pizzly](https://github.com/pmelsted/pizzly) output file. If WGS results from **[umccrise](https://github.com/umccr/umccrise)** are available then fusion genes in the **`Fusion genes`** section are ordered based on the evidence from genome-based data.


## Counts post-processing


<img src="img/counts_post-processing_scheme.pdf" width="80%"> 

>Counts post-processing scheme.


## Data scaling

### Gene-wise Z-score transformation

<img src="img/Z-score_transformation_gene_wise.pdf" width="80%"> 

>Gene-wise Z-score transformation scheme.

### Group-wise centering

<img src="img/centering_group_wise.pdf" width="80%"> 

>Group-wise centering scheme.


## Integration with WGS-based results

work in progress..

## Report generation

work in progress..

The generated html-based ***Transcriptome Patient Summary*** **report** contains several sections described more in detail in [report_structure.md](report_structure.md):
