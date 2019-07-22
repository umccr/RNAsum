## RNA-seq report workflow

The description of the main components involved in the **[read counts](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)** data processing, **[gene fusions](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv)** detection, integration with **WGS**-based data (processed using *[umccrise](https://github.com/umccr/umccrise)* pipeline) and results presentation in the *Transcriptome Patient Summary* **report**. 

## Table of contents

<!-- vim-markdown-toc GFM -->
* [Data pre-processing](#data-pre-processing)
    * [Read counts](#read-counts)
    * [Fusion genes](#fusion-genes)
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

The RNA-seq data is processed using *[bcbio-nextgen](https://github.com/bcbio/bcbio-nextgen)* *[RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* and includes per-gene **read counts** calculation with [kallisto](https://pachterlab.github.io/kallisto/about) and **gene fusions** detection wuth [pizzly](https://github.com/pmelsted/pizzly). 

### Read counts

The expected input read count data is the quantification *[abundance.tsv](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)* file from [kallisto](https://pachterlab.github.io/kallisto/about) (see example *[abundance.tsv](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)* file and its [description](https://pachterlab.github.io/kallisto/starting#results)). The per-gene abundances are reported in *estimated counts* (`est_counts`) and in *[Transcripts Per Million](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)* (`tpm`).


### Fusion genes

[Fusion genes](./fusions) are expected to be listed in the [flat table](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv) generated with [pizzly](https://github.com/pmelsted/pizzly). By default two output tables are provided: (1) *\<sample_name\>-flat.tsv* listing all gene fusions and (2) *\<sample_name\>-flat-filtered.tsv* listing only gene fusions remaining after filtering step. However, this workflow makes uses of gene fusion candidates listed in the **unfiltered** [pizzly](https://github.com/pmelsted/pizzly) output file (see example [test_sample_WTS-flat.tsv](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv)) since it was noted that some genuine fusions (based on WGS data and curation efforts) are excluded in the filtered [pizzly](https://github.com/pmelsted/pizzly) output file. If WGS results from **[umccrise](https://github.com/umccr/umccrise)** are available then fusion genes in the **`Fusion genes`** section are ordered based on the evidence from genome-based data. For more information about gene fusions and methods for their detectecion as visualisation can be found [here](./fusions/README.md).


## Counts post-processing


<img src="img/counts_post-processing_scheme.png" width="50%"> 

>Counts post-processing scheme.


## Data scaling

### Gene-wise Z-score transformation

<img src="img/Z-score_transformation_gene_wise.png" width="50%"> 

>Gene-wise Z-score transformation scheme.

### Group-wise centering

<img src="img/centering_group_wise.png" width="50%"> 

>Group-wise centering scheme.


## Integration with WGS-based results

work in progress..

## Report generation

work in progress..

The generated html-based ***Transcriptome Patient Summary*** **report** contains several sections described more in detail in [report_structure.md](report_structure.md):
