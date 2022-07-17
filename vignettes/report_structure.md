## RNAsum sections

<!-- vim-markdown-toc GFM -->
* [Input data](#input-data)
* [Clinical information](#clinical-information)
* [Findings summary](#findings-summary)
* [Mutated genes](#mutated-genes)
* [Fusion genes](#fusion-genes)
  * [Prioritisation](#prioritisation)
  * [Filtering](#filtering)
  * [Abundant transcripts](#abundant-transcripts)
* [Structural variants](#structural-variants)
* [CN altered genes](#cn-altered-genes)
* [Immune markers](#immune-markers)
* [HRD genes](#hrd-genes)
* [Cancer genes](#cancer-genes)
* [Drug matching](#drug-matching)
* [Addendum](#addendum)

<!-- vim-markdown-toc -->

<br/> 

The **`Mutated genes`**, **`Structural variants`** and **`CN altered genes`** sections will contain information about expression levels of the mutated genes, genes located within detected structural variants (SVs) and copy-number (CN) altered regions, respectively. Genes will be ordered by increasing *variants* `TIER`, *SV* `score` and `CN` *value*, resepctively, and then by decreasing absolute values in the `Patient` vs selected `dataset` column. Moreover, gene fusions detected in WTS data and reported in **`Fusion genes`** section will be first ordered based on the evidence from genome-based data (`DNA support (gene A/B)` columns).

***

### Input data

Summary of the input data

***

### Clinical information

Treatment regimen information for patient for which clinical information is available.

NOTE: for confidentiality reasons, the timeline (x-axis) projecting patient’s treatment regimens (y-axis) is set to start from 1st January 2000, but the treatments lengths are preserved.

***

### Findings summary

Plot and table summarising altered genes listed across various report sections

***

### Mutated genes

mRNA expression levels of mutated genes (containing single nucleotide variants (SNVs) or insertions/deletions (indels)) measured in patient's sample and their average mRNA expression in samples from cancer patients (from [TCGA](https://portal.gdc.cancer.gov/)). This section is available only for samples with available *[umccrise](https://github.com/umccr/umccrise) results*

***

### Fusion genes

Prioritised fusion genes based on [Arriba](https://arriba.readthedocs.io/en/latest/) results and annotated with [FusionGDB](https://ccsm.uth.edu/FusionGDB) database. If WGS results from **[umccrise](https://github.com/umccr/umccrise)** are available then fusion genes in the **`Fusion genes`** report section are ordered based on the evidence from genome-based data. For more information about gene fusions and methods for their detectecion and visualisation can be found [here](./fusions/README.md).

#### Prioritisation

Fusion genes detected in transcriptome data are prioritised based on criteria ranked in the following order:

1. Involvement of fusion gene(s) **detected in genomic data** (if [Structural variants](#structural-variants) results are available)
2. **Detected in transcriptome data** by [Arriba](https://arriba.readthedocs.io/en/latest/) tool
3. **Reported** fusion event according to [FusionGDB](https://ccsm.uth.edu/FusionGDB/) database
4. Decreasing number of **split reads**
5. Decreasing number of **pair reads**
6. Involvement of **cancer gene(s)** (see [Cancer genes](#cancer-genes) section)

#### Filtering

Fusion genes detected in transcriptome data are reported if **at least one** of the following criteria is met:

1. Involvement of fusion gene(s) **detected in genomic data** (if [Structural variants](#structural-variants) results are available)
2. **Reported** fusion event according to [FusionGDB](https://ccsm.uth.edu/FusionGDB) database
3. Involvement of **cancer gene(s)** (see [Cancer genes](#cancer-genes) section)
4. **Split reads** > 1
5. **Pair reads** > 1 and **split reads** > 1

***

### Structural variants

Similar to *Mutated genes* analysis but limited to genes located within structural variants (SVs) detected by [MANTA](https://github.com/Illumina/manta) using genomic data. This section is available only for samples with available *[MANTA](https://github.com/Illumina/manta) results*.

***

### CN altered genes

Section overlaying the mRNA expression data for [cancer genes](#cancer-genes) with per-gene somatic copy-number (CN) data (from [PURPLE](https://anaconda.org/bioconda/hmftools-purple)) and mutation status, if available.

***

### Immune markers

Similar to *Mutated genes* analysis but limited to genes considered to be immune markers. The immune markers used in the report are listed in PanelApp panel [Immune markers for WTS report](https://panelapp.agha.umccr.org/panels/243/).

***

### HRD genes

Similar to *Mutated genes* analysis but limited to genes considered to be homologous recombination deficiency (HRD) genes. The HRD genes used in the report are listed in PanelApp panel [Homologous recombination deficiency (HDR) for WTS report](https://panelapp.agha.umccr.org/panels/242/).

***

### Cancer genes

Similar to analysis above, but limited to *UMCCR cancer genes*.

***

### Drug matching

List of drugs targeting variants in detected *mutated genes*, *fusion genes*, *structural variants-affected genes*, *CN altered genes*, *HRD genes* and dysregulated *cancer genes*, which can be considered in the treatment decision making process.

###### Note

This section is not displayed as default. Set the `--drugs` argument to `TRUE` to present it in the report.

***

### Addendum

Additional information, including `Parameters`, `Reporter details` and R `Session information`,  added at the end of the report.

<br/>

