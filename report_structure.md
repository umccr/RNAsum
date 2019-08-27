## RNA-seq report sections

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



The **`Mutated genes`**, **`Structural variants`** and **`CN altered genes`** sections will contain information about expression levels of the mutated genes, genes located within detected structural variants (SVs) and copy-number (CN) altered regions, respectively. Genes will be ordered by increasing *variants* `TIER`, *SV* `score` and `CN` *value*, resepctively, and then by decreasing absolute values in the `Patient` vs selected `dataset` column. Moreover, gene fusions detected in WTS data and reported in **`Fusion genes`** section will be first ordered based on the evidence from genome-based data (`DNA support (gene A/B)` columns).

### Input data

Summary of the input data, including 

### Clinical information

Treatment regimen information for this patient is AVAILABLE.

NOTE: for confidentiality reasons, the timeline (x-axis) projecting patient’s treatment regimens (y-axis) is set to start from 1st January 2000, but the treatments lengths are preserved.

### Findings summary

### Mutated genes

### Fusion genes

If WGS results from **[umccrise](https://github.com/umccr/umccrise)** are available then fusion genes in the **`Fusion genes`** report section are ordered based on the evidence from genome-based data. For more information about gene fusions and methods for their detectecion and visualisation can be found [here](./fusions/README.md).

#### Prioritisation

Fusion genes detected in transcriptome data are prioritised based on criteria ranked in the following order:

1. Involvement of fusion gene(s) **detected in genomic data** (if [Structural variants](#structural-variants) results are available)
2. Listed as reported fusion event according to [FusionGDB](https://ccsm.uth.edu/FusionGDB/) database
3. Involvement of highly **abundant transcript(s)** (see [abundant transcripts definition](#abundant-transcripts) section)
4. Decreasing number of split counts
5. Decreasing number of pair counts
6. Involvement of **cancer gene(s)** (see [Cancer genes](#cancer-genes) section)

#### Filtering

Fusion genes detected in transcriptome data are reported if **at least one** of the following criteria is met:

1. Involvement of fusion gene(s) **detected in genomic data** (if [Structural variants](#structural-variants) results are available)
2. **Reported** fusion event according to [FusionGDB](https://ccsm.uth.edu/FusionGDB) database
3. Involvement of highly **abundant transcript(s)**  (see [abundant transcripts definition](#abundant-transcripts) section)
4. Involvement of **cancer gene(s)** (see [Cancer genes](#cancer-genes) section)
5. **Split counts** > 1
6. **Pair counts** > 1 and **split counts** > 1

#### Abundant transcripts

The following steps were performed to define abundant transcripts involved in detected fusion events:

1. Run [kallisto](https://github.com/pachterlab/kallisto) to quantify the fusion transcripts reported by [pizzly](https://github.com/pmelsted/pizzly) and select those which are supported by decent number of [Transcripts Per Kilobase Million](http://www.arrayserver.com/wiki/index.php?title=TPM) (TPM)
2. Create a new index based on the transcriptome and the fusion transcripts identified by [pizzly](https://github.com/pmelsted/pizzly)
3. Run [kallisto](https://github.com/pachterlab/kallisto) in normal quantification mode on the expanded index to quantify both normal transcripts and fusions
4. Select fusion genes involving transcripts with [TPM](http://www.arrayserver.com/wiki/index.php?title=TPM) values above 90th percentile of all [TPM](http://www.arrayserver.com/wiki/index.php?title=TPM) values (as reported by previous step)

### Structural variants

### CN altered genes

Section overlaying the mRNA expression data for [cancer genes](#cancer-genes) with per-gene somatic copy-number (CN) data (from [PURPLE](https://anaconda.org/bioconda/hmftools-purple)) and mutation status, if available.

### Immune markers

### HRD genes

### Cancer genes

### Drug matching

### Addendum


