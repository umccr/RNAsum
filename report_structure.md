## RNA-seq report sections

<!-- vim-markdown-toc GFM -->
* [Input data](#input-data)
* [Clinical information](#clinical-information)
* [Findings summary](#findings-summary)
* [Mutated genes](#mutated-genes)
* [Fusion genes](#fusion-genes)
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


### Structural variants

### CN altered genes

### Immune markers

### HRD genes

### Cancer genes

### Drug matching

### Addendum


