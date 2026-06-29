# Generate table with coloured cells indicating expression values for selected genes

Generates table with coloured cells indicating expression values for
selected genes.

## Usage

``` r
exprTable(
  data = NULL,
  genes = NULL,
  keep_all = FALSE,
  cn_data = NULL,
  sv_data = NULL,
  cn_decrease = TRUE,
  targets,
  sampleName,
  int_cancer,
  ext_cancer,
  comp_cancer,
  add_cancer = NULL,
  genes_annot = NULL,
  oncokb_annot = NULL,
  cancer_genes = NULL,
  mut_annot = NULL,
  fusion_genes = NULL,
  type = "z",
  scaling = "gene-wise",
  civic_clin_evid = NULL,
  batch_rm = FALSE
)
```

## Arguments

- data:

  Input data.

- genes:

  Selected genes.

- keep_all:

  keep all rows

- cn_data:

  Copy number data.

- sv_data:

  SV data.

- cn_decrease:

  Order of the CN data.

- targets:

  Target groups.

- sampleName:

  Sample name.

- int_cancer:

  Internal cancer group.

- ext_cancer:

  External cancer group.

- comp_cancer:

  Complete cancer group.

- add_cancer:

  Additonal cancer type.

- genes_annot:

  Geners annotation.

- oncokb_annot:

  OncoKb annotation.

- cancer_genes:

  Cancer genes.

- mut_annot:

  Mutation annotation.

- fusion_genes:

  Fusion genes.

- type:

  Type.

- scaling:

  Scaling

- civic_clin_evid:

  civic_clin_evid tibble from reference genes.

- batch_rm:

  Remove batch-associated effects between datasets.

## Value

Table with coloured cells indicating expression values for selected
genes
