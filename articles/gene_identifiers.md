# Gene Identifier Handling

RNAsum uses **Ensembl gene IDs** as its internal identifier throughout
processing, then converts to **HGNC gene symbols** for display in the
final report. The reference genome is **GRCh38** with **Ensembl
annotation version 105**.

## Input identifiers

Patient expression data (Salmon or Kallisto) must use Ensembl gene or
transcript IDs:

- **`quant.genes.sf`** (Salmon gene-level): Ensembl gene IDs are read
  directly from the `Name` column. Version suffixes
  (e.g. `ENSG00000012048.21`) are stripped automatically.
- **`quant.sf`** (Salmon transcript-level) and **`abundance.tsv`**
  (Kallisto): Ensembl transcript IDs are collapsed to gene-level Ensembl
  IDs via [`tximport`](https://bioconductor.org/packages/tximport/),
  using a pre-built transcript-to-gene mapping derived from Ensembl v105
  via [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/)
  (stored in the companion `RNAsum.data` package).

PAR_Y pseudo-autosomal loci (genes duplicated on both X and Y
chromosomes) are filtered out at input to prevent row-name conflicts
when merging with reference data.

The TCGA reference cohorts bundled with RNAsum are also keyed on Ensembl
gene IDs, so sample and reference data join on a shared, stable
identifier space.

## Internal processing

All count matrices (sample + reference) are retained with Ensembl IDs
through normalisation, filtering, batch correction, and PCA.

## Annotation and symbol conversion

After processing, genes are annotated using the same Ensembl v105
mapping table (`tx_gene_id_105.rds`), which provides:

| Field | Description |
|----|----|
| `ENSEMBL` | Ensembl gene ID (e.g. `ENSG00000012048`) |
| `SYMBOL` | HGNC gene symbol (e.g. `BRCA1`) |
| `GENEBIOTYPE` | Ensembl biotype (e.g. `protein_coding`, `lncRNA`) |
| `SEQNAME` | Chromosome |
| `GENESEQSTART` / `GENESEQEND` | Genomic coordinates |

Two deduplication steps are applied:

1.  **Duplicate Ensembl IDs** — removed (keeps first occurrence).
2.  **Duplicate gene symbols** — removed (affects multi-locus genes such
    as Y_RNAs, SNORDs, and LINC genes that map to multiple genomic
    locations).

Genes with no HGNC symbol in the Ensembl v105 annotation are silently
excluded from the reportable set.

After annotation, the processed expression matrix is re-indexed on
**HGNC gene symbols** — this is the identifier used in all HTML report
sections (Mutated genes, Fusion genes, CN altered genes, Cancer genes,
etc.) and all output tables.

## Summary

    Ensembl transcript/gene IDs (input)
            ↓  tximport (transcript-level only)
    Ensembl gene IDs  ──────────────────────────── internal processing
            ↓  Ensembl v105 annotation (tx_gene_id_105.rds)
    HGNC gene symbols (report display and output tables)

The annotation reference (`tx_gene_id_105.rds`) is versioned and fixed —
it does not update automatically. GRCh37 is not supported.
