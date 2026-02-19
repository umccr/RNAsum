# Input file formats and minimal content

This document summarises the format and minimal required content of input files used by RNAsum. For full report functionality, all listed columns are required.

## WTS (RNA-seq) inputs

### `--salmon` (or `--kallisto`)

**Salmon** (`quant.genes.sf` or `quant.sf`):
- Tab-separated file.
- For gene-level (`quant.genes.sf`): `Name`, `NumReads` columns required.
- For transcript-level (`quant.sf`): standard Salmon output, used with tx2gene mapping.

**Kallisto** (`abundance.tsv`):
- Tab-separated file with target_id and est_counts columns.

### `--arriba_tsv`

Arriba fusion detection output (`fusions.tsv`). Tab-separated; must include fusion gene columns (gene1, gene2).

### `--dragen_fusions`

DRAGEN RNA-seq `fusion_candidates.final` output. Tab-separated; first column `#FusionGene`, split by `--` into gene1 and gene2.

---

## WGS inputs

### `--pcgr_tiers_tsv`

PCGR variant calls (e.g. `snvs_indels.tiers.tsv` or `snv_indel_ann.tsv.gz`). Tab-separated.

**Required columns** (7 columns):

| Column | Alternative | Type | Purpose |
|--------|-------------|------|---------|
| `SYMBOL` | - | character | Gene symbol |
| `CONSEQUENCE` | - | character | Variant consequence (e.g. missense_variant) |
| `TIER` | `ACTIONABILITY_TIER` | character/integer | Tier (1–4); PCGR v2 uses ACTIONABILITY_TIER |
| `AF_TUMOR` | `VAF_TUMOR` | numeric | Allele fraction in tumour |
| `VARIANT_CLASS` | - | character | Variant class (e.g. SNV, indel) |
| `GENOMIC_CHANGE` | - | character | Genomic change notation |
| `PROTEIN_CHANGE` | - | character | Protein change notation |

**Minimal example**:

```tsv
SYMBOL	CONSEQUENCE	TIER	AF_TUMOR	VARIANT_CLASS	GENOMIC_CHANGE	PROTEIN_CHANGE
PIK3CA	missense_variant	2	0.215	SNV	3:g.179218303G>A	p.Glu545Lys
TP53	stop_gained	3	0.469	SNV	17:g.7674894G>A	p.Arg213Ter
```

### `--cn_gene_tsv`

Copy-number by gene file (e.g. `purple.cnv.gene.tsv`). Tab-separated.

**Required columns** (3 columns):

| Column | Type | Purpose |
|--------|------|---------|
| `gene` | character | Gene symbol |
| `minCopyNumber` | numeric | Minimum copy number |
| `maxCopyNumber` | numeric | Maximum copy number |

**Minimal example**:

```tsv
gene	minCopyNumber	maxCopyNumber
PIK3CA	2.5	2.5
TP53	1.0	1.0
```

### `--sv_tsv`

Structural variant annotations. Tab-separated. Two supported formats:

**Format 1** — Simple gene list:
- First column must be `Gene`.
- One gene symbol per row.

**Format 2** — Annotation column (e.g. umccrise SV output):
- Column `annotation` required.
- Each row: one or more comma-separated annotations.
- Each annotation: pipe-delimited `Event|Effect|Genes|Transcript|Detail|Tier`.
- The pipeline uses `Effect` (2nd field) and `Genes` (3rd field).

**Minimal example (Format 1)**:

```tsv
Gene
FAM43B
NTNG1
```

**Minimal example (Format 2)**:

```tsv
annotation
BND|no_func_effect|||unprioritized|4
DEL|frameshift_variant&start_lost|FAM43B|ENST00000332947_exon_1/1|unprioritized|4
```
