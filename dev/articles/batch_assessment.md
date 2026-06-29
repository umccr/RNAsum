# Batch Effect Assessment Guide

## Overview

Batch effects can significantly impact RNA-seq expression analysis when
comparing clinical samples to reference cohorts from different
protocols, platforms, or processing pipelines. RNAsum provides
comprehensive batch effect assessment tools with **three specialized
gene set approaches** to evaluate compatibility between your sample and
reference data before running the main analysis.

## Why Assess Batch Effects?

**Protocol Differences** between clinical and reference data can cause
systematic biases:

- **RNA extraction methods**: TRIzol vs column-based isolation
- **Library preparation**: Ribo-depletion vs poly-A selection (TCGA
  standard)
- **Sequencing platforms**: Different instruments or chemistry versions
- **Computational pipelines**: Varying alignment and quantification
  methods
- **Sample characteristics**: Fresh vs FFPE tissue, storage conditions

**Early Detection** of batch effects allows you to:

- Choose appropriate analysis parameters (`--batch_rm` flag)
- Select compatible reference cohorts
- Interpret results with proper context
- Validate findings with orthogonal methods

------------------------------------------------------------------------

## Three Gene Set Approaches

### 1. Top Variable Genes (Default)

**Best for**: Comprehensive assessment, research applications

Uses the most variable genes in the reference dataset for analysis,
providing robust statistical power and capturing global expression
patterns.

``` r

# Standard assessment with top variable genes
results <- assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "top_n",
  n_genes = 2000,  # Customizable (default: 2000)
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "batch_results_standard"
)

# Quick assessment with fewer genes
quick_batch_check(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "top_n",
  n_genes = 1000,
  save_plots = TRUE
)
```

**Advantages**: - Captures genome-wide expression patterns - Robust
statistical power with large gene sets - Protocol-agnostic assessment

**Use when**: - Protocol differences are unknown - General assessment
needed - Research contexts requiring comprehensive evaluation

### 2. Cancer Gene Sets

**Best for**: Clinical oncology samples, therapeutically focused
analysis

Uses curated cancer gene databases focusing on clinically relevant and
actionable genes. **Enhanced with comprehensive gene mapping** achieving
up to 125+ cancer genes from major databases.

``` r

# UMCCR cancer genes (comprehensive mapping)
results_umccr <- assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "cancer_genes",
  cancer_gene_source = "umccr",  # 1,248 total genes → ~125 mapped
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "batch_results_cancer_umccr"
)

# OncoKB cancer genes (clinically actionable)
results_oncokb <- assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "cancer_genes",
  cancer_gene_source = "oncokb",  # 1,019 total genes → ~125 mapped
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "batch_results_cancer_oncokb"
)

# Combined cancer gene databases (maximum coverage)
results_combined <- assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "cancer_genes",
  cancer_gene_source = "combined",  # 1,315 unique genes → ~125 mapped
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "batch_results_cancer_combined"
)
```

**Available Cancer Gene Databases:**

- **UMCCR database** (1,248 genes): Comprehensive collection from COSMIC
  CGC, OncoKB, Cancermine, and clinical panels
- **OncoKB database** (1,019 genes): Clinically actionable genes with
  therapeutic relevance
- **Combined databases** (1,315 genes): Union of UMCCR and OncoKB for
  maximum coverage

**Enhanced Gene Mapping**: Automatically maps gene symbols to Ensembl
IDs with comprehensive lookup table covering 150+ cancer genes,
achieving ~125 successfully mapped genes per analysis.

**Advantages**: - Clinically relevant assessment focusing on actionable
genes - Therapeutically informed batch effect evaluation - Multiple
database options for different clinical contexts

**Use when**: - Analyzing cancer samples - Clinical diagnostic
contexts - Therapeutic decision support workflows

### 3. Custom Gene Sets

**Best for**: Pathway-specific analysis, hypothesis-driven research

Allows specification of custom gene sets tailored to specific research
questions, pathways, or clinical contexts.

``` r

# Example: DNA repair pathway genes
repair_genes <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000139618",
                  "ENSG00000083093", "ENSG00000149311", "ENSG00000143627")
                  # TP53, BRCA1, BRCA2, PALB2, ATM, ATR

results_custom <- assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "custom",
  gene_subset = repair_genes,
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "batch_results_custom"
)

# Example: Immunotherapy targets
immuno_genes <- c("ENSG00000188389", "ENSG00000120217", "ENSG00000163599")
                  # PDCD1 (PD-1), CD274 (PD-L1), CTLA4

results_immuno <- assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "custom",
  gene_subset = immuno_genes,
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "batch_results_immunotherapy"
)
```

**Advantages**: - Targeted assessment for specific biological
questions - User-controlled gene selection - Pathway or
mechanism-specific evaluation

**Use when**: - Specific pathways are of primary interest -
Hypothesis-driven research contexts - Focused clinical applications
(e.g., specific drug targets)

------------------------------------------------------------------------

## Understanding Diagnostic Plots and Statistics

Batch effect assessment generates **four key diagnostic plots** with
accompanying statistics to evaluate protocol compatibility.

### 1. PCA Analysis Plot (`*_pca_analysis.png`)

Shows Principal Component Analysis with clinical sample (red point)
positioned relative to reference cohort samples (blue points).

**Purpose**: Visualizes overall expression profile similarity

**Interpretation**: - **Clinical sample within blue cloud**: Good
protocol compatibility - **Clinical sample outside blue cloud**:
Potential batch effects - **Distance from centroid**: Quantified as PCA
distance percentile

**Key Metric**: PCA distance percentile - `>95%`: 🔴 **High batch effect
risk** - `--batch_rm` strongly recommended - `80-95%`: 🟡 **Moderate
risk** - Consider `--batch_rm` based on use case - `<80%`: 🟢
**Generally acceptable** - TCGA comparison acceptable

### 2. Expression Distribution Comparison (`*_expression_comparison.png`)

Side-by-side density plots comparing log2(TPM+1) expression
distributions between clinical sample and reference cohort.

**Purpose**: Identifies systematic expression level shifts

**Interpretation**: - **Overlapping curves**: Similar expression
ranges - **Shifted curves**: Protocol-induced expression differences -
**Different shapes**: Distinct expression characteristics

**Key Metrics**: - Expression level shift (log2 units) - Variance ratio
between sample and reference

### 3. Z-score Distribution Analysis (`*_zscore_distribution.png`)

Histogram showing Z-score distribution of sample expression relative to
reference cohort, with ±2SD and ±3SD threshold lines.

**Purpose**: Evaluates expression ranking accuracy

**Interpretation**: - **Normal bell curve around zero**: Good ranking
preservation - **Skewed or shifted distribution**: Systematic ranking
errors - **High proportions beyond ±2SD**: Expression ranking issues

**Key Metrics**: - Mean \|Z-score\|: Overall deviation magnitude -
Extreme Z-scores (\>±2SD): Percentage of mis-ranked genes - Outlier
genes (\>±3SD): Severely affected genes

### 4. Reference Correlation Heatmap (`*_correlation_heatmap.png`)

Heatmap showing Pearson correlations between clinical sample and
individual reference samples, arranged by correlation strength.

**Purpose**: Identifies compatible reference samples and overall cohort
fit

**Interpretation**: - **High correlations (\>0.8)**: Strong protocol
compatibility - **Moderate correlations (0.5-0.8)**: Acceptable with
possible batch correction - **Low correlations (\<0.5)**: Substantial
protocol differences

**Key Metrics**: - Median Pearson correlation: Overall cohort
compatibility - Median Spearman correlation: Ranking-based
compatibility - Correlation range and variability

------------------------------------------------------------------------

## Statistical Interpretation Guidelines

### Batch Effect Risk Assessment Matrix

| PCA Distance | Median Correlation | Z-score Extremes | Risk Level | Recommendation |
|----|----|----|----|----|
| \>95% | \<0.5 | \>10% | 🔴 **High** | `--batch_rm` strongly recommended |
| 80-95% | 0.5-0.7 | 5-10% | 🟡 **Moderate** | Consider `--batch_rm` based on use case |
| \<80% | \>0.7 | \<5% | 🟢 **Low** | TCGA comparison acceptable |

### Clinical Decision Framework

**High Risk (🔴)** - Multiple warning signs:

``` bash
rnasum --sample_name clinical_sample --batch_rm --dataset BRCA [other options...]
```

**Moderate Risk (🟡)** - Context-dependent decision: - **Clinical
diagnostics**: Use `--batch_rm` for safety - **Research applications**:
Validate findings with orthogonal methods - **Exploratory analysis**:
TCGA comparison may be acceptable

**Low Risk (🟢)** - Standard workflow:

``` bash
rnasum --sample_name clinical_sample --dataset BRCA [other options...]
```

------------------------------------------------------------------------

## Practical Implementation

### Basic Setup

The batch assessment functions are exported by RNAsum, so after
installation you can call them directly:

``` r

library(RNAsum)
# assess_batch_effects() and quick_batch_check() are now available
```

> **Note**: If working from a local source clone without installation,
> you may need: `source("R/batch_assessment.R")` from the repository
> root.

### Complete Workflow Example

``` r

library(RNAsum)
library(dplyr)

# Load your sample data
sample_file <- "path/to/your/sample.quant.genes.sf"
sample_data <- read.delim(sample_file)
sample_tpm <- setNames(sample_data$TPM, sample_data$Name)

# Load reference data
ref_paths <- get_refdata(dataset = "PANCAN", batch_rm = FALSE)
counts_file <- ref_paths$ext_ref$counts
counts_data <- utils::read.table(gzfile(counts_file), header = TRUE, sep = "\t", row.names = NULL)
ref_matrix <- as.matrix(counts_data[, -1])
rownames(ref_matrix) <- counts_data[[1]]

# Handle gene ID versioning (if needed)
names(sample_tpm) <- gsub("\\..*$", "", names(sample_tpm))  # Remove versions
rownames(ref_matrix) <- gsub("\\..*$", "", rownames(ref_matrix))

# Remove duplicates and get common genes
sample_tpm <- sample_tpm[!duplicated(names(sample_tpm))]
ref_matrix <- ref_matrix[!duplicated(rownames(ref_matrix)), ]
common_genes <- intersect(names(sample_tpm), rownames(ref_matrix))
sample_final <- sample_tpm[common_genes]
ref_final <- ref_matrix[common_genes, ]

# Run comprehensive assessment with all three approaches

# 1. Top variable genes (comprehensive)
results_comprehensive <- assess_batch_effects(
  sample_data = sample_final,
  reference_data = ref_final,
  gene_set_type = "top_n",
  n_genes = 2000,
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "batch_results_comprehensive"
)

# 2. Enhanced cancer genes (clinical focus)
results_cancer <- assess_batch_effects(
  sample_data = sample_final,
  reference_data = ref_final,
  gene_set_type = "cancer_genes",
  cancer_gene_source = "combined",  # Maximum coverage
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "batch_results_cancer"
)

# 3. Custom pathway-specific genes
custom_genes <- c("ENSG00000141510", "ENSG00000133703", "ENSG00000146648")  # TP53, KRAS, EGFR
results_custom <- assess_batch_effects(
  sample_data = sample_final,
  reference_data = ref_final,
  gene_set_type = "custom",
  gene_subset = custom_genes,
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "batch_results_custom"
)

# Quick assessment for rapid evaluation
quick_batch_check(
  sample_data = sample_final,
  reference_data = ref_final,
  gene_set_type = "cancer_genes",
  cancer_gene_source = "umccr",
  save_plots = TRUE
)
```

### Gene ID Format Considerations

**Automatic Format Detection**: The enhanced system automatically
detects gene ID formats and provides helpful error messages.

**Common Formats**: - **Sample data**: Typically uses Ensembl gene IDs
(e.g., “ENSG00000141510”) - **Cancer databases**: Use gene symbols
(e.g., “TP53”, “BRCA1”) - **Custom sets**: Must match sample data format
for proper analysis

**Gene ID Versioning**: Some datasets include Ensembl version numbers
(e.g., “ENSG00000141510.14”). The system automatically handles version
compatibility:

``` r

# Automatic version normalization (if needed)
names(sample_tpm) <- gsub("\\..*$", "", names(sample_tpm))
rownames(ref_matrix) <- gsub("\\..*$", "", rownames(ref_matrix))
```

**Format Mismatch Handling**: When formats don’t match, you’ll receive
clear guidance:

    No genes from custom subset found in data (0/31). This may be due to gene ID format mismatch:
      • Sample data format: Ensembl gene IDs (with versions)
      • Custom genes format: Gene symbols (HGNC-like)
      • Suggestions:
        - Use gene_set_type='cancer_genes' for cancer-related analysis
        - Ensure custom genes match sample data ID format

------------------------------------------------------------------------

## Example with TEST Data

Using the provided TEST sample for demonstration:

``` r

library(RNAsum)
library(dplyr)

# Load TEST sample data
test_file <- system.file("rawdata/test_data/dragen/TEST.quant.genes.sf", package = "RNAsum")
sample_data <- read.delim(test_file)
sample_tpm <- setNames(sample_data$TPM, sample_data$Name)

# Get TEST reference data
ref_paths <- get_refdata(dataset = "TEST", batch_rm = FALSE)
ref_counts <- utils::read.table(gzfile(ref_paths$ext_ref$counts),
                                header = TRUE, sep = "\t", row.names = NULL)
ref_matrix <- as.matrix(ref_counts[, -1])
rownames(ref_matrix) <- ref_counts[[1]]

# Run cancer gene assessment
batch_results <- assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "cancer_genes",
  cancer_gene_source = "umccr",
  protocol_clinical = "ribo-depletion",
  protocol_reference = "TCGA_poly-A",
  output_dir = "test_batch_assessment"
)
```

**Expected Results with TEST Data**: - **PCA distance percentile**:
~50-70% (moderate positioning) - **Median correlation**: ~0.6-0.8
(reasonable compatibility) - **Z-score extremes**: ~3-7% (acceptable
ranking preservation) - **Verdict**: Likely 🟡 **Moderate risk** -
consider `--batch_rm` for clinical applications

------------------------------------------------------------------------

## Advanced Considerations

### Gene Set Selection Strategy

1.  **Start with cancer genes** for clinical samples to assess
    clinically relevant batch effects
2.  **Use top variable genes** for comprehensive validation and research
    contexts
3.  **Apply custom genes** for pathway-specific or targeted analysis
    when specific gene sets are critical

### Protocol-Specific Considerations

**Major Batch Effect Sources**: - **Ribo-depletion vs poly-A
selection**: Most significant source of systematic bias - **Different
sequencing platforms**: Can introduce platform-specific artifacts -
**RNA quality differences**: Fresh vs FFPE tissue shows distinct
patterns - **Processing time gaps**: Temporal batch effects from
different processing periods

### Validation Approaches

- **Visual inspection**: PCA plots showing sample clustering patterns
- **Correlation analysis**: Expression correlation patterns with
  reference cohorts
- **Known controls**: Validate expression patterns for
  well-characterized genes
- **Clinical relevance**: Ensure biologically plausible expression
  rankings
- **Multiple approaches**: Cross-validate findings across different gene
  set approaches

### Integration with RNAsum Workflow

**Batch Effect Detection** → **Parameter Selection** → **Analysis
Execution**

``` bash
# Based on batch assessment results:

# High batch effects detected
rnasum --sample_name clinical_sample --batch_rm --dataset BRCA [other options...]

# Low batch effects
rnasum --sample_name clinical_sample --dataset BRCA [other options...]
```

------------------------------------------------------------------------

## Enhanced Features

### Comprehensive Cancer Gene Mapping

- **Enhanced coverage**: Up to 125+ cancer genes successfully mapped
  vs. previous 63-gene limitation
- **Multiple databases**: UMCCR, OncoKB, and combined options
- **Intelligent mapping**: Automatic gene symbol to Ensembl ID
  conversion with robust error handling
- **Version compatibility**: Handles Ensembl ID versioning automatically

### Improved Error Handling

- **Format detection**: Automatic gene ID format recognition
- **Clear error messages**: Specific guidance for resolving gene ID
  mismatches
- **Fallback mechanisms**: Multiple approaches for gene mapping
- **Diagnostic reporting**: Detailed logs for troubleshooting

### Clinical Decision Support

- **Risk stratification**: Three-tier risk assessment
  (High/Moderate/Low)
- **Parameter recommendations**: Specific `--batch_rm` guidance
- **Context-aware advice**: Different recommendations for clinical vs
  research use
- **Quality metrics**: Comprehensive statistical evaluation framework

------------------------------------------------------------------------

## Summary

The enhanced batch effect assessment provides robust evaluation of
protocol compatibility with three complementary approaches:

✅ **Top Variable Genes**: Comprehensive genome-wide assessment ✅
**Cancer Gene Sets**: Clinically focused evaluation with enhanced
mapping ✅ **Custom Gene Sets**: Targeted pathway-specific analysis

Combined with detailed diagnostic plots, statistical interpretation
guidelines, and clinical decision frameworks, this provides a complete
solution for batch effect evaluation in RNA-seq analysis workflows.
