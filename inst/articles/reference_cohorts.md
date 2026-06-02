## Reference cohorts

RNAsum provides high-quality reference expression data for **33 cancer types** from two sources:

- **External cohorts**: [TCGA](https://tcga-data.nci.nih.gov/) cancer datasets
- **Internal cohort**: [CCGCM](https://genomic-cancer-medicine.unimelb.edu.au/research/clinical-cancer-genomics) pancreatic cancer samples

These reference datasets enable contextualised interpretation of patient expression profiles.

---

### External reference cohorts (TCGA)

We have curated high-quality cancer reference cohorts from The Cancer Genome Atlas (TCGA) to provide expression context for patient samples.

#### Available Datasets

**33 cancer-specific datasets** are available, each representing a distinct tumour type. Depending on the tissue source of your patient sample, select the appropriate TCGA cohort for comparison.

**Special datasets:**

- **Pan-Cancer dataset**: 330 samples (10 from each of the 33 TCGA datasets) for cross-cancer comparisons
- **Extended sets**: Available for select cohorts with larger sample sizes

**Usage:** Specify the desired dataset using the `--dataset` argument with a TCGA project ID.

#### Dataset details

| Feature | Description |
|---------|-------------|
| **Data source** | TCGA (GDC harmonized data) |
| **Processing** | Standardized pipeline ([see details](https://github.com/umccr/TCGA-data-harmonization/blob/master/expression/README.md#gdc-counts-data)) |
| **Reference** | [TCGA projects summary table](./TCGA_projects_summary.md) |

**Processing pipeline:**
All TCGA datasets have been harmonized using methods documented in the [TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/blob/master/expression/README.md#gdc-counts-data) repository, ensuring consistency and reproducibility.

---

### Internal reference cohort

To minimize batch effects and improve comparability, internal reference cohorts processed with identical pipelines to the input data are highly recommended, especially when substantial protocol differences exist.

#### Why an internal cohort?

Public TCGA datasets may show significant [batch effects](https://www.ncbi.nlm.nih.gov/pubmed/20838408) when compared to in-house data due to:

- **RNA extraction protocols**: Different RNA isolation methods (e.g., TRIzol vs column-based)
- **Library preparation**: Ribo-depletion vs poly-A selection (TCGA standard)
- **Sequencing platforms**: Different instruments or chemistry versions
- **Strandedness**: Protocol-specific strand preservation differences
- **Computational pipelines**: Varying alignment and quantification methods
- **Sample characteristics**: Fresh vs FFPE tissue, storage conditions, tumor purity

#### Batch Harmonization Strategies

**1. Protocol Matching**

- Match library preparation methods when possible
- Document and account for unavoidable protocol differences
- Use strand-specific analysis consistently across cohorts

**2. Computational Approaches**

- Enable `--batch_rm` for systematic batch correction
- Apply consistent normalization methods across datasets
- Use robust statistical approaches (median-based rather than mean-based)

**3. Quality Assessment**

- Evaluate batch correction effectiveness using PCA analysis
- Compare expression rankings before and after correction
- Validate using known positive controls when available

#### Key Advantages

✓ **Minimized batch effects** - same protocols as your input data  
✓ **High tissue quality** - stringent QC standards  
✓ **Matched processing** - [data pre-processing](./workflow.md#data-processing) identical to patient samples  
✓ **Protocol consistency** - eliminates major sources of technical variation  
✓ **Validation capability** - enables assessment of batch correction methods  

---

### Choosing the right reference

#### General Guidelines

- Use the matching TCGA **cancer-specific** dataset for tissue-appropriate comparisons
- Consider **Pan-Cancer** for cross-cancer comparisons or rare tumour types
- Use the **internal cohort** for maximum comparability and minimal batch effects

#### Batch Effect Considerations

**High batch effect risk scenarios** (prioritize internal cohorts + `--batch_rm`):

🔴 **Protocol mismatches**:

- Clinical: Ribo-depletion ↔ TCGA: Poly-A selection
- Clinical: Strand-specific ↔ TCGA: Unstranded analysis
- Clinical: Fresh tissue ↔ TCGA: Variable preservation methods

🔴 **Platform differences**:

- Different sequencing platforms (e.g., NovaSeq vs HiSeq)
- Significant time gaps between data generation
- Different computational pipelines

**Moderate batch effect risk** (consider `--batch_rm` based on assessment):

🟡 **Sample characteristics**:

- Tumor purity differences
- RNA quality variations
- Processing time differences

**Low batch effect risk** (TCGA comparison acceptable):

🟢 **Well-matched scenarios**:

- Similar protocols and platforms
- High-quality samples with good RNA integrity
- Recent TCGA processing with updated methods

#### Decision Framework

1. **Assess protocol compatibility** between clinical and reference data
2. **Evaluate sample quality** metrics and characteristics  
3. **Consider using batch assessment tools** (see [batch effects workflow](./workflow.md#batch-effects-correction-optional))
4. **Choose appropriate combination**:
   - **Internal + TCGA with `--batch_rm`**: Maximum robustness for high-risk scenarios
   - **TCGA only with `--batch_rm`**: Moderate-risk scenarios with protocol awareness
   - **TCGA only**: Low-risk scenarios with well-matched protocols

#### Validation Approaches

- **Visual inspection**: PCA plots showing sample clustering patterns
- **Correlation analysis**: Expression correlation between sample and reference cohorts
- **Known controls**: Validate expression patterns for well-characterized genes
- **Clinical relevance**: Ensure biologically plausible expression rankings

#### Enhanced Batch Assessment with Gene Set Selection

RNAsum provides **three specialized gene set options** for targeted batch effect assessment, allowing you to focus the analysis on the most relevant genes for your application.

> **Usage note:** `assess_batch_effects()` and `quick_batch_check()` are exported by RNAsum. After `library(RNAsum)` you can call them directly — no `source()` step is required. If you are working from a local clone without installing the package, you may instead `source("R/batch_assessment.R")` from the repo root.

**1. Top Variable Genes (Default)**

```r
assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "top_n",
  n_genes = 2000  # Customizable
)
```
- **Best for**: Comprehensive assessment, research applications
- **Advantages**: Captures global expression patterns, robust statistical power
- **Use when**: Protocol differences are unknown or general assessment needed

**2. Cancer Gene Sets**

```r
assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "cancer_genes",
  cancer_gene_source = "umccr"  # Options: "umccr", "oncokb", "combined"
)
```

**Available cancer gene databases:**

- **UMCCR database** (1,248 genes): Comprehensive cancer gene collection from multiple sources including COSMIC CGC, OncoKB, Cancermine, and clinical panels
- **OncoKB database** (1,019 genes): Clinically actionable cancer genes with therapeutic relevance  
- **Combined databases** (1,315 genes): Union of UMCCR and OncoKB for maximum coverage

- **Best for**: Clinical oncology samples, cancer-focused research
- **Advantages**: Clinically relevant assessment, focused on actionable genes
- **Use when**: Cancer samples require assessment of therapeutically relevant genes

**3. Custom Gene Sets**

```r
# Gene IDs must match sample data format
custom_genes <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000139618")
assess_batch_effects(
  sample_data = sample_tpm,
  reference_data = ref_matrix,
  gene_set_type = "custom",
  gene_subset = custom_genes
)
```
- **Best for**: Pathway-specific analysis, hypothesis-driven research
- **Advantages**: Targeted assessment, user-controlled gene selection
- **Use when**: Specific pathways or gene sets are of primary interest

**Gene ID Format Considerations:**

- **Sample data**: Typically uses Ensembl gene IDs (e.g., "ENSG00000141510")
- **Cancer databases**: Use gene symbols (e.g., "TP53", "BRCA1")  
- **Custom sets**: Must match sample data format for proper analysis
- **Error handling**: Enhanced system provides clear guidance for format mismatches

**Selection Strategy:**

1. **Start with cancer genes** for clinical samples to assess clinically relevant batch effects
2. **Use top variable genes** for comprehensive validation and research contexts
3. **Apply custom genes** for pathway-specific or targeted analysis when specific gene sets are critical

---

#### Batch Effect Assessment and Plot Interpretation

For comprehensive guidance on batch effect assessment including detailed plot interpretation, statistical guidelines, and clinical decision frameworks, see the **[Batch Effect Assessment Guide](batch_assessment.md)**.

**Key Topics Covered:**
- Diagnostic plot types (PCA, expression distribution, Z-score, correlation)  
- Statistical interpretation guidelines and risk assessment matrix
- Clinical decision framework for `--batch_rm` parameter usage
- Example analysis with TEST data showing expected results
- Advanced considerations for protocol matching and validation