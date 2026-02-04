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
| **Reference** | [TCGA projects summary table](./inst/articles/TCGA_projects_summary.md) |

**Processing pipeline:**
All TCGA datasets have been harmonized using methods documented in the [TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/blob/master/expression/README.md#gdc-counts-data) repository, ensuring consistency and reproducibility.

---

### Internal reference cohort

To minimize batch effects and improve comparability, one can provide an internal reference cohort processed with identical pipelines to the input data.

#### Why an internal cohort?

Public TCGA datasets may show significant [batch effects](https://www.ncbi.nlm.nih.gov/pubmed/20838408) when compared to in-house data due to:

- **Different experimental procedures** and protocols
- **Variable tissue quality** and cellularity
- **Different analytical pipelines**

The internal cohort addresses these limitations by ensuring:

- Same data generation protocols
- Consistent sample quality standards
- Identical bioinformatics processing

#### Key Advantages

✓ **Minimized batch effects** - same protocols as your input data  
✓ **High tissue quality** - stringent QC standards  
✓ **Matched processing** - [data pre-processing](./inst/articles/workflow.md#data-processing) identical to patient samples  

---

### Choosing the right reference

- Use the matching TCGA **cancer-specific** dataset
- Consider **Pan-Cancer** for cross-cancer comparisons or rare tumour types
- Use the **internal cohort** for maximum comparability and minimal batch effects

