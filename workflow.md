## RNA-seq report workflow

The description of the main workflow components involved in (**1**) *[read counts](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)* and *[gene fusions](./data/test_data/final/test_sample_WTS/pizzly/test_sample_WTS-flat.tsv)* data **[processing](#1-data-processing)**, (**2**) **[integration](#2-integration-with-wgs-based-results)** with **WGS**-based data (processed using *[umccrise](https://github.com/umccr/umccrise)* pipeline), (**3**) results **[annotation](#3-results-annotation)** and (**4**) presentation in the *Transcriptome Patient Summary* **[report](#4-report-generation)**. 

<img src="img/RNAseq_report_workflow.png" width="100%"> 

<br/>

## Table of contents

<!-- vim-markdown-toc GFM -->
* [1. Data processing](#1-data-processing)
    * [Counts processing](#counts-processing)
    	* [Data collection](#data-collection)
    	* [Transformation](#transformation)
    	* [Filtering (optional)](#filtering-optional)
    	* [Normalisation (optional)](#normalisation-optional)
    	* [Combination](#combination)
    	* [Batch-effects correction (optional)](#batch-effects-correction-optional)
    	* [Data scaling](#data-scaling)
    * [Fusion genes re-quantification](#fusion-genes-re-quantification)
* [2. Integration with WGS-based results](#2-integration-with-wgs-based-results)
	* [Somatic SNVs and small indels](#somatic-snvs-and-small-indels)
	* [Structural variants](#structural-variants)
	* [Somatic CNVs](#somatic-cnvs)
* [3. Results annotation](#3-results-annotation)
* [4. Report generation](#4-report-generation)

<!-- vim-markdown-toc -->

## 1. Data processing 

### Counts processing

The **read count** data (see [Input data](./README.md#input-data) section in the main page) in *[abundance.tsv](./data/test_data/final/test_sample_WTS/kallisto/abundance.tsv)* quantification file from [kallisto](https://pachterlab.github.io/kallisto/about) are processed following steps illustrated in [Figure 1](./img/counts_post-processing_scheme.png) and described below.

<img src="img/counts_post-processing_scheme.png" width="40%"> 

###### Figure 1
>Counts processing scheme.

#### Data collection

([Figure 1](./img/counts_post-processing_scheme.png)A)

* Load read count files from the following three sets of data:

	1. patient **sample** (see [Input data](./README.md#input-data) section in the main page)
	2. **external reference** cohort ([TCGA](https://tcga-data.nci.nih.gov/), available cancer types are listed in [TCGA projects summary table](./TCGA_projects_summary.md)) corresponding to the patient cancer sample.
	3. UMCCR **internal reference** set of in-house pancreatic cancer samples (regardless of the patient sample origin; see [Input data](./README.md#input-data) section in the main page)

#### Transformation

([Figure 1](./img/counts_post-processing_scheme.png)B)

* Subset datasets to include common genes
* Combine patient **sample** and **internal reference** dataset
* Convert counts to **[CPM](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/)** (*Counts Per Million*; default) or **[TPM](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/)** (*Transcripts Per Kilobase Million*) values in:
	1. **sample** + **internal reference** set
	2. **external reference** set

#### Filtering (optional)

([Figure 1](./img/counts_post-processing_scheme.png)C)

* Filter out genes with low counts (CPM or TPM **< 1** in more than 90% of samples) in:
	1. **sample** + **internal reference** set
	2. **external reference** set

#### Normalisation (optional)

([Figure 1](./img/counts_post-processing_scheme.png)D)

* Normalise data (see [Arguments](./README.md#arguments) section in the main page for available options) for sample-specific effects in:
	1. **sample** + **internal reference** set
	2. **external reference** set

#### Combination

([Figure 1](./img/counts_post-processing_scheme.png)E)

* Subset datasets to include common genes
* Combine **sample** + **internal reference** set with **external reference** set

#### Batch-effects correction (optional)

([Figure 1](./img/counts_post-processing_scheme.png)F)

* Consider the patient **sample** + **internal reference** (regardless of the patient sample origin) as one batch (both sets processed with the same pipeline) and corresponding [TCGA](https://tcga-data.nci.nih.gov/) dataset as another batch. The objective is to remove data variation due to technical factors.

#### Data scaling

The processed count data is scaled to facilitate expression values interpretation. The data is either scaled **[gene-wise](#gene-wise-z-scoreztransformation)** (Z-score transformation, default) or **[group-wise](#group-wise-centering)** (centering).

##### Gene-wise

Z-scores are comparable by measuring the observations in multiples of the standard deviation of given sample. The gene-wise Z-score transformation procedure is illustrated in [Figure 2](./img/Z-score_transformation_gene_wise.png) and is described below.

<img src="img/Z-score_transformation_gene_wise.png" width="30%"> 

###### Figure 2
>Gene-wise Z-score transformation scheme.

* Extract expression values across all samples for a given **gene** ([Figure 2](./img/Z-score_transformation_gene_wise.png)A)
* Compute **Z-scores** for individual samples (see equation in ([Figure 2](./img/Z-score_transformation_gene_wise.png)B)
* Compute **median Z-scores** for ([Figure 2](./img/Z-score_transformation_gene_wise.png)C):
	1.  **external reference** set
	2. **internal reference** set\*
* Present patient sample **Z-score** in the context the reference cohorts' **median Z-scores** ([Figure 2](./img/Z-score_transformation_gene_wise.png)D)

\* used only for pancreatic cancer patients

##### Group-wise

The group-wise centering apporach is presented in [Figure 3](./img/centering_group_wise.png) and is described below.


<img src="img/centering_group_wise.png" width="30%"> 

###### Figure 3
>Group-wise centering scheme.

* Extract expression values for each **group** ([Figure 3](./img/centering_group_wise.png)A)
	1. patient **sample**
	2. **external reference** set
	3. **internal reference** set\*
* For each gene compute **mean expression** value in individual groups ([Figure 3](./img/centering_group_wise.png)B)
* **Center** the mean expression values for each gene in individual groups ([Figure 3](./img/centering_group_wise.png)C)
* Present patient sample **centered** expression values in the context the reference cohorts' **centered** values ([Figure 3](./img/centering_group_wise.png)D)

\* used only for pancreatic cancer patients

### Fusion genes re-quantification

* **Add fusion** transcripts detected by [pizzly](https://github.com/pmelsted/pizzly) to the reference transcriptome sequences
* **Re-quantify** abundances of transcripts, including fusion transcripts detected by [pizzly](https://github.com/pmelsted/pizzly) using [kallisto](https://pachterlab.github.io/kallisto/about).
* Use the re-quantified fusion transcripts (see example [abundance.tsv](./data/test_data/final/test_sample_WTS/kallisto/quant_pizzly_post/abundance.tsv) file) for presenting expression levels of candidate fusion genes

## 2. Integration with WGS-based results

For patients with available [WGS](./README.md#wgs) data processed using *[umccrise](https://github.com/umccr/umccrise)* pipeline (see ```--umccrise``` [argument](README.md/#arguments)) the expression level information for [mutated](#somatic-snvs-and-small-indels) genes or genes located within detected [structural variants](#structural-variants) (SVs) or [copy-number](#somatic-cnvs) (CN) [altered regions](#somatic-cnvs), as well as the genome-based findings are incorporated and used s a primary source for expression profiles prioritisation.

### Somatic SNVs and small indels

* Check if **[PCGR](https://github.com/sigven/pcgr)** output file (see [example](./data/test_data/umccrised/test_sample_WGS/pcgr/test_sample_WGS-somatic.pcgr.snvs_indels.tiers.tsv)) is available
* **Extract** expression level **information** and genome-based findings for genes with detected genomic variants (use ```--pcgr_tier``` [argument](README.md/#arguments) to define [tier](https://pcgr.readthedocs.io/en/latest/tier_systems.html#tier-model-2-pcgr-acmg) threshold value)
* **Ordered genes** by increasing variants **[tier](https://pcgr.readthedocs.io/en/latest/tier_systems.html#tier-model-2-pcgr-acmg)** and then by decreasing absolute values representing difference between expression levels in the patient sample and the corresponding reference cohort

### Structural variants

* Check if **[Manta](https://github.com/Illumina/manta)** output file (see [example](./data/test_data/umccrised/test_sample_WGS/structural/test_sample_WGS-sv-prioritize-manta-pass.tsv)) is available
* **Extract** expression level **information** and genome-based findings for genes located within detected SVs
* **Ordered genes** by increasing **[SV score](https://github.com/vladsaveliev/simple_sv_annotation)** and then by decreasing absolute values representing difference between expression levels in the patient sample and the corresponding reference cohort
* **Compare** [gene fusions](./fusions) detected in WTS data ([pizzly](https://github.com/pmelsted/pizzly)) and WGS data ([Manta](https://github.com/Illumina/manta))
* **Priritise** WGS-supported [gene fusions](./fusions)

### Somatic CNVs

* Check if **[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator)** output file (see [example](./data/test_data/umccrised/test_sample_WGS/purple/test_sample_WGS.purple.gene.cnv)) is available
* **Extract** expression level **information** and genome-based findings for genes located within detected CNVs (use ```--cn_loss ``` and ```--cn_gain ``` [arguments](README.md/#arguments) to define CN threshold values to classify genes within lost and gained regions)
* **Ordered genes** by increasing (for genes within lost regions) or decreasing (for genes within gained regions) **[CN](https://github.com/umccr/umccrise/blob/master/workflow.md#somatic-cnv)** and then by decreasing absolute values representing difference between expression levels in the patient sample and the corresponding reference cohort

## 3. Results annotation

* [UMCCR key cancer genes set](https://github.com/umccr/umccrise/blob/master/workflow.md#somatic-cnv) build of off several sources:
	* Cancermine with at least 2 publication with at least 3 citations
	* NCG known cancer genes,
	* Tier 1 COSMIC Cancer Gene Census (CGC)
	* CACAO hotspot genes (curated from ClinVar, CiViC, cancerhotspots)
	* At least 2 matches in the following 5 sources and 8 clinical panels:
		* Cancer predisposition genes (CPSR list)
		* COSMIC Cancer Gene Census (tier 2)
		* AZ300
		* Familial Cancer
		* OncoKB annotated
		* MSKC-IMPACT
		* MSKC-Heme
		* PMCC-CCP
		* Illumina-TS500
		* TEMPUS
		* Foundation One
		* Foundation Heme
		* Vogelstein
* [OncoKB](https://oncokb.org/)
* [The Variant Interpretation for Cancer Consortium](https://cancervariants.org/) (VICC)
* [CIViC](https://civicdb.org/)
* [Cancer Genome Interpreter](https://www.cancergenomeinterpreter.org/biomarkers) (CGI) database
* [FusionGDB](https://ccsm.uth.edu/FusionGDB/)

### 4. Report generation

The generated html-based ***Transcriptome Patient Summary*** **report** contains several sections described more in detail in [report_structure.md](report_structure.md).
