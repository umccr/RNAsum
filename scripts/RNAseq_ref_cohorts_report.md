# Presenting sample data in the context of reference datasets

Set of scripts to presnet the sample expression profile in the context of genome-scale pancreas-related data from large cohorts of patients. Currenlty the sample data is combined with [TCGA](https://cancergenome.nih.gov/) data but once we optimise the *batch-effect adjustment apporach* in the [RNA-seq data combination step](https://github.com/umccr/RNA-seq-analysis/tree/master/readcount-analysis) the plan is to add more datasets, e.g. [14 metastatic PC](https://met500.path.med.umich.edu/datasets) samples from [MET500 project](https://met500.path.med.umich.edu/datasets) (published by [Dan R. Robinson et al. (2017)](https://www-nature-com.ezp.lib.unimelb.edu.au/articles/nature23306)) and 244 (?) normal pancreas samples from [recount2](https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/RNAseq/doc/recount-workshop.html).


## Table of contents

<!-- vim-markdown-toc GFM -->
* [Scripts summary](#scripts-summary)
* [Post-process, summarise and visualise sample data](#post-process,-summarise-and-visualise-sample-data)
  * [Report content](#report-content)
  * [Example report](#example-report)
  * [To-do list](#to-do-list)


<!-- vim-markdown-toc -->
<br>


## Scripts summary

Script | Description | Packages
------------ | ------------ | ------------
*[RNAseq_ref_cohorts_report.R](RNAseq_ref_cohorts_report.R)* | Collects user-defined parameters for the corresponding [RNAseq_ref_cohorts_report.Rmd](RNAseq_ref_cohorts_report.Rmd) markdown script |  *[optparse](https://cran.r-project.org/web/packages/optparse/optparse.pdf)* <br> *[knitr](https://cran.r-project.org/web/packages/knitr/knitr.pdf)*
*[RNAseq_ref_cohorts_report.Rmd](RNAseq_ref_cohorts_report.Rmd)* | Launch by *[RNAseq_ref_cohorts_report.R](RNAseq_ref_cohorts_report.R)*. Post-processes, summarises and visualises an output from *[bcbio-nextgen RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* to generate patient summary report <br> | *[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)* <br> *[preprocessCore](https://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html)* <br> *[plotly](https://plot.ly/r/)* <br> *[ClassDiscovery](https://cran.r-project.org/web/packages/ClassDiscovery/index.html)* <br> *[plotly](https://plot.ly/r/)* <br> *[made4](https://bioconductor.org/packages/release/bioc/html/made4.html)* <br> *[ade4](https://cran.r-project.org/web/packages/ade4/index.html)*
<br />


## Post-process, summarise and visualise sample data

To summarise multiple MAF files run the *[summariseMAFs.R](https://github.com/umccr/MAF-summary/tree/master/scripts/summariseMAFs.R)* script. This script catches the arguments from the command line and passes them to the *[summariseMAFs.Rmd](https://github.com/umccr/MAF-summary/tree/master/scripts/summariseMAFs.Rmd)* script to produce the html report, generate set of plots and excel spreadsheets summarising each MAF file.

**Script**: *[RNAseq_ref_cohorts_report.R](RNAseq_ref_cohorts_report.R)*

Argument | Description
------------ | ------------
--sample_name | Desired sample name to be presented in the report
--count_file | Location and name of the read count file from *[bcbio-nextgen RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)*
--report_dir | Desired location for the report
<br />

**Command line use example**:

```
Rscript  RNAseq_ref_cohorts_report.R  --sample_name ./CCR170012_MH17T001P013  --count_file ./data/CCR170012_MH17T001P013-ready.counts  --report_dir /reports
```
<br>

This will generate *[CCR170012_MH17T001P013.RNAseq_ref_cohorts_report.html](../reports/CCR170012_MH17T001P013.RNAseq_ref_cohorts_report.html)* in */reports* folder. It will also create additional output tables and plots:

Output file | Description
------------ | -----------
... | To be completed
<br />

### Report content

Currnetly this part of the report contains two sections:

1. **Comparison across tumour types**

Principal component analysis (PCA), hierarchical clustering and between-group analysis (BGA) to present the sample data in the context of patients with various tumour types derived from [TCGA](https://cancergenome.nih.gov/) cohorts. The sample expression profile is projected against the four patient cohorts representing (1) unrelated tumour types (*distant tumours*), (2) *core gastrointestinal* tumours, (3) *developmental gastrointestinal* tumours and (4) *pancreas-related* lesions/tissues, and are presented in separate tabs. Each tab has the following sub-sections:

Tab | Description
------------ | -----------
PCA | reduces the dimensionality of data while retaining most of the variation in the dataset, making it possible to visually assess similarities and differences between the investigated sample and the various patient cohorts
Dendrogram | hierarchical clustering with reference samples coloured according to corresponding tumour types
BGA | supervised classification method ([Culhane et al., (2002)](https://www.ncbi.nlm.nih.gov/pubmed/12490444))
<br />



2. **Molecular classification**

PCA, hierarchical clustering and BGA to project sample in the context of [TCGA](https://cancergenome.nih.gov/) PAAD samples classified based mRNA subtypes reported by [Bailey et al. (2016)](https://www.ncbi.nlm.nih.gov/pubmed/26909576), [Moffitt et al. (2015)](https://www.ncbi.nlm.nih.gov/pubmed/26343385) and [Collisson et al. (2011)](https://www.ncbi.nlm.nih.gov/pubmed/21460848). The molecular classification aids patient assignment into less heterogeneous and more appropriate group regarding the metastatic risk and the therapeutic response, with the consequences of better predicting evolution and better orienting the treatment. A recent report by [Birnbaum DJ1 et al. (2018)}(https://www.ncbi.nlm.nih.gov/pubmed/29499330) reviews the association between pancreatic ductal adenocarcinoma (PDAC) molecular subtypes and drugs sensitivity. Individual classification results are presented in separate tabs with the corresponding tumour subtypes description and containing sub-sections as in the [comparison across tumour types](#report-content) section.
	

### Example report

Example read count data from *[bcbio-nextgen RNA-seq pipeline](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#rna-seq)* is located on [Spartan](https://dashboard.hpc.unimelb.edu.au/) :

```
/data/cephfs/punim0010/projects/Jacek_RNA_seq_report/RNAseq-Analysis-Report/data
```

* [HTML report](https://rawgit.com/umccr/MAF-summary/master/scripts/summariseMAFs.html) - R html report for all cohorts

### To-do list

* [Cumulative distribution function](https://en.wikipedia.org/wiki/Cumulative_distribution_function) (CDF) plots for selected genes to present their expression levels in the context of the overall expression distribution in investigated sample
* Link to [umccrise](https://github.com/umccr/umccrise) report to define genes to for the  CDF plots
* ...
