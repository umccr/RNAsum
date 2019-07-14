# TCGA projects summary


The table below summarises [TCGA](https://portal.gdc.cancer.gov/) expression data available for **33 cancer types**. The dataset of interest can be specified by using one of the [TCGA](https://portal.gdc.cancer.gov/) project IDs (`Project` column) for the `--dataset` argument in *[RNAseq_report.R](./rmd_files/RNAseq_report.R)* script (see [Arguments](./README.md#arguments) section).


No | Project | Name | Tissue code\* | Samples no.\**
------------ | ------------ | ------------ | ------------ | ------------
1 | `BRCA`  | Breast Invasive Carcinoma | 1 | **1171**
2 | `THCA`  | Thyroid Carcinoma | 1 | **546**
3 | `HNSC`  | Head and Neck Squamous Cell Carcinoma | 1 | **524**
4 | `LGG`   | Brain Lower Grade Glioma | 1 | **518**
5 | `KIRC`  | Kidney Renal Clear Cell Carcinoma | 1 | **505**
6 | `LUSC`  | Lung Squamous Cell Carcinoma | 1 | **487**
7 | `LUAD`  | Lung Adenocarcinoma | 1 | **444**
8 | `PRAD`  | Prostate Adenocarcinoma | 1 | **433**
9 | `SKCM`  | Skin Cutaneous Melanoma | 1 | **429**
10 | `STAD`  | Stomach Adenocarcinoma | 1 | **403**
11 | `LIHC`  | Liver Hepatocellular Carcinoma | 1 | **400**
12 | `BLCA`  | Bladder Urothelial Carcinoma | 1 | **372**
13 | `CESC`  | Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma | 1 | **301**
14 | `COAD`  | Colon Adenocarcinoma | 1 | **289**
15 | `KIRP` | Kidney Renal Papillary Cell Carcinoma | 1 | **283**
16 | `OV`    | Ovarian Serous Cystadenocarcinoma | 1 | **254**
17 | `SARC` | Sarcoma | 1 | **218**
18 | `PCPG`  | Pheochromocytoma and Paraganglioma | 1 | **186**
19 | `UCEC`  | Uterine Corpus Endometrial Carcinoma | 1 | **184**
20 | `LAML`  | Acute Myeloid Leukaemia | 3 | **173**
21 | `ESCA`  | Esophageal Carcinoma | 1 | **170**
22 | `PAAD`  | Pancreatic Adenocarcinoma | 1 | **156**
23 | `GBM`   | Glioblastoma Multiforme | 1 | **155**
24 | `TGCT`  | Testicular Germ Cell Tumours | 1 | **155**
25 | `THYM`  | Thymoma | 1 | **120**
26 | `READ` | Rectum Adenocarcinoma | 1 | **96**
27 | `KICH`  | Kidney Chromophobe | 1 | **83**
28 | `UVM`   | Uveal Melanoma | 1 | **80**
29 | `MESO` | Mesothelioma | 1 | **78**
30 | `ACC`   | Adrenocortical Carcinoma | 1 | **78**
31 | `UCS`   | Uterine Carcinosarcoma | 1 | **57**
32 | `DLBC`  | Lymphoid Neoplasm Diffuse Large B-cell Lymphoma | 1 | **47**
33 | `CHOL`  | Cholangiocarcinoma | 1 | **43**
<br />

\* Tissue codes:

Tissue code | Letter code | Definition
------------ | ------------ | ------------
1 | TP  | Primary solid Tumour
3 | TB  | Primary Blood Derived Cancer - Peripheral Blood
<br />

\** Each dataset was cleaned based on the quality metrics provided in the *Merged Sample Quality Annotations* file **[merged_sample_quality_annotations.tsv](http://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf)** from [TCGA Pan-Cancer Clinical Data Resource](https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018) (see [TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/tree/master/expression/README.md#data-clean-up) repository for more details).
 
 