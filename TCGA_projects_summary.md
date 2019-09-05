# TCGA projects summary


The table below summarises [TCGA](https://portal.gdc.cancer.gov/) expression data available for **[33 cancer types](#primary-datasets)**. Additionally, for *Bladder Urothelial Carcinoma* and *Pancreatic Adenocarcinoma* cohorts extended sets, including neuroendocrine tumours (NETs), intraductal papillary mucinous neoplasm (IPMNs) and acinar cell carcinoma (ACC) samples, are available (see [Extended datasets](#extended-datasets) table).

The dataset of interest can be specified by using one of the [TCGA](https://portal.gdc.cancer.gov/) project IDs (`Project` column) for the `--dataset` argument in *[RNAseq_report.R](./rmd_files/RNAseq_report.R)* script (see [Arguments](./README.md#arguments) section).

###### Note

To readuce the data processing time and the size of the final html-based ***Patient Transcriptome Summary*** **report** the following datasets were restricted to inlcude expression data from 300 patients: `BRCA`, `THCA`, `HNSC`, 
`LGG`, `KIRC`, `LUSC`, `LUAD`, `PRAD`, `SKCM`, `STAD` and `LIHC`.

## Primary datasets

No | Project | Name | Tissue code\* | Samples no.\**
------------ | ------------ | ------------ | ------------ | ------------
1 | `BRCA`  | Breast Invasive Carcinoma | 1 | **300**
2 | `THCA`  | Thyroid Carcinoma | 1 | **300** (verify)
3 | `HNSC`  | Head and Neck Squamous Cell Carcinoma | 1 | **300** (verify)
4 | `LGG`   | Brain Lower Grade Glioma | 1 | **300** (verify)
5 | `KIRC`  | Kidney Renal Clear Cell Carcinoma | 1 | **300** (verify)
6 | `LUSC`  | Lung Squamous Cell Carcinoma | 1 | **300** (verify)
7 | `LUAD`  | Lung Adenocarcinoma | 1 | **300**
8 | `PRAD`  | Prostate Adenocarcinoma | 1 | **300** (verify)
9 | `SKCM`  | Skin Cutaneous Melanoma | 1 | **300** (verify)
10 | `STAD`  | Stomach Adenocarcinoma | 1 | **300** (verify)
11 | `LIHC`  | Liver Hepatocellular Carcinoma | 1 | **300** (verify)
12 | `KIRP` | Kidney Renal Papillary Cell Carcinoma | 1 | **283** (verify)
13 | `COAD`  | Colon Adenocarcinoma | 1 | **257**
14 | `BLCA`  | Bladder Urothelial Carcinoma | 1 | **246**
15 | `OV`    | Ovarian Serous Cystadenocarcinoma | 1 | **220**
16 | `SARC` | Sarcoma | 1 | **214**
17 | `PCPG`  | Pheochromocytoma and Paraganglioma | 1 | **186** (verify)
18 | `UCEC`  | Uterine Corpus Endometrial Carcinoma | 1 | **184** (verify)
19 | `LAML`  | Acute Myeloid Leukaemia | 3 | **173** (verify)
20 | `CESC`  | Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma | 1 | **171**
21 | `ESCA`  | Esophageal Carcinoma | 1 | **170** (verify)
22 | `GBM`   | Glioblastoma Multiforme | 1 | **155** (verify)
23 | `TGCT`  | Testicular Germ Cell Tumours | 1 | **155** (verify)
24 | `PAAD`  | Pancreatic Adenocarcinoma | 1 | **150**
25 | `THYM`  | Thymoma | 1 | **120** (verify)
26 | `READ` | Rectum Adenocarcinoma | 1 | **87**
27 | `KICH`  | Kidney Chromophobe | 1 | **83** (verify)
28 | `UVM`   | Uveal Melanoma | 1 | **80** (verify)
29 | `MESO` | Mesothelioma | 1 | **78** (verify)
30 | `ACC`   | Adrenocortical Carcinoma | 1 | **78** (verify)
31 | `UCS`   | Uterine Carcinosarcoma | 1 | **57** (verify)
32 | `DLBC`  | Lymphoid Neoplasm Diffuse Large B-cell Lymphoma | 1 | **47** (verify)
33 | `CHOL`  | Cholangiocarcinoma | 1 | **43** (verify)
<br />

## Extended datasets

No | Project | Name | Tissue code\* | Samples no.\**
------------ | ------------ | ------------ | ------------ | ------------
1 | `BLCA-NET`  | Bladder Urothelial Carcinoma dataset including neuroendocrine tumours (NETs, n=2) | 1 | **248**
2 | `PAAD-IPMN`  | Pancreatic Adenocarcinoma dataset including intraductal papillary mucinous neoplasm (IPMNs, n=2) | 1 | **152**
3 | `PAAD-NET`  | Pancreatic Adenocarcinoma dataset including neuroendocrine tumours (NETs, n=8) | 1 | **158**
4 | `PAAD-ACC`  | Pancreatic Adenocarcinoma dataset including acinar cell carcinoma (ACCs, n=1) | 1 | **151**
<br />

\* Tissue codes:

Tissue code | Letter code | Definition
------------ | ------------ | ------------
1 | TP  | Primary solid Tumour
3 | TB  | Primary Blood Derived Cancer - Peripheral Blood
<br />

\** Each dataset was cleaned based on the quality metrics provided in the *Merged Sample Quality Annotations* file **[merged_sample_quality_annotations.tsv](http://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf)** from [TCGA PanCanAtlas initiative webpage](https://gdc.cancer.gov/about-data/publications/pancanatlas) (see [TCGA-data-harmonization](https://github.com/umccr/TCGA-data-harmonization/tree/master/expression/README.md#data-clean-up) repository for more details).
 
 