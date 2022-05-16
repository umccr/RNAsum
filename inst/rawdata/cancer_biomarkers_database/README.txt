#
# Catalog of Cancer Genes 
#
# https://www.cancergenomeinterpreter.org/genes
#
#Citation: link to paper
#

The Catalog of Cancer Genes contains a list of genes that drive tumorigenesis in a given tumor type(s) upon a 
certain type of somatic alteration (mutation, copy number alteration and/or translocation). The Catalog includes 
genes whose role in tumorigenesis has been clinically or experimentally validated, as well as other that have 
been detected by computational methods to carry signals of positive selection across cohorts of tumors of various 
cancer types. Information of genes with validated oncogenic effect is aggregated from the Cancer Gene Census
 (license http://cancer.sanger.ac.uk/cosmic/acceptable_use) and a manual curation effort. Computational predictions
 are based on the identification of signals of positive selection in the pattern of gene alterations via an 
ensembl of bioinformatics tools applied to dozens of sequenced tumors cohorts. Most of these genes are included 
in the IntOGen resource (license https://www.intogen.org/about).


FILES:

cancer_acronyms.tsv
----------------------------
This file contains the equivalencies between cancer type acronyms used by the Cancer Genome Interpreter and 
their common names (description column). Both, acronyms and names correspond to the hierarchical organization 
of cancer types presented in the Cancer Genome Interpreter analysis page (link to analysis). 


cancer_genes_upon_mutations_or_CNAs.tsv
-------------------------------------------------------------
List of genes known (or likely, according to bioinformatics methods) to be involved in tumorigenesis across 
cancer types via different types of alterations. Each row in the file corresponds to a gene-tumor type-alteration 
type combination (defined by the first three columns of the file). The last column details the source from 
which the information in the row has been extracted. 

Sources include:

Biomarker: Gene with mutations known to be biomarkers of response to drugs (obtained from the CGI cancer biomarkers database)
CGC: Gene obtained from the Cancer Gene Census (see above)
In_silico_predicted: Gene identified by an ensemble of bioinformatics tools to have signals of positive selection in its mutational pattern across tumors (obtained from the IntOGen repository)
Predisposing: Gene with clinical evidence of predisposing to cancer (obtained from the literature)
Validated: Experimentally validated cancer driver gene (obtained from the literature)


Cancer_genes_upon_trans.tsv
-----------------------------------------
List of genes known to be involved in tumorigenesis across cancer types via translocations. Each row in the file 
corresponds to a translocation-tumor type combination (first and third columns). The second column specifies 
the effector gene (the most likely to promote the tumorigenesis between the two translocation partners). The last 
column details the source. See above description of cancer_genes_upon_mutations_or_CNAs.tsv for details on the sources.

gene_MoA.tsv
-------------------
This file contain the information on the mechanisms of action (MoA), or role of each cancer gene in tumorigenesis
 (i.e., activating, or oncogene, and loss-of-function, or tumor suppressor). The MoA is either predicted using 
computational methods or obtained from other databases or the literature (see methods in the paper 
(link to paper) for more details).

columns:
gene - gene symbol
gene_MoA - Act (Activating), LoF (Loss of Function), ambiguous

