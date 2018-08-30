#Introduction

This document records notes and links to papers and repositories, we *think* are interesting or worth exploring/incorporating into UMCCR RNAseq analysis report.

Not to forget, another Google-doc created before is [RNAseq-Workflows](https://docs.google.com/document/d/1WqQQiDgVMIrm856xdNAdqGQibHWUIFx5Ir0xzB8pOEE/edit#heading=h.mh989jx06fxl).

##https://github.com/nf-core/rnaseq

nfcore/rnaseq is a bioinformatics analysis pipeline used for RNA sequencing data, built using Nextflow.

##Explore tximport

A couple of related links:

https://cran.r-project.org/web/packages/GREP2/vignettes/vignette.html

https://www.nature.com/articles/sdata2018136

##Copying important messages from slack

- I (Oliver) didn't see anything that I wanted to incorporate (from nf-core) -- it's a good, best practice workflow. I _did_ want to look at https://bioconductor.org/packages/release/bioc/html/dupRadar.html again as it's written by an old PhD friend of mine. I also still have to check if they ship variant sites for use with HiSat2 which we could adopt.

- Right now we are mostly interested in generating expression counts and fusion calls for WTS data. These are _single samples_, so no differential expression, limited normalisation, etc.

- The biggest problem is going to be our reference sets. We'll overestimate expression of non-unique genes by a lot if we compare Salmon/Kallisto results with data from, say, Cufflinks/featureCounts.

- A couple of pointers from Oliver (what he would like to have)
	- Generate counts with Salmon, Kallisto, Star and HiSat2 and compare them to ICGC and in-house background sets.

	- Get the PCAWG RNA-Seq data set and process it here with Salmon. Yes, we'll still have tons of batch effects from all the different prep methods, libraries and sequencers but at least we process it uniformly.