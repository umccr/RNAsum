# Differentially expressed genes

This section is for discussion. Its inclusion in the report depends on the robustness of differential expression analysis methods for **single-subject (ss) transcriptome** data described below.

In order to provide patient-specific information on differentially expressed (DE) genes, it is necessary to conduct differential expression analysis at the individual level. This, however, is tricky from the statistical point of view since gene-level variance, overdispersion and/or other parameters
requiring multiple subject cannot be estimated from ss transcriptome, opposed to the conventional differential expression analysis methods focusing on inter-group comparison, which are designed to identify DE genes of ‘average responses across patients’ and rely on well-powered cohorts of both cases and controls. Several methods have been developed for identification of DE genes from transcriptome data without the large cohort requirement. The *rank*-based and *Probability of Expression (PoE)*-based approaches are of particular interest:

* [RankComp](https://academic.oup.com/bioinformatics/article/31/1/62/2364816#supplementary-data)

>It requires two inputs: (1) a disease sample and (2) a set of accumulated normal samples, which can be can be accrued during the same experiment or a priori from various external resources. RankComp begins by ranking genes within the samples (both the case and the normal) according to increasing expression values. Next, pairwise rank comparison are performed to identify (a) stable gene pairs, and (b) reversal gene pairs. Stable gene pairs are defined as those with the same ordering in 99% of the accumulated normal samples [expressiongeneA > expressiongeneB] while reversal gene pairs are identified by disruption of that ordering in the disease sample [expressiongeneA < expressiongeneB]. Fisher’s exact test is conducted to test the null hypothesis that the numbers of reversal gene pairs supporting its upregulation or downregulation are equal. This procedure enables extraction of a list of DE genes for a single subject, and interpretable results can be obtained through manual examination or by
performing gene set enrichment analyses.

* PoE scale-based from [metaArray](http://www.bioconductor.org/packages/release/bioc/vignettes/metaArray/inst/doc/metaArray.pdf)

>The expression data is ranked based on the probability of expression scores computed using *three component normal-uniform mixture model* under Bayesian hierarchical analysis (involves Markov Chain Monte Carlo (MCMC) techniques, [MCMC](https://rdrr.io/bioc/metaArray/man/poe.mcmc.html) in [metaArray](https://rdrr.io/bioc/metaArray/) package) or  *two-component normal-uniform mixture distribution* algorithm (a faster algorithm based on the expectation-maximization (EM), [EM](https://rdrr.io/bioc/metaArray/man/poe.em.html) using [fit.em](https://rdrr.io/bioc/metaArray/man/fit.em.html) function) ([1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2246152/pdf/1471-2105-8-364.pdf)).

* [NOISeq-sim](https://bioconductor.org/packages/devel/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf)

>It simulates replicates under the assumption that gene expression counts follow multinomial distribution, in which the probability of each gene corresponds to the probability of a read mapping to that gene. The probability of each gene is estimated by the proportion of its read counts relative to the total number of mapped reads from the only sample under the corresponding experimental condition. The method generates a joint null distribution of fold-changes (M) and absolute differences (D) of the expression counts from the replicates within the same condition. This joint null distribution is then used to assess differential expression by gene‘s (M, D) pair computed between conditions (plots on page 23 in [NOISeq manual](https://bioconductor.org/packages/devel/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf)).


A neat summary of methods applicable to ss transcriptome data are described in [this](http://www.lussierlab.net/publications/Developing%20a%20personalome%20for%20precision%20medicine.pdf) review by [Lussier Group](http://www.lussiergroup.org/).

Work in progress...

## Table of contents

<!-- vim-markdown-toc GFM -->
* [Scripts summary](#scripts-summary)
* [Post-process, summarise and visualise sample data](#post-process-summarise-and-visualise-sample-data)
  * [Report content](#report-content)
  * [Example report](#example-report)
  * [To-do list](#to-do-list)


<!-- vim-markdown-toc -->
<br>


## Scripts summary

To be added

Script | Description | Packages
------------ | ------------ | ------------
... | ... | To be completed

<br />


## Post-process summarise and visualise sample data


**Script**: 

Argument | Description
------------ | ------------
... | To be completed
<br />

**Command line use example**:

```
```
<br>

Output file | Description
------------ | -----------
... | To be completed
<br />

### Report content

To be added

### Example report

To be added

```
```


### To-do list

* volcano plot (interactive)
* per-chromosome Manhattan plot (similar to the one on page 24 in [NOISeq manual](https://bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf)) (interactive)