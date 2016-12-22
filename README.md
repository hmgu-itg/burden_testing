# Rare-variant burden testing

A collection a scripts to run a gene-based genome-wide burden tests using MONSTER v1.2 (released on June 26, 2015). For more information about MONSTER see: [website](http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/) and the [publication](http://onlinelibrary.wiley.com/doi/10.1002/gepi.21775/abstract) of the McPeak group. 


This robust and highly customizable pipeline allows users to define which variants should be included in the association test based on the overlapping genomic feature (eg. GENCODE annotation, if the annotation belongs to a canoncial transcript, overlap with associated regulatory feature etc.), variant feature (MAF threshold, missingness threshold) and adds custom weights (CADD, phred-scaled CADD, Eigen, phred-scaled Eigen, GWAVA).

**Steps**

* Generating a file with annotated genomic features.
* Selecting variants.
* Adding scores if desired.
* Prepare suitable input format for MONSTER.
* Call MONSTER, process outputs, catch and report issues.

## Annotate genomic regions using **prepare_regions.sh**


As an initial step an annotation file is generated based on which we can select regions of interest for we consider overlapping variants.  

This script downoads and annotate genomic regions to allow precise feature selection. 
Sources:
GENCODE
APPRIS
Ensembl Regulation


The annotated features are 



