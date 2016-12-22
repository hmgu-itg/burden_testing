# Rare-variant burden testing

A collection a scripts to run a gene-based genome-wide burden tests using MONSTER v1.2 (released on June 26, 2015). For more information about MONSTER see: [website](http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/) and the [publication](http://onlinelibrary.wiley.com/doi/10.1002/gepi.21775/abstract) of the McPeak group.


This robust and highly customizable pipeline allows users to define which variants should be included in the association test based on the overlapping genomic feature (eg. GENCODE annotation, if the annotation belongs to a canoncial transcript, overlap with associated regulatory feature etc.), variant feature (MAF threshold, missingness threshold) and adds custom weights (CADD, phred-scaled CADD, Eigen, phred-scaled Eigen, GWAVA).

**WARNING: This pipeline is designed for GRCh38!!**

__Steps:__

* Generating a file with annotated genomic features.
* Selecting variants.
* Adding scores if desired.
* Prepare suitable input format for MONSTER.
* Call MONSTER, process outputs, catch and report issues.

## Annotate genomic regions using **prepare_regions.sh**

As an initial step, an annotation file is generated based on which we can select regions of interest in we consider overlapping variants. The annotation is gene based: the output (bed)[http://www.ensembl.org/info/website/upload/bed.html] file has the coordinates of the genes, and in the 5th column, there is the annotation of all possible associated features for that gene. The association is based on GENCODE annotation (transcript, exon etc. of the given gene), or associated Ensembl regulatory features. A regulatory feautre is considered to be associated with a gene if

* it overlaps with the gene.
* a variant overlaps with the regualtory feature that has been found to be an eQTL of that gene in the GTEx databse.

**usage:** *./prepare_regions.sh \<target directory\>*

**Requirements:**

* loftOver
* Chainfile to lift over from hg18 to hg38 stored in the script dir.
* bgzip and tabix in path
