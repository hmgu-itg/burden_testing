# Rare-variant burden testing

A collection a scripts to run a gene-based genome-wide burden tests using MONSTER v1.2 (released on June 26, 2015). For more information about MONSTER see: [website](http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/) and the [publication](http://onlinelibrary.wiley.com/doi/10.1002/gepi.21775/abstract) of the McPeak group.


This robust and highly customizable pipeline allows users to define which variants should be included in the association test based on the overlapping genomic feature (eg. GENCODE annotation, if the annotation belongs to a canoncial transcript, overlap with associated regulatory feature etc.), variant feature (MAF threshold, missingness threshold) and adds custom weights (CADD, phred-scaled CADD, Eigen, phred-scaled Eigen, GWAVA).

<b>WARNING: This pipeline is designed for GRCh38!!</b>

### Prepartion

Before starting the runs, an annotation file is generated that allows the precise selection of genomic features.

__Usage:__
```bash
./prepare_regions.sh <target directory>
```
This file is created based on various sources

As an initial step, an annotation file is generated based on which we can select regions of interest in we consider overlapping variants. The annotation is gene based: the output [bed](http://www.ensembl.org/info/website/upload/bed.html) file has the coordinates of the genes, the 5th column contains the JSON formatted annotation of all possible associated features for that gene. These associated regions can be GENCODE features (transcripts, exons etc. of the given gene), or associated Ensembl regulatory features. A regulatory feautre is considered to be associated with a gene if

* it overlaps with the gene.
* a variant overlaps with the regualtory feature that has been found to be an eQTL of that gene in the GTEx databse.

(The annotation has further information that allows fine selection of the features.)

<b>usage:</b> *./prepare_regions.sh \<target directory\>*

<b>Requirements:</b>

* [liftOver](http://genome.sph.umich.edu/wiki/LiftOver), [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html), [bgzip and tabix](http://www.htslib.org/doc/tabix.html) in path
* [Chainfile](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/) to lift over from hg18 to hg38 stored in the script dir.

### Parsing the output file

The 5th column of the *bed* file contains the annotation in a regularly formatted [JSON](https://en.wikipedia.org/wiki/JSON) string. In this section the stored information is explained. Each feature has a *source* key describing from which source a given feature is linked to the gene. The containing information depend on the source tag.


<b>Features Sources:</b>

* <b>GENCODE</b>: the region is a GENCODE feature of the gene, based on the "class" it could be exon/gene/transcript etc. The "appris" tag tells if the feature belongs to a principal transcript or not. The "chr", "start" and "end" tags give the GRCh38 coordinates of the feature.
* <b>GTEx</b>: This is an Ensembl regulatory feature (Ensembl stabil ID is stored in the "regulatory_ID" tag). The "class" tag tells if it is an enhancer, promoter or something else. "Tissues" tag lists all the cell types and tissues in which the regulatory feature was found to be active. This regulatory region overlaps with a variant that was found to be an eQTL for this gene in the [GTEx](gtexportal.org/home/) database. The eQTLs are listed in the "GTEx_rsIDs". The list of tissues in which the eQTL was found is kept in "GTEx_tissues".
* <b>Overlap</b>:  This is also an Ensembl regulatory feature that overlaps with the gene. The "class" tag tells if it is an enhancer, promoter or something else. "Tissues" tag lists all the cell types and tissues in which the regulatory feature was found to be active.
