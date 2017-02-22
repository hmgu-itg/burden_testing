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

### Running the burden test genome-wide

Once the preparation steps are finished, paths to input files have to be adjusted in the following files:
`config.txt`, `MONSTERgenome-wide.sh`. 

Once everything is set up, run the wrapper script:

```bash
./MONSTERgenome-wide.sh -g exon -d 1 -p TG -m 0.05 -s Eigen
```

In the above example the MONSTER will be run for exons of all protein coding genes, applying a 5% upper MAF cutoff on the selected variants, then adding raw Eigen scores as weights, and using TG (triglycerides) as phenotype. The script is highly customizable in terms of the selected genomic features, variant filters and weights. For more information see:

```bash
./MONSTERgenome-wide.sh -h
```
If a gene has p-value less than 1e-5, the script considers it as a potential hit, and performs extra checks to test if the association is driven by a single variant: the script repeats the test with one snp removed at a time. It will show that the signal is attenuated if a single signal is removed. This method however does not account for variants in LD.  

### Summarizing the results

As the test is run for potentially a large number of genes, a script is included to process and summarize the results.

```bash
./Summarize.sh <output of the burden runs> <output directory for the summary>
```

This script pools all p-values together for a given trait, and for the hits it reports the highest, less significant p-values it got for the above mentioned testing.
