# Mummy (the wrapped MONSTER)

A pipeline to run genome-wide burdent tests using [MONSTER](http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/) (v1.3 released on December 17, 2016). MONSTER ( (Minimum P-value Optimized Nuisance parameter Score Test Extended to Relatives) is a method for detecting association between a set of rare variants and a quantitative trait in samples that contain related individuals. MONSTER is based on a mixed-effects model that includes additive and environmental components of variance and adjustment for covariates. It can handle essentially arbitrary combinations of related and unrelated individuals, including small outbred pedigrees and unrelated individuals, as well as large, complex inbred pedigrees (Duo Jiang and Mary Sara McPeek, [DOI: 10.1002/gepi.21775](http://onlinelibrary.wiley.com/doi/10.1002/gepi.21775/full)).

Mummy robust and highly customizable pipeline that allows users to define which variants should be included in the association test based on the overlapping genomic feature (eg. GENCODE annotation, if the annotation belongs to a canoncial transcript, overlap with associated regulatory feature etc.), variant feature (MAF threshold, missingness threshold) and adds custom weights (CADD, phred-scaled CADD, Eigen, phred-scaled Eigen or Linsight scores).

<b>WARNING: This pipeline is designed for GRCh38!!</b>

## Requirements

The following programs have to be in the path:

* [liftOver](https://genome.sph.umich.edu/wiki/LiftOver)
* [MONSTER](https://www.stat.uchicago.edu/~mcpeek/software/MONSTER/)
* [tabix](http://www.htslib.org/doc/tabix.html)
* [bedtools](http://bedtools.readthedocs.io/en/latest/)

The following items should be installed/downloaded:

* [bigWigTools](https://genome.ucsc.edu/goldenpath/help/bigWig.html)
* [Linsight genome-wide scores](http://compgen.cshl.edu/~yihuang/LINSIGHT/)
* [CADD genome-wide scores](http://cadd.gs.washington.edu/)
* [Eigen scores computed genome-wide](http://www.columbia.edu/~ii2135/eigen.html) (The downloaded Eigen scores have to be processed... see details below)
* linked features file (preparation of the file is detailed below.)

The filtering of variants partially based on the annotation found in the vcf files. The following INFO fields have to be present:

* **consequence** - the most severe [consequence](http://www.ensembl.org/info/genome/variation/predicted_data.html) assigned to the variant based on Ensembl [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html)
* **lof** - [loftee](https://github.com/konradjk/loftee) loss-of-function annotation (HC and LC for the high confidence and low confidence variants respectively)
* **AN**, **AC** and **AF** for filtering for allele frequency and missingness.

## Setting up the pipeline:

### Generating linked features file:

This file contains information which genomic regions can be linked to a gene. Eg: a gene is linked to its exons, CDs, transcript and also the regulatory features that overlap with the gene plus other regulatory features that overlap with variants that are known eQTLs of the gene (based on GTEx data). The following script takes all these information and combines it together using various sources: [GENCODE](http://www.gencodegenes.org/), [APPRIS](http://appris.bioinfo.cnio.es/#/), [Ensembl regulation](http://uswest.ensembl.org/info/genome/funcgen/index.html), [GTEx](gtexportal.org). Except GTEx, the data is accessed directly from the web, but the GTEx data has to be downloaded by the user and point to it when calling the script.

**Requirements for this scrip:** bgzip, tabix, bedtools, liftOver.

**Usage:** `./prepare_regions.sh -G <path to GTEx file> -o <Output folder>`

**For more information:** `./prepare_regions.sh -h`

``

### Generating the Phred-scaled Eigen scores:

...

### Adjusting the config file used by the pipeline:

The location of the score files, linked feature file, path to bigWigTools etc. 

