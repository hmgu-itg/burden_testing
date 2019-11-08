#!/bin/bash

## Description:
# Creating a gene based genome annotation file using information from various sources,
# that allows linking regulatory features to genes. The link between the regulatory
# feature and the gene is based on simple overlap between regulatory features and
# the gene, or a regulatory feature is considered to be associated with a gene
# if the regulatory feature overlaps with a variant that has been found to be
# and eQTL for the gene in the GTEx dataset.

# The resulting file is also contains information if the given transcript, or exon, CDS
# belong to a principal or minor transcript based on APPRIS annotation.

# This script relies on a series of online sources:
## 1. Gencode
## 2. Ensembl Regulation
## 3. Appris

# The GTEx eQTL data should be specified uing the -G switch.

# For reproductibility, both the GENCODE and the Ensembl versions are hardcoded.

# Requirements: tabix, bgzip, bedtools in the path.

##
## Warning: this version intentionally use the older (V84) Ensembl release:
#### The V85 ensembl release contained a faulty regulatory build, which although was rich in
#### cell types, the number of features were much fewer, there were only promoters and
#### enhancers were annotated as active. After a discussion with the Ensembl team we decided
#### to regress to V84. To use a newer version, the script has to be adjusted
#### as data format of subsequent releases are quite different.

##
## Date: 2019.09.16 by Andrei Barysenka. andrei.barysenka@helmholtz-muenchen.de
##
script_version=3.0
last_modified=2019.09.16

## Built in versions:
GENCODE_release=32
Ensembl_release=97 
GTExRelease=8

# Get script dir:
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## printing out information if no parameter is provided:
function usage {
    echo ""
    echo "Usage: $0 -G <GTEx file> -t <targetdir>"
    echo ""
    echo " This script was written to prepare input file for the burden testing pipeline."
    echo ""
    echo ""
    echo "Version: ${script_version}, Last modified: ${last_modified}"
    echo ""
    echo "Requirements:"
    echo "  bgzip, tabix in PATH"
#    echo "  liftOver in path"
#    echo "  hg19ToHg38.over.chain chain file in script dir"
    echo "  bedtools in PATH"
    echo "  downloaded GTEx datafile with the single eQTLs (eg. GTEx_Analysis_V6_eQTLs.tar.gz)"
    echo ""
    echo ""
    echo "Workflow:"
    echo "  1: Downloads v32 GENCODE release."
    echo "  2: Downloads V97 Ensembl Regulation release."
    echo "  3: Downloads newest APPRIS release"
    echo "  4: Adds Appris annotation to Gencode transcripts."
    echo "  5: Creates cell-specific regulatory features."
    echo "  6: Lifts over GTEx coordinates to GRCh38."
    echo "  7: Links regulatory features to genes based on GTEx data."
    echo "  8: Links regulatory features to genes based on overlapping."
    echo "  9: Combined GENCODE, GTEx and Overlap data together into a single bedfile."
    echo "  10: Tabix output, cleaning up."
    echo ""
    echo ""
    echo "The output is a bed file, where the first 4 columns are the chromosome, start/end
coordinates and the stable ID of the gene respectively. The 5th column is a json
formatted string describing one genomic region associated to the given gene. This
line contains all information of the association."
    echo ""
    echo "JSON tags:"
    echo "  -source: from which source the given region is coming from (GENCODE, GTEx, Overlap)."
    echo "  -class: class of the given feature (eg. exon, CDS, gene, enhancer, promoter etc.)"
    echo "  -chr, start, end: GRCh38 coordintes of the feature."
    echo "  -other sources contain information about the evidence. (linked rsID, tissue in
    which the feature in active etc.)"
    echo ""
    echo "WARNINGS: ALL COORDINATES ARE BASED ON GRCH38 BUILD!"
    echo ""
    exit 0
}

# Function to test if a given file exists or not in which case it reports and terminates the
# execution.
function testFile {
    if [[ ! -e "$1"  ]]; then
        echo "[Error] At this step something failed. The file was not created! $1"
        echo "[Error] Exiting."
        exit 1
    fi
}

# Checking if tabix, liftOver and bedtools are installed...
function checkCommand {
    isCommand=$( which $1 | wc -l )

    # exit program is not in path:
    if [[ $isCommand == 0 ]]; then
        echo "[Error] $1 is not in path. Install program before proceeding. Exiting.";
        exit 1;
    fi
}

# We also run a test to check if the number of lines of a temporary file is zero or not.
# If it is zero, the script exits, because it indicates there were some problems.
function testFileLines {

    # Check if file is zipped:
    IsCompressed=$( file $1 | grep compressed | wc -l)

    # Check the number of lines:
    if [[ $IsCompressed -ne 0 ]]; then
        lines=$( zcat $1 | wc -l )
    else
        lines=$( cat $1 | wc -l )
    fi

    # exit if lines are zero:
    if [[ $lines == 0 ]]; then
        echo "[Error] file ($1) contains no lines. There were errors. Exiting.";
        exit 1;
    fi
}

# This function prints out all the reports that were generated during the run (with time stamp!):
function info {
    hourMin=$(date +"%T" | awk 'BEGIN{FS=OFS=":"}{print $1, $2}')
    echo -ne "[Info ${hourMin}] $1"
}

# Printing help message if no parameters are given:
if [[ $# == 0 ]]; then usage; fi

# Processing command line options:
OPTIND=1
while getopts "G:t:h" optname; do
    case "$optname" in
        "G" ) GTExFile="${OPTARG}" ;;
        "t" ) targetDir="${OPTARG}" ;;
        "h" ) usage ;;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

# Checking the provided working directory:
if [[  ! -d "${targetDir}" ]]; then
    echo "[Error] The provided directory does not exists: $targetDir"
    exit 1
# Checking if the defined working directory is writable:
elif [[ ! -w "${targetDir}" ]]; then
    echo "[Error] The provided working directory is not writable: ${targetDir}"
    exit 1
fi

# Checking if GTEx file exists:
if [[ -z "${GTExFile}" ]]; then
    echo "[Error] The compressed GTEx file is needed! eg. GTEx_Analysis_V8_eQTLs.tar.gz"
    exit 1
elif [[ ! -e "${GTExFile}" ]]; then
    echo "[Error] The provided GTEx file (${GTExFile}) does not exist."
    exit 1
fi

# Checking required commands:
checkCommand tabix
checkCommand liftOver
checkCommand bedtools

# Checking chainfile in the scriptDir:
testFile ${scriptDir}/hg19ToHg38.over.chain

# Last step in setup:
today=$(date "+%Y.%m.%d")
info "Current date: ${today}\n"
info "Working directory: ${targetDir}/${today}\n\n"

#=================================== GENCODE =================================================

# Get the most recent version of the data:
mkdir -p ${targetDir}/${today}/GENCODE
info "Downloading GENCODE annotation from ftp://ftp.ebi.ac.uk/. Release version: ${GENCODE_release}... "
wget -q ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_release}/gencode.v${GENCODE_release}.annotation.gtf.gz -O ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz
echo -e "done."

# Testing if the file is exists or not:
testFile "${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz"

# Counting genes in the dataset:
genes=$(zcat ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | awk 'BEGIN{FS="\t";}$3 == "gene"{print $3;}' | wc -l )
info "Total number of genes in the GENCODE file: ${genes}\n\n"

#=================================== REGULATION ==============================================

# prepare target directory:
mkdir -p ${targetDir}/${today}/EnsemblRegulation
info "Downloading cell specific regulatory features from Ensembl.\n"

# Get list of all cell types:
cells=$(curl -s ftp://ftp.ensembl.org/pub/release-${Ensembl_release}/regulation/homo_sapiens/RegulatoryFeatureActivity/ | perl -lane 'print $F[-1]')

# If there are no cell types present in the downloaded set, it means there were some problems. We are exiting.
if [ -z "${cells}" ]; then
    echo "[Error] No cell types were found in the Ensembl regulation folder."
    echo "[Error] URL: ftp://ftp.ensembl.org/pub/release-${Ensembl_release}/regulation/homo_sapiens/regulatory_features/RegulatoryFeatureActivity/"
    echo "Exiting."
    exit 1
fi

#Download all cell types:
for cell in ${cells}; do
    echo -n "."
    wget -q ftp://ftp.ensembl.org/pub/release-${Ensembl_release}/regulation/homo_sapiens/RegulatoryFeatureActivity/${cell}/homo_sapiens.*Regulatory_Build.regulatory_activity.*.gff.gz -O ${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz
    testFile "${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz"
done
echo "Done."

# Printing out report of the downloaded cell types:
cellTypeCount=$(ls -la ${targetDir}/${today}/EnsemblRegulation/*gff.gz | wc -l)
info "Number of cell types downloaded: ${cellTypeCount}.\n\n"

#=================================== APPRIS ==================================================

#Downloading APPRIS data:
mkdir -p ${targetDir}/${today}/APPRIS
info "Downloading APPRIS isoform data.\n"
info "Download from the current release folder. Build: GRCh38, for GENCODE version: ${GENCODE_release}\n"
wget -q http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt \
    -O ${targetDir}/${today}/APPRIS/appris_data.principal.txt

# Testing if the file exists or not:
testFile "${targetDir}/${today}/APPRIS/appris_data.principal.txt"

info "Download complete.\n\n"

#=============================================================================================
#Combining APPRIS and GENCODE data

info "Combining APPRIS and GENCODE data.. "
mkdir -p ${targetDir}/${today}/processed
export APPRIS_FILE=${targetDir}/${today}/APPRIS/appris_data.principal.txt
zcat ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | grep -v "#" | awk '$3 != "Selenocysteine" && $3 != "start_codon" && $3 != "stop_codon"' \
                | perl -MJSON -M"Data::Dumper"  -F"\t" -lane '
                BEGIN {
                    $af = $ENV{APPRIS_FILE};
                    open($APP, "<", $af);
                    while ($line = <$APP>){
                        chomp $line;
                        $line =~ /(ENSG\S+)\s+(ENST\S+)\s+\S+\s+(\S+):/;
                        $h{$1}{$2} = $3;
                    }
                }{
                    if ( $_ =~ /(ENSG.+?)\./){ # Gene related annotation
                        $geneID = $1;

                        $F[0] =~ s/chr//;
                        $start = $F[3];
                        $end = $F[4];
                        $strand = $F[6];
                        $class = $F[2];

                        ($transcriptID) = $F[8] =~ /(ENST.+?)\./ ? $F[8] =~ /(ENST.+?)\./ : "NA";
                        ($exonID) = $F[8] =~ /(ENSE.+?)\./ ? $F[8] =~ /(ENSE.+?)\./ : "NA";

                        $appris = "NA";
                        # Check if current feature belongs to an APPRIS annotated feature:
                        if( exists $h{$geneID} && $transcriptID ne "NA" ){
                            if (exists $h{$geneID}{$transcriptID}){
                                $appris = $h{$geneID}{$transcriptID};
                            }
                            else {
                                $appris = "Minor";
                            }
                        }

                        # Saving output in json format:
                        %hash = (
                            "chr" => $F[0],
                            "start" => $start,
                            "end" => $end,
                            "source" => "GENCODE",
                            "strand" => $strand,
                            "class" => $class,
                            "gene_ID" => $geneID,
                            "appris" => $appris
                        );
                        $hash{"transcript_ID"} = $transcriptID if $transcriptID ne "NA";
                        $hash{"exon_ID"} = $exonID if $exonID ne "NA";
                        print JSON->new->utf8->encode(\%hash);
                    }
                }' | gzip > ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz

# Test if output is empty or not:
testFileLines  ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz

# OUTPUT:
#
#{"source":"GENCODE","gene_ID":"ENSG00000186092","appris":"NA","start":"65419","chr":"1","strand":"+","class":"gene","end":"71585"}
#{"transcript_ID":"ENST00000641515","end":"71585","class":"transcript","strand":"+","appris":"Minor","start":"65419","chr":"1","source":"GENCODE","gene_ID":"ENSG00000186092"}
#{"strand":"+","start":"65419","chr":"1","appris":"Minor","gene_ID":"ENSG00000186092","exon_ID":"ENSE00003812156","source":"GENCODE","transcript_ID":"ENST00000641515","end":"65433","class":"exon"}
#{"appris":"Minor","chr":"1","start":"65520","strand":"+","source":"GENCODE","exon_ID":"ENSE00003813641","gene_ID":"ENSG00000186092","transcript_ID":"ENST00000641515","end":"65573","class":"exon"}
#{"class":"CDS","transcript_ID":"ENST00000641515","end":"65573","gene_ID":"ENSG00000186092","exon_ID":"ENSE00003813641","source":"GENCODE","strand":"+","start":"65565","chr":"1","appris":"Minor"}
#{"gene_ID":"ENSG00000186092","source":"GENCODE","exon_ID":"ENSE00003813949","chr":"1","start":"69037","appris":"Minor","strand":"+","class":"exon","transcript_ID":"ENST00000641515","end":"71585"}
#{"end":"70005","transcript_ID":"ENST00000641515","class":"CDS","strand":"+","appris":"Minor","start":"69037","chr":"1","exon_ID":"ENSE00003813949","source":"GENCODE","gene_ID":"ENSG00000186092"}
#{"gene_ID":"ENSG00000186092","exon_ID":"ENSE00003812156","source":"GENCODE","strand":"+","start":"65419","chr":"1","appris":"Minor","class":"UTR","end":"65433","transcript_ID":"ENST00000641515"}
#{"gene_ID":"ENSG00000186092","exon_ID":"ENSE00003813641","source":"GENCODE","strand":"+","start":"65520","chr":"1","appris":"Minor","class":"UTR","end":"65564","transcript_ID":"ENST00000641515"}
#{"end":"71585","transcript_ID":"ENST00000641515","class":"UTR","strand":"+","chr":"1","start":"70006","appris":"Minor","gene_ID":"ENSG00000186092","exon_ID":"ENSE00003813949","source":"GENCODE"}
#{"source":"GENCODE","gene_ID":"ENSG00000186092","strand":"+","appris":"PRINCIPAL","chr":"1","start":"69055","class":"transcript","transcript_ID":"ENST00000335137","end":"70108"}
#{"appris":"PRINCIPAL","chr":"1","start":"69055","strand":"+","source":"GENCODE","exon_ID":"ENSE00002319515","gene_ID":"ENSG00000186092","transcript_ID":"ENST00000335137","end":"70108","class":"exon"}
#{"appris":"PRINCIPAL","chr":"1","start":"69091","strand":"+","source":"GENCODE","exon_ID":"ENSE00002319515","gene_ID":"ENSG00000186092","transcript_ID":"ENST00000335137","end":"70005","class":"CDS"}
#{"class":"UTR","end":"69090","transcript_ID":"ENST00000335137","gene_ID":"ENSG00000186092","source":"GENCODE","exon_ID":"ENSE00002319515","chr":"1","start":"69055","appris":"PRINCIPAL","strand":"+"}
#{"source":"GENCODE","exon_ID":"ENSE00002319515","gene_ID":"ENSG00000186092","appris":"PRINCIPAL","chr":"1","start":"70006","strand":"+","class":"UTR","transcript_ID":"ENST00000335137","end":"70108"}

echo "Done."

# Print out report:
appris_lines=$(zcat ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz | wc -l | awk '{print $1}')
info "Number of Appris annotated GENCODE annotations: ${appris_lines}\n\n"

#=============================================================================================

##
## Step 7. Pre-processing cell specific regulatory data
##
info "Aggregate cell specific information of regulatory features... "
#CellTypes=$( ls -la ${targetDir}/${today}/EnsemblRegulation/ | perl -lane 'print $1 if  $F[-1] =~ /RegulatoryFeatures_(.+).gff.gz/ ' )
CellTypes=$( ls -la ${targetDir}/${today}/EnsemblRegulation/ | perl -lane 'print $1 if  $F[-1] =~ /(.+).gff.gz/ ' )
for cell in ${CellTypes}; do
    export cell
    # parsing cell specific files (At this point we only consider active features. Although repressed regions might also be informative.):
    zcat ${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz | grep -i "activity=active" \
        | perl -F"\t" -lane 'next unless length($F[0]) < 3 || $F[0]=~/^chr/; # We skip irregular chromosome names.
                $F[0]=~s/^chr//;
                $cell_type = $ENV{cell};
                $start = $F[3];
                $type = $F[2];
                $end = $F[4];
                ($ID) = $_ =~ /regulatory_feature_stable_id=(ENSR\d+)/;
                ($bstart) = $F[8] =~ /bound_start=(.+?);/;
                ($bend) = $F[8] =~ /bound_end=(.+?);/;
                print join "\t", $cell_type, $F[0], $start-1, $end, $ID, $type, $bstart, $bend;' # $start - 1 as bed is 0-based
# Now combining these lines in a fashion that each line will contain all active tissues:
done | perl -F"\t" -lane '
    $x =shift @F;
    $h{$F[3]}{line} = [@F];
    push(@{$h{$F[3]}{cells}}, $x);
    END {
        foreach $ID (keys %h){
            $cells = join "|", @{$h{$ID}{cells}};
            printf "%s\t%s\t%s\t%s\tchr=%s;start=%s;end=%s;class=%s;regulatory_ID=%s;Tissues=%s\n",
                $h{$ID}{line}[0], $h{$ID}{line}[1], $h{$ID}{line}[2], $h{$ID}{line}[3], $h{$ID}{line}[0],
                $h{$ID}{line}[1], $h{$ID}{line}[2], $h{$ID}{line}[4], $h{$ID}{line}[3], $cells
        }
    }
' | sort -k1,1 -k2,2n | bgzip -f > ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz


# OUTPUT:
#
#1    9800    10400   ENSR00000344264 chr=1;start=9800;end=10400;class=CTCF_binding_site;regulatory_ID=ENSR00000344264;Tissues=HUVEC|HeLa_S3|mammary_epithelial_1
#1    13400   13600   ENSR00000344265 chr=1;start=13400;end=13600;class=CTCF_binding_site;regulatory_ID=ENSR00000344265;Tissues=HUVEC|NHLF|keratinocyte
#1    15400   16600   ENSR00000344266 chr=1;start=15400;end=16600;class=CTCF_binding_site;regulatory_ID=ENSR00000344266;Tissues=A549|B|H1_hESC_3|HCT116|HSMM|HeLa_S3|HepG2|K562|MCF_7|MM_1S|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|mammary_epithelial_1|osteoblast
#1    16102   16451   ENSR00000918273 chr=1;start=16102;end=16451;class=TF_binding_site;regulatory_ID=ENSR00000918273;Tissues=DND_41|K562|PC_3|bipolar_neuron
#1    24600   25000   ENSR00000344268 chr=1;start=24600;end=25000;class=CTCF_binding_site;regulatory_ID=ENSR00000344268;Tissues=A549|H1_hESC_3|HepG2|MCF_7
#1    25600   26400   ENSR00000344269 chr=1;start=25600;end=26400;class=CTCF_binding_site;regulatory_ID=ENSR00000344269;Tissues=A549|A673|B|GM12878|H1_hESC_3|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MCF_7|MM_1S|NHLF|PC_3|astrocyte|cardiac_muscle|dermal_fibroblast|mammary_epithelial_1|osteoblast


# Test if output is empty or not:
testFileLines ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz

tabix -f -p bed ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz
echo  "Done."

# Print out report:
cellSpecFeatLines=$(zcat ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz | wc -l | awk '{print $1}')
info "Number of cell specific regulatory features: $cellSpecFeatLines\n\n"


##
## Step 8. Adding GRCh38 coordinates to GTEx data. (based on rsID)
##

# Instead of the single step we can generate a bedfile and run liftover
# This step takes around 7 minutes.
# info "Mapping GTEx variants to GRCh38 build.\n"
# info "Creating temporary bed file (~9 minutes)... "
# zcat ${GTExFile}  | perl -F"\t" -lane '
#         if ($_ =~ /snpgenes/){
#             ($tissue) = $_ =~ /([A-Z]+.+)_Analysis.snpgenes/;
#             next;
#         }
#         ($chr, $pos, $ref, $alt, $build) = split("_", $F[0]);
#         ($gene) = $F[1] =~ /(ENS.+)\./;
#         $rsID = $F[22];

#         $h{$rsID}{chr}= $chr;
#         $h{$rsID}{pos}= $pos;
#         push( @{$h{$rsID}{genes}{$gene}}, $tissue ) if $tissue;

#         END {
#             foreach $rsID ( keys %h){
#                 $chr = $h{$rsID}{chr};
#                 $pos = $h{$rsID}{pos};

#                 foreach $gene ( keys %{$h{$rsID}{genes}}){
#                     $tissues = join "|", @{$h{$rsID}{genes}{$gene}};

#                     # Reporting problem if something comes upon:

#                     printf "chr$chr\t%s\t$pos\tgene=$gene;rsID=$rsID;tissue=$tissues\n", $pos - 1 if $chr and $pos;
#                 }
#             }
#         }
#     '  | sort -k1,1 -k2,2n > ${targetDir}/${today}/processed/GTEx_temp.bed

# # Testing if output file has lines:
# testFileLines ${targetDir}/${today}/processed/GTEx_temp.bed

# echo "Done."

# info "Running liftOver (~2 minutes).... "
# liftOver ${targetDir}/${today}/processed/GTEx_temp.bed ${scriptDir}/hg19ToHg38.over.chain \
#     ${targetDir}/${today}/processed/GTEx_temp_GRCh38.bed \
#     ${targetDir}/${today}/processed/GTEx_temp_failed_to_map.bed
# echo "Done."

# # Generate report:
# failedMap=$(wc -l ${targetDir}/${today}/processed/GTEx_temp_failed_to_map.bed | awk '{print $1}')
# Mapped=$(wc -l ${targetDir}/${today}/processed/GTEx_temp_GRCh38.bed | awk '{print $1}')
# info "Successfully mapped GTEx variants: ${Mapped}, failed variants: ${failedMap}.\n\n"

#=============================================================================================

##
## Step 8. Combine individual files from GTEx tar.gz file into one BED file
##

tmpGTEx=${targetDir}/${today}/processed/GTEx_tmp.bed
info "Creating GTEx bed file ... "
listOfGTExFiles=$(tar -ztf ${GTExFile} | grep "signif_variant")
for f in ${listOfGTExFiles};do
    g=$(basename ${f})
    tissue=$(echo ${g}|perl -lne '$x="NA";if (/^([^.]+)\./){$x=$1;} print $x;')
    export tissue
    tar -zxf ${GTExFile} ${f} -O | zcat - | tail -n +2 | perl -F"\t" -lane '($chr, $pos, $ref, $alt, $build) = split("_", $F[0]);($gene) = $F[1] =~ /(ENS.+)\./;$tissue=$ENV{tissue};$,="\t";$chr=~s/^chr//;print $tissue,$chr,$pos,$F[0],$gene;'
done > ${tmpGTEx}

cat ${tmpGTEx} | perl -F"\t" -lane '$tissue=$F[0];$chr=$F[1];$pos=$F[2];$ID=$F[3];$gene=$F[4];$H{$ID}{chr}=$chr;$H{$ID}{pos}=$pos;push( @{$H{$ID}{genes}{$gene}}, $tissue ); END {foreach $id (keys %H){$chr=$H{$id}{chr};$pos=$H{$id}{pos};foreach $gene (keys %{$H{$id}{genes}}){$tissues = join "|", @{$H{$id}{genes}{$gene}};printf "$chr\t%s\t$pos\tgene=$gene;rsID=$id;tissue=$tissues\n", $pos - 1;}}}' | sort -k1,1 -k2,2 > ${targetDir}/${today}/processed/GTEx.bed

echo "Done."
rm -f ${tmpGTEx}

# OUTPUT:
#
#1    10000043        10000044        gene=ENSG00000054523;rsID=chr1_10000044_A_T_b38;tissue=Adipose_Visceral_Omentum|Artery_Tibial|Breast_Mammary_Tissue|Cells_Cultured_fibroblasts|Nerve_Tibial|Skin_Not_Sun_Exposed_Suprapubic|Skin_Sun_Exposed_Lower_leg
#1    10000043        10000044        gene=ENSG00000130939;rsID=chr1_10000044_A_T_b38;tissue=Esophagus_Mucosa
#1    10000043        10000044        gene=ENSG00000175279;rsID=chr1_10000044_A_T_b38;tissue=Artery_Aorta|Artery_Tibial|Esophagus_Mucosa|Lung|Skin_Not_Sun_Exposed_Suprapubic
#1    100000722       100000723       gene=ENSG00000117620;rsID=chr1_100000723_G_A_b38;tissue=Brain_Cortex|Cells_Cultured_fibroblasts|Lung|Nerve_Tibial|Skin_Not_Sun_Exposed_Suprapubic|Skin_Sun_Exposed_Lower_leg|Stomach|Testis
#1    100000722       100000723       gene=ENSG00000122435;rsID=chr1_100000723_G_A_b38;tissue=Skin_Sun_Exposed_Lower_leg
#1    100000722       100000723       gene=ENSG00000122477;rsID=chr1_100000723_G_A_b38;tissue=Adipose_Subcutaneous|Adipose_Visceral_Omentum|Brain_Cerebellar_Hemisphere|Brain_Cerebellum|Breast_Mammary_Tissue|Colon_Transverse|Esophagus_Muscularis|Liver|Thyroid


#=============================================================================================

##
## Step 9. Using intersectbed. Find overlap between GTEx variations and regulatory regions
##
info "Linking genes to regulatory features using GTEx data... "
intersectBed -wb -a ${targetDir}/${today}/processed/GTEx.bed -b ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz \
    | perl -MData::Dumper -MJSON -F"\t" -lane '
        # Name of the source is GTEx
        $source= "GTEx";

# OUTPUT OF INTERSECTBED
#
#1    100034322       100034323       gene=ENSG00000122435;rsID=chr1_100034323_TG_T_b38;tissue=Muscle_Skeletal        1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034322       100034323       gene=ENSG00000122477;rsID=chr1_100034323_TG_T_b38;tissue=Adipose_Subcutaneous|Adipose_Visceral_Omentum|Artery_Aorta|Artery_Tibial|Brain_Cerebellar_Hemisphere|Brain_Cerebellum|Brain_Cortex|Breast_Mammary_Tissue|Cells_Cultured_fibroblasts|Colon_Transverse|Esophagus_Gastroesophageal_Junction|Esophagus_Muscularis|Liver|Lung|Nerve_Tibial|Ovary|Pancreas|Skin_Not_Sun_Exposed_Suprapubic|Skin_Sun_Exposed_Lower_leg|Spleen|Stomach|Thyroid|Whole_Blood     1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034326       100034327       gene=ENSG00000122435;rsID=chr1_100034327_T_A_b38;tissue=Muscle_Skeletal 1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034326       100034327       gene=ENSG00000122477;rsID=chr1_100034327_T_A_b38;tissue=Adipose_Subcutaneous|Adipose_Visceral_Omentum|Artery_Aorta|Artery_Tibial|Brain_Cerebellar_Hemisphere|Brain_Cerebellum|Brain_Cortex|Breast_Mammary_Tissue|Cells_Cultured_fibroblasts|Colon_Transverse|Esophagus_Gastroesophageal_Junction|Esophagus_Muscularis|Liver|Lung|Nerve_Tibial|Ovary|Pancreas|Skin_Not_Sun_Exposed_Suprapubic|Skin_Sun_Exposed_Lower_leg|Spleen|Stomach|Thyroid|Whole_Blood      1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034329       100034330       gene=ENSG00000122435;rsID=chr1_100034330_A_G_b38;tissue=Muscle_Skeletal 1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034329       100034330       gene=ENSG00000122477;rsID=chr1_100034330_A_G_b38;tissue=Adipose_Subcutaneous|Adipose_Visceral_Omentum|Artery_Aorta|Artery_Tibial|Brain_Cerebellar_Hemisphere|Brain_Cerebellum|Brain_Cortex|Breast_Mammary_Tissue|Cells_Cultured_fibroblasts|Colon_Transverse|Esophagus_Gastroesophageal_Junction|Esophagus_Muscularis|Liver|Lung|Nerve_Tibial|Ovary|Pancreas|Skin_Not_Sun_Exposed_Suprapubic|Skin_Sun_Exposed_Lower_leg|Spleen|Stomach|Thyroid|Whole_Blood      1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast


        # Parsing input:
        ($gene) = $F[3] =~ /gene=(ENSG.+?);/;
        ($G_rsID) = $F[3] =~ /rsID=(.+?);/;
        ($G_tissues) = $F[3] =~ /tissue=(\S+)/;
        $E_chr = $F[4];
        $E_start = $F[5];
        $E_end = $F[6];
        $E_ID = $F[7];
        ($E_class) = $F[8] =~ /class=(.+?);/;
        ($E_tissues) = $F[8] =~ /Tissues=(\S+)/;

        # Building hash:
        $h{$gene."_".$E_ID}{gene_ID} = $gene;
        $h{$gene."_".$E_ID}{class} = $E_class;
        $h{$gene."_".$E_ID}{source} = $source;
        $h{$gene."_".$E_ID}{chr} = $E_chr;
        $h{$gene."_".$E_ID}{start} = $E_start;
        $h{$gene."_".$E_ID}{end} = $E_end;
        $h{$gene."_".$E_ID}{regulatory_ID} = $E_ID;
        $h{$gene."_".$E_ID}{Tissues} = [split /\|/, $E_tissues];

        # Adding GTEx details:
        push(@{$h{$gene."_".$E_ID}{GTEx_rsIDs}}, $G_rsID);
        push(@{$h{$gene."_".$E_ID}{GTEx_tissues}}, (split /\|/, $G_tissues));

        # Saving results when reading data has finished:
        END {
            # Looping through all gene/reg feature pairs:
            for $key ( keys %h){
                # Looping through all GTEx tissues and keep only the unique ones.
                %a = ();
                foreach $tissue (@{$h{$key}{GTEx_tissues}}){
                    $a{$tissue} = 1;
                }
                $h{$key}{GTEx_tissues} = [keys %a];

                # Saving json:
                print JSON->new->utf8->encode($h{$key})
            }
        }
    ' | gzip > ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz

# OUTPUT
# we don't know which of the rsIDs corresponds to which GTEx tissue
#{"chr":"chr2","GTEx_rsIDs":[...],"gene_ID":"ENSG00000235584","GTEx_tissues":["Thyroid","Lung"],"class":"open_chromatin_region","Tissues":["MM_1S"],"source":"GTEx","start":"96131175","regulatory_ID":"ENSR00000613314","end":"96131863"}

# Testing if output file has lines:
testFileLines ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz
echo "Done."


# Generate report:
GTExLinkedFeatures=$( zcat ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz | wc -l | awk '{print $1}')
info "Number of GTEx linked regulatory features: ${GTExLinkedFeatures}\n\n"


#=============================================================================================


##
## Step 10. Using intersectbed. Find overlap between genes and regulatory regions
##
info "Linking genes to regulatory features based on overlap... "
# generating a file.
zcat ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | awk '$3 == "gene"' | perl -lane '
        ($g_name) = $_ =~ /gene_name "(.+?)";/;
        ($g_ID) = $_ =~ /gene_id "(.+?)\.*";/;
        $F[0]=~s/^chr//;
        printf "$F[0]\t$F[3]\t$F[4]\tID:$g_ID;Name:$g_name\n";
    ' | sort -k1,1 -k2,2n | bgzip -f > ${targetDir}/${today}/processed/genes.bed.gz

# OUTPUT:
#
#1    11869   14409   ID:ENSG00000223972;Name:DDX11L1


# Intersect bed run.
# 1	16048	29570	ID:ENSG00000227232;Name:WASH7P	1	16048	30847	ENSR00000528774	chr=1;start=16048;end=30847;class=CTCF_binding_site;regulatory_ID=ENSR00000528774;Tissues=DND-41|HMEC|HSMMtube|IMR90|K562|MultiCell|NHDF-AD
intersectBed -wb -a ${targetDir}/${today}/processed/genes.bed.gz -b ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz -sorted \
    | perl -MData::Dumper -MJSON -F"\t" -lane '
        # Parsing gene info:
        ($g_ID) = $F[3] =~ /ID:(ENSG\d+?)/;
        ($g_name) = $F[3] =~ /Name:(\S+)/;

        # Parsing regulatory feature info:
        ($r_chr, $r_start, $r_end, $r_ID) = ($F[4], $F[5], $F[6], $F[7]);
        ($r_class) = $F[8] =~ /class=(.+?);/;
        ($r_tissue_string) = $F[8] =~ /Tissues=(.+)/;
        @r_tissues = split(/\|/, $r_tissue_string);

        # Saving JSON formatted string:
        print JSON->new->utf8->encode({
            "chr"   => $r_chr,
            "start" => $r_start,
            "end"   => $r_end,
            "class" => $r_class,
            "gene_ID"   => $g_ID,
            "gene_name" => $g_name,
            "Tissues"   => \@r_tissues,
            "regulatory_ID" => $r_ID,
            "source"    => "overlap",
        })
    ' | bgzip -f > ${targetDir}/${today}/processed/overlapping_features.txt.gz
echo "Done."

# OUTPUT:
#
#{"start":"13400","chr":"1","gene_name":"DDX11L1","gene_ID":"ENSG00000223972","regulatory_ID":"ENSR00000344265","class":"CTCF_binding_site","source":"overlap","Tissues":["HUVEC","NHLF","keratinocyte"],"end":"13600"}

# Generate report:
OverlapLinkedFeatures=$( zcat ${targetDir}/${today}/processed/overlapping_features.txt.gz | wc -l | awk '{print $1}')
info "Number of regulatory features linked by overlap: ${OverlapLinkedFeatures}\n\n"


#=============================================================================================

##
## Step 11. Merging all the components together create compressed, sorted bedfile.
##
info "Merging GENCODE, GTEx and overlap data together into an indexed bedfile. "
export gene_file=${targetDir}/${today}/processed/genes.bed.gz # make sure file readable from within the perl script

#GENES
#
#1    11869   14409   ID:ENSG00000223972.5;Name:DDX11L1

#INPUT:
#
#{"start":"13400","chr":"chr1","gene_name":"DDX11L1","gene_ID":"ENSG00000223972","regulatory_ID":"ENSR00000344265","class":"CTCF_binding_site","source":"overlap","Tissues":["HUVEC","NHLF","keratinocyte"],"end":"13600"}
#{"chr":"chr2","GTEx_rsIDs":[null],"gene_ID":"ENSG00000235584","GTEx_tissues":["Thyroid","Lung"],"class":"open_chromatin_region","Tissues":["MM_1S"],"source":"GTEx","start":"96131175","regulatory_ID":"ENSR00000613314","end":"96131863"}
#{"appris":"NA","start":"11869","chr":"chr1","strand":"+","source":"GENCODE","gene_ID":"ENSG00000223972","end":"14409","class":"gene"}


zcat ${targetDir}/${today}/processed/overlapping_features.txt.gz \
     ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz \
     ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz  \
     | perl -lane 'BEGIN {
            open  $GF, "zcat $ENV{gene_file} |";
            while ($line = <$GF>){
                chop $line;
                @a = split "\t", $line;
                ($ID) = $a[3] =~ /ID:(ENSG\d+)/;
                $h{$ID} = [$a[0], $a[1], $a[2], $ID];

            }
        }{

            ($ID) = $_ =~ /"gene_ID":"(ENSG\d+)"/;
            exists $h{$ID} ? print join "\t", @{$h{$ID}}, $_ : print STDERR "$ID : gene was notfound in gencode! line: $_"
        }'  2> ${targetDir}/${today}/failed | sort -k1,1 -k2,2n > ${targetDir}/${today}/Linked_features.bed

echo -e "Done.\n"

# Creating header for the final output:
cat <(echo -e "# Regions file for burden testing. Created: ${today}
#
# GENCODE version: v.${GENCODE_release}
# Ensembl version: v.${Ensembl_release}
# GTEx version: v.${GTExRelease}
#
# CHR\tSTART\tEND\tGENEID\tANNOTATION" ) ${targetDir}/${today}/Linked_features.bed | sponge ${targetDir}/${today}/Linked_features.bed

# Compressing and indexing output file:
bgzip -f ${targetDir}/${today}/Linked_features.bed > ${targetDir}/${today}/Linked_features.bed.gz
tabix -f -p bed ${targetDir}/${today}/Linked_features.bed.gz

# Final report and we are done.
info "Output file was saved as: ${targetDir}/${today}/Linked_features.bed.gz\n"
totalLines=$(zcat ${targetDir}/${today}/Linked_features.bed.gz | wc -l | awk '{print $1}')
info "Total number of lines in the final files: ${totalLines}\n"

# Report failed associations:
FailedAssoc=$(wc -l ${targetDir}/${today}/failed | awk '{print $1}')
FailedGenes=$( cat ${targetDir}/${today}/failed | perl -lane '$_ =~ /(ENSG\d+)/; print $1' | sort | uniq | wc -l )
FailedSources=$( cat ${targetDir}/${today}/failed | perl -lane '$_ =~ /"source":"(.+?)"/; print $1' | sort | uniq | tr "\n" ", " )
info "Number of lost associations: ${FailedAssoc}, belonging to ${FailedGenes} genes in the following sournces: ${FailedSources}\n\n"

# Backing up intermediate files:
#tar czf ${targetDir}/${today}/${today}_annotation.backup.tar.gz --remove-file   ${targetDir}/${today}/APPRIS  \
#    ${targetDir}/${today}/EnsemblRegulation  ${targetDir}/${today}/failed  ${targetDir}/${today}/GENCODE  \
#    ${targetDir}/${today}/processed
#
#info "Intermediate files are backed in in ${targetDir}/${today}/${today}_annotation.backup.tar.gz\n" 

# Exit.
info "Program finished.\n"
exit 0
