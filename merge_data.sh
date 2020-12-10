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

# This merges data downloaded using \"download_data.sh\" script

script_version=1.0
last_modified=2020.Dec.10
today=$(date "+%Y.%m.%d")

# Get script dir:
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## printing out information if no parameter is provided:
function usage {
    echo ""
    echo "Usage: $0 -i <input directory> : required, input directory containing data downloaded using \"download_data.sh\""
    echo "          -x : optional, create backup of the downloaded data"
    echo ""
    echo " This script was written to prepare input file for the burden testing pipeline."
    echo ""
    echo ""
    echo "Version: ${script_version}, Last modified: ${last_modified}"
    echo ""
    echo ""
    echo "Workflow:"
    echo "  1: Adds Appris annotation to Gencode transcripts"
    echo "  2: Creates cell-specific regulatory features"
    echo "  3: Links regulatory features to genes based on GTEx data"
    echo "  4: Links regulatory features to genes based on overlapping"
    echo "  5: Combined GENCODE, GTEx and Overlap data together into a single bedfile"
    echo "  6: Tabix output"
    echo ""
    echo ""
    echo "This script produces two output files, both in the same directory as the input GTEx file."
    echo "1) the first file  is \"Linked_features.bed.gz\"; its first 4 columns are the chromosome, start/end
coordinates and the stable ID of the gene respectively. The 5th column is a json
formatted string describing one genomic region associated with the given gene. This
line contains all information of the association."
    echo ""
    echo "2) the second file is \"gencode.basic.annotation.tsv.gz\"; it's a trimmed version of the downloaded GENCODE file."
    echo ""
    echo "JSON tags:"
    echo "  -source: from which source the given region is coming from (GENCODE, GTEx, Overlap)."
    echo "  -class: class of the given feature (eg. exon, CDS, gene, enhancer, promoter etc.)"
    echo "  -chr, start, end: GRCh38 coordintes of the feature."
    echo "  -other information. (linked rsID, tissue in which the feature in active etc.)"
    echo ""
    echo "IMPORTANT: ALL COORDINATES ARE BASED ON GRCh38 BUILD"
    echo ""
    exit 0
}

function checkdir {
    echo -n "Checking if directory $1 exists ... "
    if [[ ! -d "$1"  ]]; then
        echo "[Error] Directory does not exist: $1"
        echo "[Error] Exit"
        exit 1
    else
	echo "OK"
    fi
}

function checkfile {
    echo -n "Checking if file $1 exists ... "
    if [[ ! -f "$1"  ]]; then
        echo "[Error] File does not exist: $1"
        echo "[Error] Exit"
        exit 1
    else
	echo "OK"
    fi
}

# check if GZ file is OK
function checkGZfile {
    echo -n "Checking GZ file integrity: $1 ... "
    if ! gzip -q -t "$1";then
	echo "[Error] Integrity check failed for $1"
        echo "[Error] Exit"
        exit 1
    else
	echo "OK"
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
        echo "[Error] File ($1) contains no lines. Exit.";
        exit 1;
    fi
}

# This function prints out all the reports that were generated during the run (with time stamp!):
function info {
    hourMin=$(date +"%T" | awk 'BEGIN{FS=OFS=":"}{print $1, $2}')
    echo -e "[Info ${hourMin}] $1"
}

# Printing help message if no parameters are given:
if [[ $# == 0 ]]; then usage; fi

# Processing command line options:
OPTIND=1
indir=""
backup="no"
while getopts "hi:x" optname; do
    case "$optname" in
        "h" ) usage ;;
        "i" ) indir="${OPTARG}" ;;
        "x" ) backup="yes" ;;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ -z ${indir} ]];then
    echo "[Error] no input directory specified"
    exit 1
fi

# full dirname
indir=`readlink -f $indir`
indir=${indir%/}

# output config file
configfile=${indir}/config.txt
rm -f ${configfile}

ENSVEPDIR="${indir}/ensembl-vep"
VEPDIR="${indir}/vep"
targetDir="${indir}/downloads"
GENCODEDIR="${targetDir}/GENCODE"
ENSDIR="${targetDir}/ENSEMBL"
APPRISDIR="${targetDir}/APPRIS"
GTEXDIR="${targetDir}/GTEx"

GENCODEFILE="${GENCODEDIR}/gencode.annotation.gtf.gz"
GTEXFILE="${GTEXDIR}/GTEx.tar.gz"
APPRISFILE="${APPRISDIR}/appris_data.principal.txt"

checkdir ${ENSVEPDIR}
checkdir ${VEPDIR}
checkdir ${targetDir}
checkdir ${GENCODEDIR}
checkdir ${ENSDIR}
checkdir ${APPRISDIR}
checkdir ${GTEXDIR}

checkfile ${GENCODEFILE}
checkfile ${GTEXFILE}
checkfile ${APPRISFILE}

PROCDIR="${targetDir}/processed"
mkdir -p ${PROCDIR}
checkdir ${PROCDIR}

#===================================== VEP ===================================================

cd ${ENSVEPDIR}
PATH=$PATH:${VEPDIR}/htslib PERL5LIB=$PERL5LIB:${VEPDIR} perl INSTALL.pl -a ac -n --ASSEMBLY GRCh38 -s homo_sapiens -c ${VEPDIR} -d ${VEPDIR}
echo "VEPdir=${VEPDIR}" >>  ${configfile}
echo "VEPexec=${ENSVEPDIR}/vep" >>  ${configfile}
cd ${indir}

#=================================== GENCODE =================================================

# Get the most recent version of the data:
genes=$(zcat ${GENCODEFILE} | awk 'BEGIN{FS="\t";}$3 == "gene"{print $3;}' | wc -l )
info "Total number of genes in the GENCODE file: ${genes}\n"

#=================================== REGULATION ==============================================

# Printing out report of the downloaded cell types:
cellTypeCount=$(ls -la ${ENSDIR}/*gff.gz | wc -l)
info "Number of downloaded cell types: ${cellTypeCount}\n"

#=============================================================================================
#Combining APPRIS and GENCODE data

info "Combining APPRIS and GENCODE data.. "
export APPRIS_FILE=${APPRISFILE}
zcat ${GENCODEFILE} | grep -v "#" | awk '$3 != "Selenocysteine" && $3 != "start_codon" && $3 != "stop_codon"' \
                | perl -MJSON -M"Data::Dumper"  -F"\t" -lane '
                BEGIN {
                    $af = $ENV{APPRIS_FILE};
                    open($APP, "<", $af);
                    while ($line = <$APP>){
                        chomp $line;
                        $line =~ /(ENSG\S+)\s+(ENST\S+)\s+\S+\s+(\S+):/;
                        $geneID=$1;
                        $transcriptID=$2;
                        $status=$3;
                        $h{$geneID}{$transcriptID} = $status;
                    }
                }{
                    if ( $_ =~ /gene_id\s+\"(ENSG.+?)\"/){ # Gene related annotation
                        $geneID = $1;
                        # if the gene is from PAR region on Y then we keep the full gene ID and transcript ID 
                        # otherwise, we remove the suffix
                        if ($geneID !~ /_PAR_Y/){
                        $geneID =~ /(ENSG\d+)\./;
                        $geneID=$1;
                        }

			$transcriptID="NA";
			if ($F[8] =~ /transcript_id\s+\"(ENST.+?)\"/){
			$tmp=$1;
			if ($tmp !~ /_PAR_Y/){
			$tmp =~ /(ENST\d+)\./;
			$transcriptID=$1;
			}
			else{
			$transcriptID=$tmp;
			}
			}

			$exonID="NA";
			if ($F[8] =~ /exon_id\s+\"(ENSE.+?)\./){
			$exonID=$1;
			}

                        $F[0] =~ s/chr//;
                        $start = $F[3];
                        $end = $F[4];
                        $strand = $F[6];
                        $class = $F[2];

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
                        # start/end are zero based
                        %hash = (
                            "chr" => $F[0],
                            "start" => $start-1,
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
                }' | gzip > ${PROCDIR}/Appris_annotation_added.txt.gz

# Test if output is empty:
checkfile ${PROCDIR}/Appris_annotation_added.txt.gz
testFileLines ${PROCDIR}/Appris_annotation_added.txt.gz # 0-based

# OUTPUT:
# start/end are 0-based coordinates of "class"
#
#{"source":"GENCODE","gene_ID":"ENSG00000186092","appris":"NA","start":"65419","chr":"1","strand":"+","class":"gene","end":"71585"}
#{"transcript_ID":"ENST00000641515","end":"71585","class":"transcript","strand":"+","appris":"Minor","start":"65419","chr":"1","source":"GENCODE","gene_ID":"ENSG00000186092"}
#{"strand":"+","start":"65419","chr":"1","appris":"Minor","gene_ID":"ENSG00000186092","exon_ID":"ENSE00003812156","source":"GENCODE","transcript_ID":"ENST00000641515","end":"65433","class":"exon"}

echo "Done"

# Print out report:
appris_lines=$(zcat ${PROCDIR}/Appris_annotation_added.txt.gz | wc -l | awk '{print $1}')
info "Number of Appris annotated GENCODE annotations: ${appris_lines}\n"

#=============================================================================================

##
## Step 7. Pre-processing cell specific regulatory data
##
info "Aggregating cell specific information of regulatory features... "
#CellTypes=$( ls -la ${targetDir}/${today}/EnsemblRegulation/ | perl -lane 'print $1 if  $F[-1] =~ /RegulatoryFeatures_(.+).gff.gz/ ' )
CellTypes=$( ls -la ${ENSDIR} | perl -lane 'print $1 if  $F[-1] =~ /(.+).gff.gz/ ' )
for cell in ${CellTypes}; do
    export cell
    fn=${ENSDIR}/${cell}.gff.gz
    
    # parsing cell specific files (At this point we only consider active features. Although repressed regions might also be informative):

    # Check integrity
    checkGZfile ${fn}
    
    zcat ${fn} | grep -i "activity=active" \
        | perl -F"\t" -lane 'next unless length($F[0]) < 3 || $F[0]=~/^chr/; # We skip irregular chromosome names.
                $F[0]=~s/^chr//;
                $cell_type = $ENV{cell};
                $start = $F[3];

		# $start - 1 as bed is 0-based
		$start=$start-1;

                $type = $F[2];
                $end = $F[4];
                ($ID) = $_ =~ /regulatory_feature_stable_id=(ENSR\d+)/;
                print join "\t", $cell_type, $F[0], $start, $end, $ID, $type;' 
# Now combining these lines in a fashion that each line will contain all active tissues:
done | perl -F"\t" -lane '
    $x =shift @F;
    $h{$F[3]}{line} = [@F];
    push(@{$h{$F[3]}{cells}}, $x);
    END {
        foreach $ID (keys %h){
            $cells = join "|", @{$h{$ID}{cells}};
            printf "%s\t%s\t%s\t%s\tchr=%s;start=%s;end=%s;class=%s;regulatory_ID=%s;Tissues=%s\n",$h{$ID}{line}[0], $h{$ID}{line}[1], $h{$ID}{line}[2], $ID, $h{$ID}{line}[0],$h{$ID}{line}[1], $h{$ID}{line}[2], $h{$ID}{line}[4], $ID, $cells;
        }
    }
' | sort -k1,1 -k2,2n | bgzip -f > ${PROCDIR}/Cell_spec_regulatory_features.bed.gz # 0-based coordinates here


# OUTPUT:
#
#1    9800    10400   ENSR00000344264 chr=1;start=9800;end=10400;class=CTCF_binding_site;regulatory_ID=ENSR00000344264;Tissues=HUVEC|HeLa_S3|mammary_epithelial_1
#1    13400   13600   ENSR00000344265 chr=1;start=13400;end=13600;class=CTCF_binding_site;regulatory_ID=ENSR00000344265;Tissues=HUVEC|NHLF|keratinocyte

# Test if output is empty:
checkfile ${PROCDIR}/Cell_spec_regulatory_features.bed.gz
testFileLines ${PROCDIR}/Cell_spec_regulatory_features.bed.gz
tabix -f -p bed ${targetDir}/processed/Cell_spec_regulatory_features.bed.gz
echo  "Done"

# Print out report:
cellSpecFeatLines=$(zcat ${PROCDIR}/Cell_spec_regulatory_features.bed.gz | wc -l | awk '{print $1}')
info "Number of cell specific regulatory features: $cellSpecFeatLines\n"

#=============================================================================================

##
## Step 8. Combine individual files from GTEx tar.gz file into one BED file
##

# no PAR_Y IDs in GTEx (v.8) files, so we don't check

tmpGTEx=${PROCDIR}/GTEx_tmp.bed
info "Creating GTEx bed file ... "
listOfGTExFiles=$(tar -ztf ${GTEXFILE} | grep "signif_variant")
for f in ${listOfGTExFiles};do
    g=$(basename ${f})
    tissue=$(echo ${g}|perl -lne '$x="NA";if (/^([^.]+)\./){$x=$1;} print $x;')
    export tissue

    # taking care of deletions as well
    tar -zxf ${GTEXFILE} ${f} -O | zcat - | tail -n +2 | perl -F"\t" -lane '($chr, $pos, $ref, $alt, $build) = split("_", $F[0]);($gene) = $F[1] =~ /(ENSG\d+)\./;$tissue=$ENV{tissue};$,="\t";$chr=~s/^chr//;$start=$pos-1;$end=$pos;if (length($ref)>length($alt)){$end=$start+length($ref)-1;}  print $tissue,$chr,$start,$end,$F[0],$gene;'
done > ${tmpGTEx}

cat ${tmpGTEx} | perl -F"\t" -lane '$tissue=$F[0];$chr=$F[1];$start=$F[2];$end=$F[3];$ID=$F[4];$gene=$F[5];$H{$ID}{chr}=$chr;$H{$ID}{start}=$start;$H{$ID}{end}=$end;push( @{$H{$ID}{genes}{$gene}}, $tissue ); END {foreach $id (keys %H){$chr=$H{$id}{chr};$start=$H{$id}{start};$end=$H{$id}{end};foreach $gene (keys %{$H{$id}{genes}}){$tissues = join "|", @{$H{$id}{genes}{$gene}};print "$chr\t$start\t$end\tgene=$gene;rsID=$id;tissue=$tissues";}}}' | sort -k1,1 -k2,2n > ${PROCDIR}/GTEx.bed # 0-based

echo "Done"
checkfile ${PROCDIR}/GTEx.bed
testFileLines ${PROCDIR}/GTEx.bed
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
## Step 9. Using intersectBed. Overlap between GTEx variants and regulatory regions
##
info "Linking genes to regulatory features using GTEx data ... "

# OUTPUT OF INTERSECTBED BELOW
#
#1    100034322       100034323       gene=ENSG00000122435;rsID=chr1_100034323_TG_T_b38;tissue=Muscle_Skeletal        1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034322       100034323       gene=ENSG00000122477;rsID=chr1_100034323_TG_T_b38;tissue=Adipose_Subcutaneous|Adipose_Visceral_Omentum|Artery_Aorta|Artery_Tibial|Brain_Cerebellar_Hemisphere|Brain_Cerebellum|Brain_Cortex|Breast_Mammary_Tissue|Cells_Cultured_fibroblasts|Colon_Transverse|Esophagus_Gastroesophageal_Junction|Esophagus_Muscularis|Liver|Lung|Nerve_Tibial|Ovary|Pancreas|Skin_Not_Sun_Exposed_Suprapubic|Skin_Sun_Exposed_Lower_leg|Spleen|Stomach|Thyroid|Whole_Blood     1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034326       100034327       gene=ENSG00000122435;rsID=chr1_100034327_T_A_b38;tissue=Muscle_Skeletal 1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034326       100034327       gene=ENSG00000122477;rsID=chr1_100034327_T_A_b38;tissue=Adipose_Subcutaneous|Adipose_Visceral_Omentum|Artery_Aorta|Artery_Tibial|Brain_Cerebellar_Hemisphere|Brain_Cerebellum|Brain_Cortex|Breast_Mammary_Tissue|Cells_Cultured_fibroblasts|Colon_Transverse|Esophagus_Gastroesophageal_Junction|Esophagus_Muscularis|Liver|Lung|Nerve_Tibial|Ovary|Pancreas|Skin_Not_Sun_Exposed_Suprapubic|Skin_Sun_Exposed_Lower_leg|Spleen|Stomach|Thyroid|Whole_Blood      1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034329       100034330       gene=ENSG00000122435;rsID=chr1_100034330_A_G_b38;tissue=Muscle_Skeletal 1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast
#1    100034329       100034330       gene=ENSG00000122477;rsID=chr1_100034330_A_G_b38;tissue=Adipose_Subcutaneous|Adipose_Visceral_Omentum|Artery_Aorta|Artery_Tibial|Brain_Cerebellar_Hemisphere|Brain_Cerebellum|Brain_Cortex|Breast_Mammary_Tissue|Cells_Cultured_fibroblasts|Colon_Transverse|Esophagus_Gastroesophageal_Junction|Esophagus_Muscularis|Liver|Lung|Nerve_Tibial|Ovary|Pancreas|Skin_Not_Sun_Exposed_Suprapubic|Skin_Sun_Exposed_Lower_leg|Spleen|Stomach|Thyroid|Whole_Blood      1    100034200       100035200       ENSR00000253249 chr=1;start=100034200;end=100035200;class=CTCF_binding_site;regulatory_ID=ENSR00000253249;Tissues=A549|A673|B|CD14_monocyte_1|DND_41|GM12878|HCT116|HSMM|HUVEC|HeLa_S3|HepG2|K562|MM_1S|NHLF|PC_3|PC_9|SK_N_|astrocyte|bipolar_neuron|cardiac_muscle|dermal_fibroblast|keratinocyte|myotube|osteoblast

intersectBed -wb -a ${PROCDIR}/GTEx.bed -b ${PROCDIR}/Cell_spec_regulatory_features.bed.gz 2>/dev/null | perl -MData::Dumper -MJSON -F"\t" -lane '
        $source= "GTEx";

        ($gene) = $F[3] =~ /gene=(ENSG\d+)/;
        ($G_rsID) = $F[3] =~ /rsID=(.+?);/;
        ($G_tissues) = $F[3] =~ /tissue=(\S+)/;
        $E_chr = $F[4];
        $E_start = $F[5];
        $E_end = $F[6];
        $E_ID = $F[7];
        ($E_class) = $F[8] =~ /class=(.+?);/;
        ($E_tissues) = $F[8] =~ /Tissues=(\S+)/;

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

        END {
            # Looping through all gene/reg feature pairs:
            for $key ( keys %h){
                # Looping through all GTEx tissues and keep only the unique ones.
                %a = ();
                foreach $tissue (@{$h{$key}{GTEx_tissues}}){
                    $a{$tissue} = 1;
                }
                $h{$key}{GTEx_tissues} = [keys %a];

                print JSON->new->utf8->encode($h{$key})
            }
        }
    ' | gzip > ${PROCDIR}/GTEx_Regulation_linked.txt.gz

# OUTPUT
# we don't know which of the rsIDs corresponds to which GTEx tissue
#{"chr":"chr2","GTEx_rsIDs":[...],"gene_ID":"ENSG00000235584","GTEx_tissues":["Thyroid","Lung"],"class":"open_chromatin_region","Tissues":["MM_1S"],"source":"GTEx","start":"96131175","regulatory_ID":"ENSR00000613314","end":"96131863"}

# Testing if output file has lines:
checkfile ${PROCDIR}/GTEx_Regulation_linked.txt.gz
testFileLines ${PROCDIR}/GTEx_Regulation_linked.txt.gz # start/end are 0-based
echo "Done"

# Generate report:
GTExLinkedFeatures=$( zcat ${PROCDIR}/GTEx_Regulation_linked.txt.gz | wc -l | awk '{print $1}')
info "Number of GTEx linked regulatory features: ${GTExLinkedFeatures}\n"

#=============================================================================================


##
## Step 10. Using intersectbed. Find overlap between genes and regulatory regions
##
info "Linking genes to regulatory features based on overlap ... "
# generating a file.
# 0-based
zcat ${GENCODEFILE} | awk '$3 == "gene"' | perl -lane '
     $g_name="NA";
        if ($_ =~ /gene_name\s+\"(.+?)\";/){
	$g_name=$1;
	}

	$g_ID="NA";
        if ($_ =~ /gene_id\s+\"(.+?)\";/){
	$g_ID=$1;
	}

	# if gene ID contains PAR_Y then keep it
	# otherwise remove the suffix

	if ($g_ID !~ /_PAR_Y/){
	 $g_ID =~ /(ENSG\d+)\./;
	 $g_ID=$1;
	}

        $F[0]=~s/^chr//;
        $start=$F[3]-1;
        print "$F[0]\t$start\t$F[4]\tID:$g_ID;Name:$g_name";
    ' | sort -k1,1 -k2,2n | bgzip -f > ${PROCDIR}/genes.bed.gz # 0-based

# OUTPUT:
#
#1    11868   14409   ID:ENSG00000223972;Name:DDX11L1

checkfile ${PROCDIR}/genes.bed.gz
testFileLines ${PROCDIR}/genes.bed.gz

# Intersect bed output:
# 1	16048	29570	ID:ENSG00000227232;Name:WASH7P	1	16048	30847	ENSR00000528774	chr=1;start=16048;end=30847;class=CTCF_binding_site;regulatory_ID=ENSR00000528774;Tissues=DND-41|HMEC|HSMMtube|IMR90|K562|MultiCell|NHDF-AD
intersectBed -wb -a ${PROCDIR}/genes.bed.gz -b ${PROCDIR}/Cell_spec_regulatory_features.bed.gz -sorted 2>/dev/null | perl -MData::Dumper -MJSON -F"\t" -lane '
        # Parsing gene info:
        ($g_ID) = $F[3] =~ /ID:(ENSG[^;]+);/;
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
    ' | bgzip -f > ${PROCDIR}/overlapping_features.txt.gz # 0-based coordinates
echo "Done"

checkfile ${PROCDIR}/overlapping_features.txt.gz
testFileLines ${PROCDIR}/overlapping_features.txt.gz

# OUTPUT:
#
#{"start":"13400","chr":"1","gene_name":"DDX11L1","gene_ID":"ENSG00000223972","regulatory_ID":"ENSR00000344265","class":"CTCF_binding_site","source":"overlap","Tissues":["HUVEC","NHLF","keratinocyte"],"end":"13600"}

# Generate report:
OverlapLinkedFeatures=$( zcat ${PROCDIR}/overlapping_features.txt.gz | wc -l | awk '{print $1}')
info "Number of regulatory features linked by overlap: ${OverlapLinkedFeatures}\n"

#=============================================================================================

##
## Step 11. Merging all the components together create sorted, compressed bedfile.
##
info "Merging GENCODE, GTEx and overlap data"
export gene_file=${PROCDIR}/genes.bed.gz

#GENES
#
#1    11869   14409   ID:ENSG00000223972.5;Name:DDX11L1

#INPUT:
#
#{"start":"13400","chr":"chr1","gene_name":"DDX11L1","gene_ID":"ENSG00000223972","regulatory_ID":"ENSR00000344265","class":"CTCF_binding_site","source":"overlap","Tissues":["HUVEC","NHLF","keratinocyte"],"end":"13600"}
#{"chr":"chr2","GTEx_rsIDs":[...],"gene_ID":"ENSG00000235584","GTEx_tissues":["Thyroid","Lung"],"class":"open_chromatin_region","Tissues":["MM_1S"],"source":"GTEx","start":"96131175","regulatory_ID":"ENSR00000613314","end":"96131863"}
#{"appris":"NA","start":"11869","chr":"chr1","strand":"+","source":"GENCODE","gene_ID":"ENSG00000223972","end":"14409","class":"gene"}


zcat ${PROCDIR}/overlapping_features.txt.gz \
     ${PROCDIR}/GTEx_Regulation_linked.txt.gz \
     ${PROCDIR}/Appris_annotation_added.txt.gz  \
     | perl -lane 'BEGIN {
            open  $GF, "zcat $ENV{gene_file} |";
            while ($line = <$GF>){
                chop $line;
                @a = split "\t", $line;
                ($ID) = $a[3] =~ /ID:(ENSG[^;]+);/;
                $h{$ID} = [$a[0], $a[1], $a[2], $ID];
            }
        }{
            ($ID) = $_ =~ /\"gene_ID\":\"(ENSG[^"]+)\"/;
            exists $h{$ID} ? print join "\t", @{$h{$ID}}, $_ : print STDERR "$ID : gene was not found in gencode! line: $_"
        }'  2> ${indir}/failed | sort -k1,1 -k2,2n > ${indir}/Linked_features.bed # 0-based

echo -e " Done.\n"

# source == GENCODE => chr,start,end in the 5th field are those of transcript,gene,exon,UTR,CDS
# source == GTEx => chr,start,end in the 5th field are those of regulatory element
# source == overlap => chr,start,end in the 5th field are those of regulatory element

# Creating header for the final output:
cat <(echo -e "# Regions file for burden testing. Created: ${today}
# CHR\tSTART\tEND\tGENEID\tANNOTATION" ) ${indir}/Linked_features.bed | sponge ${indir}/Linked_features.bed

# Compressing and indexing output file:
bgzip -f ${indir}/Linked_features.bed > ${indir}/Linked_features.bed.gz
tabix -f -p bed ${indir}/Linked_features.bed.gz

info "Output file was saved as: ${indir}/Linked_features.bed.gz\n"
totalLines=$(zcat ${indir}/Linked_features.bed.gz | wc -l | awk '{print $1}')
info "Total number of lines in the final file: ${totalLines}\n"

#=============================================================================================

# GENCODE basic annotation, IDs are not modified
# coordinates are 1-based

zcat  ${GENCODEFILE} | grep -v "^#"| perl -F"\t" -lane 'next if $F[2] ne "gene";$x=$F[8];$id="NA";$id=$1 if ($x=~/gene_id\s+\"(ENSG[^"]+)\"/); $gn="NA"; $gn=$1 if $x=~/gene_name\s+\"([^"]+)\"/;$,="\t";$F[0]=~s/^chr//;print $F[0],$F[3],$F[4],$gn,$id;' | gzip > ${indir}/gencode.basic.annotation.tsv.gz
zcat  ${GENCODEFILE} | grep -v "^#"| perl -F"\t" -lane 'next if $F[2] ne "gene";$x=$F[8];next unless $x=~/gene_type\s+\"protein_coding\"/;$id="NA";$id=$1 if ($x=~/gene_id\s+\"(ENSG[^"]+)\"/); $gn="NA"; $gn=$1 if $x=~/gene_name\s+\"([^"]+)\"/;$,="\t";$F[0]=~s/^chr//;print $F[0],$F[3],$F[4],$gn,$id;' | gzip > ${indir}/gencode.basic.annotation.protein_coding.tsv.gz

#==================================== OUTPUT config.txt ======================================

echo "Linked_features=${indir}/Linked_features.bed.gz" >> ${configfile}
echo "gencode_file=${indir}/gencode.basic.annotation.tsv.gz" >> ${configfile}

#=============================================================================================

# Report failed associations:
FailedAssoc=$(wc -l ${indir}/failed | awk '{print $1}')
FailedGenes=$( cat ${indir}/failed | perl -lane '$_ =~ /(ENSG\d+)/; print $1' | sort | uniq | wc -l )
FailedSources=$( cat ${indir}/failed | perl -lane '$_ =~ /"source":"(.+?)"/; print $1' | sort | uniq | tr "\n" ", " )
info "Number of lost associations: ${FailedAssoc}, belonging to ${FailedGenes} genes in the following sournces: ${FailedSources}\n"

if [[ $backup == "yes" ]];then
    info "Backing up intermediate files ...\n"
    tar czf ${targetDir}/${today}_annotation.backup.tar.gz --remove-file ${APPRISDIR} ${ENSDIR} ${indir}/failed ${GENCODEDIR} ${PROCDIR}
    info "Intermediate files are saved in ${targetDir}/${today}_annotation.backup.tar.gz\n"
else
    rm -rf ${targetDir}
fi

if [[ -f "${indir}/scores/eigen.phred_v2.dat" ]];then
    echo "EigenPath=${indir}/scores/eigen.phred_v2.dat" >> ${configfile}
fi
if [[ -f "${indir}/scores/whole_genome_SNVs.tsv.gz" ]];then
    echo "caddPath=${indir}/scores/whole_genome_SNVs.tsv.gz" >> ${configfile}
fi

info "Merging finished\n"
exit 0
