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
    echo "  bgzip, tabix in path"
    echo "  liftOver in path"
    echo "  hg19ToHg38.over.chain chain file in script dir"
    echo "  bedtools in path"
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

# Get the most recent version of the data:
mkdir -p ${targetDir}/${today}/GENCODE
info "Downloading GENCODE annotation from ftp://ftp.ebi.ac.uk/. Release version: ${GENCODE_release}... "
wget -q ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_release}/gencode.v${GENCODE_release}.annotation.gtf.gz \
        -O ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz
echo -e "done."

# Testing if the file is exists or not:
testFile "${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz"

# Counting genes in the dataset:
genes=$(zcat ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | awk 'BEGIN{FS="\t";}$3 == "gene"{print $3;}' | wc -l )
info "Total number of genes in the GENCODE file: ${genes}\n\n"

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

# Download all cell types:
for cell in ${cells}; do
    echo -n "."
#    echo $cell
    # Download all cell type:
    wget -q ftp://ftp.ensembl.org/pub/release-${Ensembl_release}/regulation/homo_sapiens/RegulatoryFeatureActivity/${cell}/homo_sapiens.*Regulatory_Build.regulatory_activity.*.gff.gz \
        -O ${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz

    # Testing if the file is exists or not:
    testFile "${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz"

done
echo "Done."

# Printing out report of the downloaded cell types:
cellTypeCount=$(ls -la ${targetDir}/${today}/EnsemblRegulation/*gff.gz | wc -l)
info "Number of cell types downloaded: ${cellTypeCount}.\n\n"

# Downloading APPRIS data:
mkdir -p ${targetDir}/${today}/APPRIS
info "Downloading APPRIS isoform data.\n"
info "Download from the current release folder. Build: GRCh38, for GENCODE version: ${GENCODE_release}\n"
wget -q http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt \
    -O ${targetDir}/${today}/APPRIS/appris_data.principal.txt

# Testing if the file is exists or not:
testFile "${targetDir}/${today}/APPRIS/appris_data.principal.txt"

info "Download complete.\n\n"

## Combining APPRIS and GENCODE data
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
                        $line =~ /(ENSG.+?)\s+(ENST.+?)\s.+?\s+(\S+)\:/;
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
                        # Check if the given feature is belong to an annotated feature:
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
                            "chr" => "chr".$F[0],
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

echo "Done."

# Print out report:
appris_lines=$(zcat ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz | wc -l | awk '{print $1}')
info "Number of Appris annotated GENCODE annotations: ${appris_lines}\n\n"

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
        | perl -F"\t" -lane 'next unless length($F[0]) < 3; # We skip irregular chromosome names.
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
            printf "chr%s\t%s\t%s\t%s\tchr=%s;start=%s;end=%s;class=%s;regulatory_ID=%s;Tissues=%s\n",
                $h{$ID}{line}[0], $h{$ID}{line}[1], $h{$ID}{line}[2], $h{$ID}{line}[3], $h{$ID}{line}[0],
                $h{$ID}{line}[1], $h{$ID}{line}[2], $h{$ID}{line}[4], $h{$ID}{line}[3], $cells
        }
    }
' | sort -k1,1 -k2,2n | bgzip -f > ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz

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
info "Mapping GTEx variants to GRCh38 build.\n"
info "Creating temporary bed file (~9 minutes)... "
zcat ${GTExFile}  | perl -F"\t" -lane '
        if ($_ =~ /snpgenes/){
            ($tissue) = $_ =~ /([A-Z]+.+)_Analysis.snpgenes/;
            next;
        }
        ($chr, $pos, $ref, $alt, $build) = split("_", $F[0]);
        ($gene) = $F[1] =~ /(ENS.+)\./;
        $rsID = $F[22];

        $h{$rsID}{chr}= $chr;
        $h{$rsID}{pos}= $pos;
        push( @{$h{$rsID}{genes}{$gene}}, $tissue ) if $tissue;

        END {
            foreach $rsID ( keys %h){
                $chr = $h{$rsID}{chr};
                $pos = $h{$rsID}{pos};

                foreach $gene ( keys %{$h{$rsID}{genes}}){
                    $tissues = join "|", @{$h{$rsID}{genes}{$gene}};

                    # Reporting problem if something comes upon:

                    printf "chr$chr\t%s\t$pos\tgene=$gene;rsID=$rsID;tissue=$tissues\n", $pos - 1 if $chr and $pos;
                }
            }
        }
    '  | sort -k1,1 -k2,2n > ${targetDir}/${today}/processed/GTEx_temp.bed

# Testing if output file has lines:
testFileLines ${targetDir}/${today}/processed/GTEx_temp.bed

echo "Done."

info "Running liftOver (~2 minutes).... "
liftOver ${targetDir}/${today}/processed/GTEx_temp.bed ${scriptDir}/hg19ToHg38.over.chain \
    ${targetDir}/${today}/processed/GTEx_temp_GRCh38.bed \
    ${targetDir}/${today}/processed/GTEx_temp_failed_to_map.bed
echo "Done."

# Generate report:
failedMap=$(wc -l ${targetDir}/${today}/processed/GTEx_temp_failed_to_map.bed | awk '{print $1}')
Mapped=$(wc -l ${targetDir}/${today}/processed/GTEx_temp_GRCh38.bed | awk '{print $1}')
info "Successfully mapped GTEx variants: ${Mapped}, failed variants: ${failedMap}.\n\n"

##
## Step 9. Using intersectbed. Find overlap between GTEx variations and regulatory regions
##
info "Linking genes to regulatory features using GTEx data... "
intersectBed -wb -a ${targetDir}/${today}/processed/GTEx_temp_GRCh38.bed -b ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz \
    | perl -MData::Dumper -MJSON -F"\t" -lane '
        # Name of the source is GTEx
        $source= "GTEx";

        # Parsing input:
        ($gene) = $F[3] =~ /gene=(ENSG.+?);/;
        ($G_rsID) = $F[3] =~ /rsID=(rs.+?);/;
        ($G_tissues) = $F[3] =~ /tissue=(.+)/;
        $E_chr = $F[4];
        $E_start = $F[5];
        $E_end = $F[6];
        $E_ID = $F[7];
        ($E_class) = $F[8] =~ /class=(.+?);/;
        ($E_tissues) = $F[8] =~ /Tissues=(.+)/;

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

# Testing if output file has lines:
testFileLines ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz

echo "Done."

# Generate report:
GTExLinkedFeatures=$( zcat ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz | wc -l | awk '{print $1}')
info "Number of GTEx linked regulatory features: ${GTExLinkedFeatures}\n\n"

##
## Step 10. Using intersectbed. Find overlap between genes and regulatory regions
##
info "Linking genes to regulatory features based on overlap... "
# generating a file.
zcat ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | awk '$3 == "gene"' | perl -lane '
        ($g_name) = $_ =~ /gene_name "(.+?)";/;
        ($g_ID) = $_ =~ /gene_id "(.+?)";/;
        printf "$F[0]\t$F[3]\t$F[4]\tID:$g_ID;Name:$g_name\n";
    ' | sort -k1,1 -k2,2n | bgzip -f > ${targetDir}/${today}/processed/genes.bed.gz

# Intersect bed run.
# chr1	16048	29570	ID:ENSG00000227232.5;Name:WASH7P	chr1	16048	30847	ENSR00000528774	chr=1;start=16048;end=30847;class=CTCF_binding_site;regulatory_ID=ENSR00000528774;Tissues=DND-41|HMEC|HSMMtube|IMR90|K562|MultiCell|NHDF-AD
intersectBed -wb -a ${targetDir}/${today}/processed/genes.bed.gz -b ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz -sorted \
    | perl -MData::Dumper -MJSON -F"\t" -lane '
        # Parsing gene info:
        ($g_ID) = $F[3] =~ /ID:(ENSG\d+)/;
        ($g_name) = $F[3] =~ /Name:(.+)/;

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

# Generate report:
OverlapLinkedFeatures=$( zcat ${targetDir}/${today}/processed/overlapping_features.txt.gz | wc -l | awk '{print $1}')
info "Number of regulatory features linked by overlap: ${OverlapLinkedFeatures}\n\n"

##
## Step 11. Merging all the components together create compressed, sorted bedfile.
##
info "Merging GENCODE, GTEx and overlap data together into an indexed bedfile. "
export gene_file=${targetDir}/${today}/processed/genes.bed.gz # make sure file readable from within the perl script
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
# GTEx version: ${GTExRelease}
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

# Backing up intermedier files:
tar czf ${targetDir}/${today}/${today}_annotation.backup.tar.gz --remove-file   ${targetDir}/${today}/APPRIS  \
    ${targetDir}/${today}/EnsemblRegulation  ${targetDir}/${today}/failed  ${targetDir}/${today}/GENCODE  \
    ${targetDir}/${today}/GTEx  ${targetDir}/${today}/processed

info "Intermedier files are backed in in ${targetDir}/${today}/${today}_annotation.backup.tar.gz\n" 

# Exit.
info "Program finished.\n"
exit 0
