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
## 4. GTEx (v.8)

# The GTEx eQTL data (v.8) will be downloaded

# For reproducibility, both the GENCODE, Ensembl and GTEx versions are hardcoded


script_version=4.1.1

## Built in versions:
GENCODE_release=32
Ensembl_release=98
GTExRelease=8

# Get script dir:
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## printing out information if no parameter is provided:
function usage {
    echo ""
    echo "Usage: $0 -o <output directory> : required, output directory"
    echo "          -n : optional, do not download Eigen scores"
    echo "          -c : optional, do not download CADD scores"
    echo "          -x : optional, create backup of the downloaded data"
    echo "          -r : optional, re-use previous downloads if present"
    echo "          -d : optional, just download data and exit, do not process them"
    echo "          -s : optional, do not perform checksums (unsafe)"
    echo "          -t : optional, directory to store temporary data, default: \"/tmp\""
    echo "          -m : optional, debug mode"
    echo "          -e <ftp://ftp.ensembl.org> : optional, Ensembl FTP server"
    echo ""
    echo " This script was written to prepare input file for the burden testing pipeline."
    echo ""
    echo ""
    echo "Version: ${script_version}, Last modified: ${last_modified}"
    echo ""
    echo "  downloaded GTEx datafile with the single eQTLs (eg. GTEx_Analysis_v8_eQTLs.tar.gz)"
    echo ""
    echo ""
    echo "Workflow:"
    echo "  0. Downloads GTEx (v.8) dataset."
    echo "  1: Downloads v32 GENCODE release."
    echo "  2: Downloads V97 Ensembl Regulation release."
    echo "  3: Downloads newest APPRIS release"
    echo "  4: Adds Appris annotation to Gencode transcripts."
    echo "  5: Creates cell-specific regulatory features."
    echo "  6: Links regulatory features to genes based on GTEx data."
    echo "  7: Links regulatory features to genes based on overlapping."
    echo "  8: Combined GENCODE, GTEx and Overlap data together into a single bedfile."
    echo "  9: Tabix output"
    echo "  10: optionally download Egen Phred scores (if \"-n\" is not sepcified)"
    echo "  11: optionally download CADD scores (if \"-c\" is not sepcified)"
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
    echo "Optionally, a \"scores\" folder will be created and Eigen Phred scores will be saved in it"
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

# Function to test if a given file exists or not in which case it reports and terminates the
# execution.
function testFile {
    if [[ ! -e "$1"  ]]; then
        echo "[Error] At this step something failed. The file was not created! $1"
        echo "[Error] Exiting."
        exit 1
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

# check if GZ file is OK
function checkGZfile {
    # echo -n "Checking GZ file integrity: $1 ... "
    if ! gzip -q -t "$1";then
    echo "[Error] Integrity check failed for $1"
        echo "[Error] Exit"
        exit 1
    # else
    # echo "OK"
    fi
}

# This function prints out all the reports that were generated during the run (with time stamp!):
function info {
    hourMin=$(date +"%T" | awk 'BEGIN{FS=OFS=":"}{print $1, $2}')
    echo -e "[Info ${hourMin}] $1"
}

function debug {
    if [[ "$debug_mode" -eq 1 ]] ; then
      hourMin=$(date +"%T" | awk 'BEGIN{FS=OFS=":"}{print $1, $2}')
      echo -e "[Debug ${hourMin}] $1"
    fi
}

# Printing help message if no parameters are given:
if [[ $# == 0 ]]; then usage; fi

# Processing command line options:
getScores="yes"
ensftp="ftp://ftp.ensembl.org"
OPTIND=1
outdir=""
getCadd="yes"
backup="no"
reuse=0
justdl=0
noSums=0
debug_mode=0
tempdir="/tmp"
while getopts "hncxst:e:rdmo:" optname; do
    case "$optname" in
        "h" ) usage ;;
        "n" ) getScores="no" ;;
        "c" ) getCadd="no" ;;
        "x" ) backup="yes" ;;
        "s" ) noSums=1 ;;
        "t" ) tempdir="${OPTARG}" ;;
        "e" ) ensftp="${OPTARG}" ;;
        "r" ) reuse=1 ;;
        "d" ) justdl=1 ;;
        "m" ) debug_mode=1 ;;
        "o" ) outdir="${OPTARG}" ;;
        "?" ) usage ;;
        *) usage ;;
    esac;
done

if [[ -z ${outdir} ]];then
    echo "[Error] no output directory specified"
    exit 1
fi

mkdir -p ${outdir}
if [ $? -ne 0 ] ; then
    echo "[Error] Could not create ${outdir}"
    exit 1
fi

# full dirname
outdir=`realpath $outdir`
outdir=${outdir%/}

# output config file
configfile=${outdir}/config.txt
if [ "$reuse" -eq 0 ]; then
  rm -f ${configfile}
else
  if [[ -s "${configfile}" ]]; then
    info "Previous config file found at ${configfile}."
  else
    info "Reuse flag is set but no or empty config file found at ${configfile}. Starting from scratch."
  fi
fi

# temp dir
if [ ! -d ${tempdir} ]; then
    mkdir -p ${tempdir}
    if [ $? -ne 0 ] ; then
    echo "[Error] Could not create temp dir ${tempdir}"
    exit 1
    fi
fi
tempdir=`realpath $tempdir`

# report arguments
echo ""
echo "Provided arguments:"
echo ""
echo "Output directory         : $outdir"
echo "Temp directory           : $tempdir"
echo "Get Eigen scores         : $getScores"
echo "Get CADD scores          : $getCadd"
echo "Make backup              : $backup"
echo "Do not perform checksums : $noSums"
echo "Ensembl FTP server       : $ensftp"
echo "Re-use previous downloads: $reuse"
echo "Download only            : $justdl"
echo "Debug mode               : $debug_mode"
echo ""

## Write a metadata.txt file which contains info about the executed run
metafile=${outdir}/metadata.txt
_script_loc="`dirname \"$0\"`"
# We use --git-dir instead of -C, because -C isn't available in older versions of git
commit_id=$(git --git-dir ${_script_loc}/.git rev-parse HEAD)
if [ -z "$commit_id" ] ; then
  info "Commit ID couldn't be retrieved. Falling back to an unknown value."
  commit_id='???'
fi
echo "repo-commit-id: $commit_id" > $metafile
echo "prepare-regions-ver: $script_version" >> $metafile
echo "GENCODE_release: $GENCODE_release" >> $metafile
echo "Ensembl_release: $Ensembl_release" >> $metafile
echo "GTEx_release: $GTExRelease" >> $metafile
echo "Arguments: 
  - outdir: $outdir
  - tempdir: $tempdir
  - getScores: $getScores
  - getCadd: $getCadd
  - backup: $backup
  - noSums: $noSums
  - ensftp: $ensftp
  - reuse: $reuse
  - justdl: $justdl
  - debug_mode: $debug_mode
" >> $metafile
echo "Started: $(date +'%Y-%m-%d %H:%M:%S')" >> $metafile

#===================================== VEP ===================================================

if (( "$reuse" > 0 )) && [[ ! -z "$(grep VEPdir ${configfile})" ]]; then
  info "VEP information found in config file. Skipping VEP download (unsafe - no checks performed)..."
else
  cd ${outdir}

  git clone https://github.com/Ensembl/ensembl-vep.git
  cd ${outdir}/ensembl-vep
  git checkout release/$Ensembl_release

  mkdir -p ${outdir}/vep \
    && cd ${outdir}/vep \
    && axel -a ftp://ftp.ensembl.org/pub/release-${Ensembl_release}/variation/indexed_vep_cache/homo_sapiens_vep_${Ensembl_release}_GRCh38.tar.gz \
    && echo Unpacking ... \
    && tar -xzf homo_sapiens_vep_${Ensembl_release}_GRCh38.tar.gz \
    && rm homo_sapiens_vep_${Ensembl_release}_GRCh38.tar.gz \
    && cd ${outdir}/ensembl-vep

  sed 's/ensembl\.org/ebi\.ac\.uk\/ensemblorg/g' INSTALL.pl | sponge INSTALL.pl

  PATH=$PATH:${outdir}/vep/htslib PERL5LIB=$PERL5LIB:${outdir}/vep perl INSTALL.pl -a ac -n --ASSEMBLY GRCh38 -s homo_sapiens -c ${outdir}/vep -d ${outdir}/vep

  echo "VEPdir=${outdir}/vep" >>  ${configfile}
  echo "VEPexec=${outdir}/ensembl-vep/vep" >>  ${configfile}
fi

#=============================================================================================

cd ${outdir}
targetDir=${outdir}"/prepare_regions_tempfiles"
mkdir -p ${targetDir}
if [ $? -ne 0 ] ; then
    echo "[Error] Could not create ${targetDir}"
    exit 1
fi

GTExFile=$outdir/GTEx_Analysis_v8_eQTL.tar
#echo $reuse $noSums $GTExFile
#if [[ -s "$GTExFile" ]];then echo lol; fi
#exit 0
if (( "$reuse" > 0 )) \
  && [[ -s "$GTExFile" ]] \
  && [[ "$noSums" == "1" || $(md5sum $GTExFile | cut -d' ' -f1) == "d35b32152bdb21316b2509c46b0af998" ]]; then
    info "GTEx file found and has the right checksum. Skipping download..."
else
    cd ${outdir}
    axel -a https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar
    # gzip -f GTEx_Analysis_v8_eQTL.tar

    if [[ ! -e "${GTExFile}" ]]; then
        echo "[Error] GTEx file (${GTExFile}) does not exist."
        exit 1
    fi
fi

if [[ "$noSums" == "0" && $(md5sum $GTExFile | cut -d' ' -f1) != "d35b32152bdb21316b2509c46b0af998" ]]; then
    echo "[Error] Checksum invalid ($(md5sum $GTExFile | cut -d' ' -f1)). The download probably failed. Please rerun with the reuse option (-r) to retry."
    rm $GTExFile
    exit 1
fi

# Last step in setup:
today=$(date "+%Y.%m.%d")
info "Current date: ${today}\n"
info "Working directory: ${targetDir}/${today}\n\n"

#=================================== GENCODE =================================================

# Get the most recent version of the data:
mkdir -p ${targetDir}/${today}/GENCODE

info "Getting MD5 hash of gencode.v${GENCODE_release}.annotation.gtf.gz file."
checksum=$(wget -q -O- http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_release}/MD5SUMS | grep -w gencode.v${GENCODE_release}.annotation.gtf.gz | cut -d' ' -f1)

if (( "$reuse" > 0 )) \
  && [[ ! -z "$(find $targetDir -name gencode.v${GENCODE_release}.annotation.gtf.gz | head -1)" ]] \
  && [[ "$noSums" == "1" || "$(md5sum $(find $targetDir -name gencode.v${GENCODE_release}.annotation.gtf.gz | head -1) | cut -d' ' -f1)" == "$checksum" ]] ; then
    info "GENCODE file found and has the right checksum. Skipping download..."
    fn=$(find $targetDir -name gencode.v${GENCODE_release}.annotation.gtf.gz | head -1)
    if [[ ! "$fn" -ef ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz ]];then
        debug 'mv "$fn" ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz'
        mv "$fn" ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz
    fi
else
    info "Downloading GENCODE annotation. Release version: ${GENCODE_release}... "
    gencode_output_file="${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz"
    axel -a ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_release}/gencode.v${GENCODE_release}.annotation.gtf.gz \
      -o $gencode_output_file
    echo -e "done."

    if [[ ! -f "$gencode_output_file" ]] ; then
        info "GENCODE download with FTP failed. Retrying with HTTP"
        axel -a http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_release}/gencode.v${GENCODE_release}.annotation.gtf.gz \
          -o $gencode_output_file
        echo -e "done."
    fi

    # Testing if the file exists:
    testFile "$gencode_output_file"

    if [[ "$noSums" == "0" && "$(md5sum $gencode_output_file | cut -d' ' -f1)" != "$checksum" ]]; then
        echo "[Error] Checksum invalid ($(md5sum $gencode_output_file | cut -d' ' -f1)). The download probably failed. Please rerun with the reuse option (-r) to retry."
        rm $gencode_output_file
        exit 1
    fi
fi
# Counting genes in the dataset:
genes=$(zcat ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | awk 'BEGIN{FS="\t";}$3 == "gene"{print $3;}' | wc -l )
info "Total number of genes in the GENCODE file: ${genes}\n\n"

#=================================== REGULATION ==============================================

# prepare target directory:
mkdir -p ${targetDir}/${today}/EnsemblRegulation
info "Downloading cell specific regulatory features from Ensembl.\n"

# Get list of all cell types:
localCellFile=$(find $targetDir -name RegCellList.txt)
if [[ ! -z "$localCellFile" && "$reuse" == "1" ]]; then
    echo Found cell type file at $localCellFile with $(cat $localCellFile | wc -l) types.
else
    localCellFile=$targetDir/$today/RegCellList.txt
    curl -s ${ensftp}/pub/release-${Ensembl_release}/regulation/homo_sapiens/RegulatoryFeatureActivity/ | perl -lane 'print $F[-1]' > $localCellFile
fi
cells=$(cat $localCellFile)


# If there are no cell types present in the downloaded set, it means there were some problems. We are exiting.
if [ -z "${cells}" ]; then
    echo "[Error] No cell types were found in the Ensembl regulation folder."
    echo "[Error] URL: ${ensftp}/pub/release-${Ensembl_release}/regulation/homo_sapiens/regulatory_features/RegulatoryFeatureActivity/"
    echo "Exiting."
    exit 1
fi

#Download all cell types:
#GFF is 1-based
for cell in ${cells}; do
    echo "Downloading cell type : $cell"
    checksum=$(wget -O- -q ${ensftp}/pub/release-${Ensembl_release}/regulation/homo_sapiens/RegulatoryFeatureActivity/${cell}/CHECKSUM | cut -f1 -d' ')
    if (( "$reuse" > 0 )) && [[ ! -z "$(find $targetDir -name ${cell}.gff.gz | head -1)" ]] && [[ "$noSums" == "1" || "$(md5sum $(find $targetDir -name ${cell}.gff.gz | head -1) | cut -d' ' -f1)" == "$checksum" ]]; then
        info "File found and has the right checksum. Skipping download..."
        fn=$(find $targetDir -name ${cell}.gff.gz | head -1)
        if [[ ! "$fn" -ef ${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz ]];then
            mv "$fn" ${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz
        fi
    else
        axel \
          -q \
          ${ensftp}/pub/release-${Ensembl_release}/regulation/homo_sapiens/RegulatoryFeatureActivity/${cell}/homo_sapiens.*Regulatory_Build.regulatory_activity.*.gff.gz \
          -o ${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz
        testFile "${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz"
        checkGZfile "${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz"
        if [[ "$noSums" == "0" && "$(md5sum ${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz | cut -d' ' -f1)" != "$checksum" ]]; then
            echo "[Error] Checksum invalid ($(md5sum ${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz | cut -d' ' -f1)). The download probably failed. Please rerun with the reuse option (-r) to retry."
            rm ${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz
            exit 1
        fi
    fi
done
echo "Done"

# Printing out report of the downloaded cell types:
cellTypeCount=$(ls -la ${targetDir}/${today}/EnsemblRegulation/*gff.gz | wc -l)
info "Number of downloaded cell types: ${cellTypeCount}\n\n"

#=================================== APPRIS ==================================================

#Downloading APPRIS data:
mkdir -p ${targetDir}/${today}/APPRIS
appris_checksum=$(wget -q -O- http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/md5checksums.txt | grep appris_data.principal.txt | cut -d' ' -f1)
if (( "$reuse" > 0 )) \
  && [[ ! -z "$(find $targetDir -name appris_data.principal.txt | head -1)" ]] \
  && [[ "$(md5sum $(find $targetDir  -name appris_data.principal.txt | head -1) | cut -d' ' -f1)" == "$appris_checksum" ]]; then
    info "Appris file found and has the right checksum. Skipping download..."
    fn=$(find $targetDir -name appris_data.principal.txt | head -1)
    if [[ ! "$fn" -ef ${targetDir}/${today}/APPRIS/appris_data.principal.txt ]]; then
        mv "$fn" ${targetDir}/${today}/APPRIS/appris_data.principal.txt
    fi
else
    info "Downloading APPRIS isoform data\n"
    info "Download from the current release folder. Build: GRCh38, for GENCODE version: ${GENCODE_release}\n"
    axel -a \
      http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt \
      -o ${targetDir}/${today}/APPRIS/appris_data.principal.txt
    
    downloaded_appris_checksum=$(md5sum ${targetDir}/${today}/APPRIS/appris_data.principal.txt | cut -d' ' -f1)
    if [[ "$noSums" == "0" && "$appris_checksum" != "$downloaded_appris_checksum" ]]; then
        echo "[Error] Downloaded checksum ($downloaded_appris_checksum) differs from expected ($appris_checksum). Download probably failed. Rerun with reuse (-r) option."
        exit 1
    fi
    # Testing if the file exists or not:
    testFile "${targetDir}/${today}/APPRIS/appris_data.principal.txt"
    info "Download complete\n\n"
fi

#=================================== moved UP to do downloads first
# Downloading scores
if [[ $getScores == "yes" ]]; then
    info "Processing scores.\n"
    if (( "$reuse" > 0 )) \
      && [[ -s "$outdir/scores/eigen.phred_v2.dat" ]] \
      && [[ "$noSums" == "1" || $(md5sum $outdir/scores/eigen.phred_v2.dat | cut -d' ' -f1) == "2346005c1cd457bb5ff48c64667736b2" ]]; then
        info "Eigen scores file found and has the right checksum. Skipping download..."
        cat <(grep -v "EigenPath=" ${configfile}) <(echo "EigenPath=${outdir}/scores/eigen.phred_v2.dat") | sponge ${configfile}
    else
        info "Downloading Eigen Phred scores\n"
        mkdir -p $outdir/scores
        cd $outdir/scores

        axel -a ftp://anonymous@ftpexchange.helmholtz-muenchen.de:21021/ticketnr_3523523523525/eigen.phred_v2.dat
        axel -a ftp://anonymous@ftpexchange.helmholtz-muenchen.de:21021/ticketnr_3523523523525/eigen.phred_v2.dat.tbi
    
        if [[ $? -ne 0 ]]; then
            echo "Error: could not download Eigen scores (ftp://anonymous@ftpexchange.helmholtz-muenchen.de:21021/ticketnr_3523523523525/eigen.phred_v2.dat)\n"
            echo "Try downloading later\n"
        fi
        localcksm=$(md5sum $outdir/scores/eigen.phred_v2.dat | cut -d' ' -f1)
        if [[ "$noSums" == "0" && "$localcksm" != "2346005c1cd457bb5ff48c64667736b2" ]]; then
            echo "[Error] Downloaded checksum ($localcksm) differs from expected (2346005c1cd457bb5ff48c64667736b2). Download probably failed. Rerun with reuse (-r) option."
            rm  $outdir/scores/eigen.phred_v2.dat
            exit 1
        else
            cat <(grep -v "EigenPath=" ${configfile}) <(echo "EigenPath=${outdir}/scores/eigen.phred_v2.dat") | sponge ${configfile}
        fi
        cd ..
    fi
fi

if [[ $getCadd == "yes" ]]; then
    if (( "$reuse" > 0 )) \
      && [[ -s "$outdir/scores/whole_genome_SNVs.tsv.gz" ]] \
      && [[ "$noSums" == "1" || $(md5sum $outdir/scores/whole_genome_SNVs.tsv.gz | cut -d' ' -f1) == "cb3856be4c3bb969ff8f0a6139ca226f" ]]; then
        info "CADD scores file found and has the right checksum. Skipping download..."
        cat <(grep -v "caddPath=" ${configfile}) <(echo "caddPath=${outdir}/scores/whole_genome_SNVs.tsv.gz") | sponge ${configfile}
    else
        info "Downloading CADD scores\n"
        mkdir -p scores
        cd scores
        axel -a https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
        axel -a https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz.tbi

        if [[ $? -ne 0 ]]; then
            echo "Error: could not download CADD scores (https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz)\n"
            echo "Try downloading later\n"
        fi
        localcksm=$(md5sum $outdir/scores/whole_genome_SNVs.tsv.gz | cut -d' ' -f1)
        if [[ "$noSums" == "0" && "$localcksm" != "cb3856be4c3bb969ff8f0a6139ca226f" ]]; then
            echo "[Error] Downloaded checksum ($localcksm) differs from expected (cb3856be4c3bb969ff8f0a6139ca226f). Download probably failed. Rerun with reuse (-r) option."
            rm  $outdir/scores/whole_genome_SNVs.tsv.gz
            exit 1
        else
            cat <(grep -v "caddPath=" ${configfile}) <(echo "caddPath=${outdir}/scores/whole_genome_SNVs.tsv.gz") | sponge ${configfile}
        fi
    cd ..
    fi
fi

#=============================================================================================
# Stopping early if just downloading

if (( "$justdl" > 0 )); then
  info "Download flag enabled, stopping early."
  exit 0
fi

#=============================================================================================
#Combining APPRIS and GENCODE data

info "Combining APPRIS and GENCODE data.. "
mkdir -p ${targetDir}/${today}/processed
export APPRIS_FILE=${targetDir}/${today}/APPRIS/appris_data.principal.txt

## Warning: we cannot check the contents since they are dynamic

#if (( "$reuse" > 0 )) && [[ ! -z "$(find $targetDir -name Appris_annotation_added.txt.gz | head -1)" ]]; then
#    echo "[Info] Appris has already been combined with Gencode. Skipping merge, but this is unsafe."

zcat ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | grep -v "#" | awk '$3 != "Selenocysteine" && $3 != "start_codon" && $3 != "stop_codon"' \
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
                }' | gzip > ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz

# Test if output is empty:
testFile ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz
testFileLines ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz # 0-based

# OUTPUT:
# start/end are 0-based coordinates of "class"
#
#{"source":"GENCODE","gene_ID":"ENSG00000186092","appris":"NA","start":"65419","chr":"1","strand":"+","class":"gene","end":"71585"}
#{"transcript_ID":"ENST00000641515","end":"71585","class":"transcript","strand":"+","appris":"Minor","start":"65419","chr":"1","source":"GENCODE","gene_ID":"ENSG00000186092"}
#{"strand":"+","start":"65419","chr":"1","appris":"Minor","gene_ID":"ENSG00000186092","exon_ID":"ENSE00003812156","source":"GENCODE","transcript_ID":"ENST00000641515","end":"65433","class":"exon"}

echo "Done"

# Print out report:
appris_lines=$(zcat ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz | wc -l | awk '{print $1}')
info "Number of Appris annotated GENCODE annotations: ${appris_lines}\n\n"

#=============================================================================================

##
## Step 7. Pre-processing cell specific regulatory data
##
info "Aggregating cell specific information of regulatory features... "
#CellTypes=$( ls -la ${targetDir}/${today}/EnsemblRegulation/ | perl -lane 'print $1 if  $F[-1] =~ /RegulatoryFeatures_(.+).gff.gz/ ' )
CellTypes=$( ls -la ${targetDir}/${today}/EnsemblRegulation/ | perl -lane 'print $1 if  $F[-1] =~ /(.+).gff.gz/ ' )
for cell in ${CellTypes}; do
    export cell
    fn=${targetDir}/${today}/EnsemblRegulation/${cell}.gff.gz
    checkGZfile ${fn}

    # parsing cell specific files (At this point we only consider active features. Although repressed regions might also be informative):

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
' | sort -k1,1 -k2,2n -T ${tempdir} | bgzip -f > ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz # 0-based coordinates here


# OUTPUT:
#
#1    9800    10400   ENSR00000344264 chr=1;start=9800;end=10400;class=CTCF_binding_site;regulatory_ID=ENSR00000344264;Tissues=HUVEC|HeLa_S3|mammary_epithelial_1
#1    13400   13600   ENSR00000344265 chr=1;start=13400;end=13600;class=CTCF_binding_site;regulatory_ID=ENSR00000344265;Tissues=HUVEC|NHLF|keratinocyte

# Test if output is empty or not:
testFileLines ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz

tabix -f -p bed ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz
echo  "Done."

# Print out report:
cellSpecFeatLines=$(zcat ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz | wc -l | awk '{print $1}')
info "Number of cell specific regulatory features: $cellSpecFeatLines\n\n"


#=============================================================================================

##
## Step 8. Combine individual files from GTEx tar.gz file into one BED file
##

# no PAR_Y IDs in GTEx (v.8) files, so we don't check

tmpGTEx=${targetDir}/${today}/processed/GTEx_tmp.bed
info "Creating GTEx bed file ... "
listOfGTExFiles=$(tar -tf ${GTExFile} | grep "signif_variant")
for f in ${listOfGTExFiles};do
    g=$(basename ${f})
    tissue=$(echo ${g}|perl -lne '$x="NA";if (/^([^.]+)\./){$x=$1;} print $x;')
    export tissue

    # taking care of deletions as well
    tar -xf ${GTExFile} ${f} -O | zcat - | tail -n +2 | perl -F"\t" -lane '($chr, $pos, $ref, $alt, $build) = split("_", $F[0]);($gene) = $F[1] =~ /(ENSG\d+)\./;$tissue=$ENV{tissue};$,="\t";$chr=~s/^chr//;$start=$pos-1;$end=$pos;if (length($ref)>length($alt)){$end=$start+length($ref)-1;}  print $tissue,$chr,$start,$end,$F[0],$gene;'
done > ${tmpGTEx}

cat ${tmpGTEx} | perl -F"\t" -lane '$tissue=$F[0];$chr=$F[1];$start=$F[2];$end=$F[3];$ID=$F[4];$gene=$F[5];$H{$ID}{chr}=$chr;$H{$ID}{start}=$start;$H{$ID}{end}=$end;push( @{$H{$ID}{genes}{$gene}}, $tissue ); END {foreach $id (keys %H){$chr=$H{$id}{chr};$start=$H{$id}{start};$end=$H{$id}{end};foreach $gene (keys %{$H{$id}{genes}}){$tissues = join "|", @{$H{$id}{genes}{$gene}};print "$chr\t$start\t$end\tgene=$gene;rsID=$id;tissue=$tissues";}}}' | sort -k1,1 -k2,2n -T ${tempdir} > ${targetDir}/${today}/processed/GTEx.bed # 0-based

echo "Done"
#rm -f ${tmpGTEx}

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

intersectBed -wb -a ${targetDir}/${today}/processed/GTEx.bed -b ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz 2>/dev/null | perl -MData::Dumper -MJSON -F"\t" -lane '
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
    ' | gzip > ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz

# OUTPUT
# we don't know which of the rsIDs corresponds to which GTEx tissue
#{"chr":"chr2","GTEx_rsIDs":[...],"gene_ID":"ENSG00000235584","GTEx_tissues":["Thyroid","Lung"],"class":"open_chromatin_region","Tissues":["MM_1S"],"source":"GTEx","start":"96131175","regulatory_ID":"ENSR00000613314","end":"96131863"}

# Testing if output file has lines:
testFileLines ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz # start/end are 0-based
echo "Done"

# Generate report:
GTExLinkedFeatures=$( zcat ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz | wc -l | awk '{print $1}')
info "Number of GTEx linked regulatory features: ${GTExLinkedFeatures}\n\n"

#=============================================================================================


##
## Step 10. Using intersectbed. Find overlap between genes and regulatory regions
##
info "Linking genes to regulatory features based on overlap ... "
# generating a file.
# 0-based
zcat ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | awk '$3 == "gene"' | perl -lane '
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
    ' | sort -k1,1 -k2,2n -T ${tempdir} | bgzip -f > ${targetDir}/${today}/processed/genes.bed.gz # 0-based

# OUTPUT:
#
#1    11868   14409   ID:ENSG00000223972;Name:DDX11L1


# Intersect bed output:
# 1    16048    29570    ID:ENSG00000227232;Name:WASH7P    1    16048    30847    ENSR00000528774    chr=1;start=16048;end=30847;class=CTCF_binding_site;regulatory_ID=ENSR00000528774;Tissues=DND-41|HMEC|HSMMtube|IMR90|K562|MultiCell|NHDF-AD
intersectBed -wb -a ${targetDir}/${today}/processed/genes.bed.gz -b ${targetDir}/${today}/processed/Cell_spec_regulatory_features.bed.gz -sorted 2>/dev/null | perl -MData::Dumper -MJSON -F"\t" -lane '
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
    ' | bgzip -f > ${targetDir}/${today}/processed/overlapping_features.txt.gz # 0-based coordinates
echo "Done"

# OUTPUT:
#
#{"start":"13400","chr":"1","gene_name":"DDX11L1","gene_ID":"ENSG00000223972","regulatory_ID":"ENSR00000344265","class":"CTCF_binding_site","source":"overlap","Tissues":["HUVEC","NHLF","keratinocyte"],"end":"13600"}

# Generate report:
OverlapLinkedFeatures=$( zcat ${targetDir}/${today}/processed/overlapping_features.txt.gz | wc -l | awk '{print $1}')
info "Number of regulatory features linked by overlap: ${OverlapLinkedFeatures}\n\n"


#=============================================================================================

##
## Step 11. Merging all the components together create sorted, compressed bedfile.
##
info "Merging GENCODE, GTEx and overlap data ..."
export gene_file=${targetDir}/${today}/processed/genes.bed.gz

#GENES
#
#1    11869   14409   ID:ENSG00000223972.5;Name:DDX11L1

#INPUT:
#
#{"start":"13400","chr":"chr1","gene_name":"DDX11L1","gene_ID":"ENSG00000223972","regulatory_ID":"ENSR00000344265","class":"CTCF_binding_site","source":"overlap","Tissues":["HUVEC","NHLF","keratinocyte"],"end":"13600"}
#{"chr":"chr2","GTEx_rsIDs":[...],"gene_ID":"ENSG00000235584","GTEx_tissues":["Thyroid","Lung"],"class":"open_chromatin_region","Tissues":["MM_1S"],"source":"GTEx","start":"96131175","regulatory_ID":"ENSR00000613314","end":"96131863"}
#{"appris":"NA","start":"11869","chr":"chr1","strand":"+","source":"GENCODE","gene_ID":"ENSG00000223972","end":"14409","class":"gene"}


zcat ${targetDir}/${today}/processed/overlapping_features.txt.gz \
     ${targetDir}/${today}/processed/GTEx_Regulation_linked.txt.gz \
     ${targetDir}/${today}/processed/Appris_annotation_added.txt.gz  \
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
        }'  2> ${targetDir}/${today}/failed | sort -k1,1 -k2,2n -T ${tempdir} > ${targetDir}/${today}/Linked_features.bed # 0-based

echo -e "Done\n"

# source == GENCODE => chr,start,end in the 5th field are those of transcript,gene,exon,UTR,CDS
# source == GTEx => chr,start,end in the 5th field are those of regulatory element
# source == overlap => chr,start,end in the 5th field are those of regulatory element

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

info "Output file was saved as: Linked_features.bed.gz\n"
totalLines=$(zcat ${targetDir}/${today}/Linked_features.bed.gz | wc -l | awk '{print $1}')
info "Total number of lines in the final files: ${totalLines}\n"

# FOR LATER USE
mv -f ${targetDir}/${today}/Linked_features.bed.gz ${outdir}
mv -f ${targetDir}/${today}/Linked_features.bed.gz.tbi ${outdir}

#=============================================================================================

# GENCODE basic annotation, IDs are not modified
# coordinates are 1-based

zcat  ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | grep -v "^#"| perl -F"\t" -lane 'next if $F[2] ne "gene";$x=$F[8];$id="NA";$id=$1 if ($x=~/gene_id\s+\"(ENSG[^"]+)\"/); $gn="NA"; $gn=$1 if $x=~/gene_name\s+\"([^"]+)\"/;$,="\t";$F[0]=~s/^chr//;print $F[0],$F[3],$F[4],$gn,$id;' | gzip > ${outdir}/gencode.basic.annotation.tsv.gz
zcat  ${targetDir}/${today}/GENCODE/gencode.v${GENCODE_release}.annotation.gtf.gz | grep -v "^#"| perl -F"\t" -lane 'next if $F[2] ne "gene";$x=$F[8];next unless $x=~/gene_type\s+\"protein_coding\"/;$id="NA";$id=$1 if ($x=~/gene_id\s+\"(ENSG[^"]+)\"/); $gn="NA"; $gn=$1 if $x=~/gene_name\s+\"([^"]+)\"/;$,="\t";$F[0]=~s/^chr//;print $F[0],$F[3],$F[4],$gn,$id;' | gzip > ${outdir}/gencode.basic.annotation.protein_coding.tsv.gz

#==================================== OUTPUT config.txt ======================================

cat <(grep -v "Linked_features=" ${configfile}) <(echo "Linked_features=${outdir}/Linked_features.bed.gz") | sponge ${configfile}
cat <(grep -v "gencode_file=" ${configfile}) <(echo "gencode_file=${outdir}/gencode.basic.annotation.tsv.gz") | sponge ${configfile}

#=============================================================================================

# Report failed associations:
FailedAssoc=$(wc -l ${targetDir}/${today}/failed | awk '{print $1}')
FailedGenes=$( cat ${targetDir}/${today}/failed | perl -lane '$_ =~ /(ENSG\d+)/; print $1' | sort -T ${tempdir} | uniq | wc -l )
FailedSources=$( cat ${targetDir}/${today}/failed | perl -lane '$_ =~ /"source":"(.+?)"/; print $1' | sort -T ${tempdir} | uniq | tr "\n" ", " | sed 's/,$/\n/')
info "Number of lost associations: ${FailedAssoc}, belonging to ${FailedGenes} genes in the following sources: ${FailedSources}\n\n"

if [[ $backup == "yes" ]]; then
    info "Backup enabled, will not delete $targetDir."
    # tar czf ${targetDir}/${today}/${today}_annotation.backup.tar.gz --remove-file ${targetDir}/${today}/APPRIS  \
    # ${targetDir}/${today}/EnsemblRegulation ${targetDir}/${today}/failed ${targetDir}/${today}/GENCODE  \
    # ${targetDir}/${today}/processed
    # info "Intermediate files are saved in ${today}_annotation.backup.tar.gz\n"
else
    rm -rf ${targetDir}
fi


info "The End\n"

echo "Ended: $(date +'%Y-%m-%d %H:%M:%S')" >> $metafile