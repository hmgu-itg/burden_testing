#!/bin/bash

version="v1.0 Last modified: 2020.Nov.30"
today=$(date "+%Y.%b.%d")

# Folder with the variant selector script:
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
regionSelector=${scriptDir}/"selectVariantsSMMAT.py"

configFile=""
chunksTotal=1
chunkNo=""

# SMMAT Python selector wrapper

# --- print out help message and exit:
function display_help() {
    echo "$1"
    echo ""
    echo "Create group file"
    echo "version: ${version}"
    echo ""
    echo "Usage: $0 <parameters>"
    echo ""
    echo "Variant filter options:"
    echo "     -C  - config file (reguired, no default)"
    echo "     -i  - input variant list, tab-separated, bgzipped and TABIX indexed (required, no default, each line has to have 5 fields: chr,pos,ID,ref,alt)"
    echo "     -g  - comma separated list of GENCODE features (gene, exon, transcript, CDS or UTR)"
    echo "     -e  - comma separated list of GTEx features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)"
    echo "     -l  - comma separated list of overlap features (promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind or allreg)"
    echo "     -x  - extend genomic regions by this amount (bp) (default: 0)"
    echo "     -o  - include variants with severe consequences only (more severe than missense)"
    echo ""
    echo "Parameters to set up scores for variants:"
    echo "     -s  - apply weighting; accepted values: CADD, EigenPhred"
    echo ""
    echo "Gene list and chunking:"
    echo "     -L  - file with gene IDs (if not specified all genes will be analyzed)."
    echo "     -d  - total number of chunks (default: 1)."
    echo "     -c  - chunk number (default: 1). Takes precedence over SLURM_ARRAY_TASK_ID etc."
    echo ""
    echo "General options:"
    echo "     -w  - output directory (required, no default)"
    echo ""
    echo "Other options:"
    echo "     -h  - print this message and exit"
#    echo "     -v  - increase verbosity from from INFO to DEBUG"
    echo ""
    echo ""

    exit 1
}

# --- Capture command line options --------------------------------------------

if [ $# == 0 ]; then display_help; fi

OPTIND=1
score=""
geneListFile=""
verbosity="debug"

#while getopts ":hL:c:d:bg:m:s:l:e:x:k:t:ofw:jC:V:" optname; do
while getopts ":hL:c:d:bg:s:l:e:x:ow:C:i:v" optname; do
    case "$optname" in
      # Gene list related parameters:
        "L") geneListFile=${OPTARG} ;;
        "c") chunkNo=${OPTARG} ;;
        "d") chunksTotal=${OPTARG} ;;

      # variant filter parameters:
        "g") gencode=${OPTARG} ;;
        "s") score=${OPTARG} ;;
        "l") overlap=${OPTARG} ;;
        "e") gtex=${OPTARG} ;;
        "x") xtend=${OPTARG} ;;
        "o") lof=1 ;;
        "C") configFile=${OPTARG} ;;
        "i") inputFile=${OPTARG} ;;
#        "v") verbosity="debug" ;;

      # Other parameters:
        "w") outputDir=${OPTARG} ;;
        "h") display_help ;;
        "?") display_help "[Error] Unknown option $OPTARG" ;;
        ":") display_help "[Error] No argument value for option $OPTARG" ;;
        *) display_help "[Error] Unknown error while processing options" ;;
    esac
done

if [[ -z "${inputFile}" ]]; then
    "[Error] input file not specified!"
    exit 1
fi

# -------------------------------- TESTING INPUT FILE FORMAT -------------------------------------

echo "Checking compression integrity"
if ! gzip -t ${inputFile};then
    echo "[Error]: GZIP compression integrity failed for $inputFile"
    exit 1
fi

echo "Checking format"
zcat $inputFile | awk 'BEGIN{FS="\t";x=0;}{if (NF!=5){print "Wrong number of fields (should be 5) in line" NR;print $0;print "";x=1;} if ($2 !~ /^[1-9][0-9]*$/){print "Coordinate (" $2 ") in line " NR " should be an integer";print "";x=1;} if ($4 !~ /^[ACGT]+$/){print "REF allele (" $4 ")in line " NR " has wrong format";print "";x=1;} if ($5 !~ /^[ACGT]+$/){print "ALT allele (" $5 ")in line " NR " has wrong format";print "";x=1;}}END{exit x;}'

if [[ $? != 0 ]];then
    echo "[Error]: $inputFile has wrong format"
    exit 1
fi

# ------------------------------------------------------------------------------------------------


if [[ -z "${outputDir}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Output directory not specified"
    exit 1
fi

# full path
outputDir=`readlink -f $outputDir`
outputDir=${outputDir%/}

if [[ -z "${configFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Config file was not specified"
    exit 1
fi

if [[ ! -e "${configFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Config file does not exist: $configFile"
    exit 1
fi

no_list_warning=""
if [[ -z ${geneListFile} ]];then
    gencode_file=$(grep "^gencode_file" ${configFile} | cut -f 2 -d '=')
    if [[ ! -e ${gencode_file} ]];then
	echo "[Error] GENCODE file (${gencode_file}) specified in the config file (${configFile} does not exist)"
	exit 1
    fi

    totalGenes=$(zcat ${gencode_file} | cut -f 4 | wc -l)
    no_list_warning="[Warning] No gene list specified; using all genes from $gencode_file"
else
    totalGenes=$(cat ${geneListFile} | wc -l)    
fi

if [ $totalGenes -lt $chunksTotal ];then
    echo "[Error] Number of chunks ($chunksTotal) is larger than the number of genes in the gene list ($totalGenes) "
    exit 1
fi

commandOptions=" --config ${configFile} --smmat ${inputFile} "

# -----------------------------------------------------------------------------------------------------------------------------

# setting chunkNo

chunk_warning=""
if [[ ! -z ${SLURM_ARRAY_TASK_ID} ]];then
    if [[ ! -z ${chunkNo} ]];then
	chunk_warning="WARNING: both SLURM_ARRAY_TASK_ID and chunkNo ( -c ) are defined; using chunkNo"
    else
	chunkNo=${SLURM_ARRAY_TASK_ID}
    fi
elif [[ ! -z ${LSB_JOBINDEX} ]];then
    if [[ ! -z ${chunkNo} ]];then
	chunk_warning="WARNING: both LSB_JOBINDEX and chunkNo ( -c ) are defined; using chunkNo"
    else
	chunkNo=${LSB_JOBINDEX}
    fi
else
    if [[ -z ${chunkNo} ]];then
	chunkNo=1 # default
    fi
fi

if [[ ${chunkNo} -gt ${chunksTotal} ]];then
    echo "ERROR: current chunk number ($chunkNo) is greater than the total number of chunks ($chunksTotal). EXIT"
    exit 1
fi

# -----------------------------------------------------------------------------------------------------------------------------

outprefix="group_file"
# GENCODE features
if [[ ! -z "${gencode}" ]]; then
    commandOptions="${commandOptions} --gencode ${gencode}"
    str=$( echo "${gencode}" | perl -lane '$_ =~ s/^\.//;$_ =~ s/,/_/g; print $_;')
    outprefix=${outprefix}"_GENCODE_"${str}
fi

# GTEx - expecting a comma separated list of regulatory features.
if [[ ! -z "$gtex" ]]; then
    commandOptions="${commandOptions} --gtex ${gtex}"
    str=$( echo "${gtex}" | perl -lane '$_ =~ s/^\.//;$_ =~ s/,/_/g; print $_;')
    outprefix=${outprefix}"_GTEx_"${str}
fi

# Overlap - expecting a comma separated list of regulatory features.
if [[ ! -z "${overlap}" ]]; then
    commandOptions="${commandOptions} --overlap ${overlap}"
    str=$( echo "${overlap}" | perl -lane '$_ =~ s/^\.//;$_ =~ s/,/_/g; print $_;')
    outprefix=${outprefix}"_overlap_"${str}
fi

# If lof is set, only variants with severe consequences will be selected.
if [[ ! -z "$lof" ]]; then
    commandOptions="${commandOptions} --lof "
    outprefix=${outprefix}"_severe"
fi

warning1=""
warning2=""
score_tmp=""
# Score - If score is not given we apply no score. Otherwise we test the submitted value:
# Accepted scores:
if [[ ! -z "${score}" ]]; then
    score="${score^^}"
    case "${score}" in
        EIGENPHRED )   score="EigenPhred";;
        CADD )         score="CADD";;
        * )            score_tmp="noweight";;
    esac
else
    score="noweight"
fi

if [[ ! -z ${score_tmp} ]];then
    warning1="[Warning] Submitted score name ($score) is not recognized! Accepted scores: CADD, EigenPhred."
    warning2="[Warning] No scoring will be applied."
    score="noweight"
fi

# Only adding score to command line if score is requested:
if [[ "${score}" != "noweight" ]]; then
    commandOptions="${commandOptions} --score ${score}";
fi
outprefix=${outprefix}"_score_"${score}

if [[ ! -z "${xtend}" ]]; then
    commandOptions="${commandOptions} --extend ${xtend}"
    outprefix=${outprefix}"_extend_"${xtend}
fi

outputDir2=${outputDir}"/"${outprefix}
mkdir -p ${outputDir2}
if [[ ! -d ${outputDir2} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Could not create ${outputDir2}"
    exit 1
fi

# -----------------------------------------------------------------------------------------------------------------------------

outFile="group_file_gene_set."${chunkNo}
outputDir3=${outputDir2}"/chunk_${chunkNo}"
mkdir -p ${outputDir3}
if [[ ! -d ${outputDir3} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Could not create ${outputDir3}"
    exit 1
fi

commandOptions="${commandOptions} --output-dir ${outputDir3}"
LOGFILE=${outputDir3}/"make_group_file_gene_set.${chunkNo}.log"

#--------------------------------------------------------------------------------------------
if [[ -z ${geneListFile} ]];then
    gencode_file=$(grep "^gencode_file" ${configFile} | cut -f 2 -d '=')
    geneListFile="${outputDir3}/gencode_gene_list.txt"
    zcat ${gencode_file} | cut -f 4 > "${geneListFile}"
fi

if [[ ! -e "${geneListFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene list file could not be opened: $geneListFile" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] MAKE GROUP FILE DONE" >> ${LOGFILE}
    exit 1
fi

# -----------------------------------------------------------------------------------------------------------------------------
# creating gene list

rem=$(( totalGenes % chunksTotal ))
chunkSize=$(( totalGenes / chunksTotal ))
lastChunkSize=$(( chunkSize + rem ))
inputGeneList=${outputDir3}"/input_gene.list"

if [[ ${chunkNo} -eq ${chunksTotal} ]];then
    tail -n ${lastChunkSize} ${geneListFile} > ${inputGeneList}
else
    awk -v cn="${chunkNo}" -v cs="${chunkSize}" 'NR > (cn-1)*cs && NR <= cn*cs' ${geneListFile} > ${inputGeneList}
fi

n=$( cat ${inputGeneList} | wc -l)
if [[ $n -eq 0 ]];then
    echo "[Error]: Chunk ${chunkNo} is empty" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] MAKE GROUP FILE DONE" >> ${LOGFILE}
    exit 1
fi

# -----------------------------------------------------------------------------------------------------------------------------

if [[ ! -z ${warning1} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` ${warning1} >> ${LOGFILE}
fi

if [[ ! -z ${warning2} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` ${warning2} >> ${LOGFILE}
fi

if [[ ! -z ${chunk_warning} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` ${chunk_warning} >> ${LOGFILE}
fi

if [[ ! -z ${no_list_warning} ]];then
    echo `date "+%Y.%b.%d_%H:%M"` ${no_list_warning} >> ${LOGFILE}
fi

# --- Reporting parameters ------------------------------------------------------

echo `date "+%Y.%b.%d_%H:%M"` "##"  >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "## Version ${version}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "## Date: ${today}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "##" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] General options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Variant selector: ${regionSelector}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Output directory: ${outputDir}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Gene list options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Gene list file: ${geneListFile}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Number of chunks the gene list is split into: ${chunksTotal}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Current chunk: ${chunkNo}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Number of genes in one chunk: ${chunkSize}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Variant filtering options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "input file: ${inputFile}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "GENCODE features: ${gencode:--}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "GTEx features: ${gtex:--}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Overlapping reg.features: ${overlap:-NA}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Features are extended by ${xtend:-0}bp" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] Weighting options:" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"`  "Weighting: ${score}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"`  "[Info] Verbosity: ${verbosity}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

echo `date "+%Y.%b.%d_%H:%M"` "[Info] command line options for the selector script: ${commandOptions}" >> ${LOGFILE}
echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

# -----------------------------------------------------------------------------------------------------------------------------

selectorLog=${outputDir3}/"selector_chunk_"${chunkNo}.log
echo `date "+%Y.%b.%d_%H:%M"` "Calling ${regionSelector} --input ${inputGeneList} ${commandOptions} --verbose ${verbosity} --log ${selectorLog}"  >> ${LOGFILE}
${regionSelector} --input ${inputGeneList} ${commandOptions} --verbose ${verbosity} --log ${selectorLog}
echo `date "+%Y.%b.%d_%H:%M"` "[Info] MAKE GROUP FILE DONE" >> ${LOGFILE}
