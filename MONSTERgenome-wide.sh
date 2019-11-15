#!/bin/bash

# TODO: add output dir option
# TODO: floor ?

# A wrapper script to automate genome-wide burden testing using MONSTER.
# For more information on the applied method see: http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/

# In this version when the script finds a significant gene, it tests if the association is
## driven by a single variant or more by repeating the test with removing a variant
## each time.

# For more information on how to use the wrapper, please see the README.md file or check the repository:
# https://github.com/wtsi-team144/burden_testing

# Todo for the next generation of the wrapper:
## An older verision of MONSTER had an issue why we had to wrap genes individually for

version="v11 Last modified: 2017.08.01" # This version mainly changes in the documentation.

today=$(date "+%Y.%m.%d") # Get the date

# The variant selector script, that generates snp and genotype input for MONSTER:
regionSelector=Burden_testing.pl

# Folder with all the scripts:
export scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Folder with all the phenotypes: # No longer hardwired. Accept as a command line parameter.
# phenotypeDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/phenotypes
# phenotypeDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/phenotypes/correct_names.andmissing

# Default file with all the gene names:
#geneListFile="${scriptDir}/gene_list.txt"

# Kinship matrix file: # No longer hardwired, accepted as command line parameter.
# kinshipMatrix=/nfs/team144/ds26/burden_testing/kinship/2016.10.20_fix_diagonal/kinship.fixdiag.txt

MONSTER=$(which MONSTER)
missing_cutoff=1 # Above the set missingness threshold, variants will be excluded. Below the missing genotype will be imputed.
imputation_method='-A' # The default imputation method is BLUP, slowest, but the most accurate. For other options, see documentation.
configFile=""

# Testing if MONSTER is in the path, exiting if not:
if [[ -z "${MONSTER}" ]]; then echo "[Error] MONSTER is not in the path. Exiting."; exit; fi

##
## Defining default values:
##
export signifThreshold=1e-5 # By default the hit threshold is 1e-5
export keep_temp="No" # By default we are not keeping temporary files:
export chunkCount=1 # By default we process all genes at one chunk.
export chunkNo=1 # By default we are processing the first chunk.
export MAF=0.05 # By default this is the upper minor allele frequency.
if [[ $LSB_JOBINDEX > 0 ]]; then chunkNo=$LSB_JOBINDEX; fi # If jobindex is available, we set that as chunk number, it might be overridden by -c

# --- print out help message and exit:
function display_help() {
    echo "$1"
    echo ""
    echo "Genome-wide Monster wrapper"
    echo "version: ${version}"
    echo ""
    echo "This script was written to run MONSTER genome wide. This script takes a series of arguments
based on which it calls downstream helper scripts, and generates specific directory for the output files.
This is supposed to be called by bsub, where the chunk number is read from LSB_JOBINDEX.
It pools results together within one chunk."
    echo ""
    echo "Usage: $0 <parameters>"
    echo ""
    echo "Variant filter options:"
    echo "     -g  - list of gencode features."
    echo "     -e  - list of linked GTEx featuress"
    echo "     -l  - list of linked overlapping features."
    echo "     -m  - upper maf thresholds"
    echo "     -x  - extend genomic regions (bp)."
    echo "     -o  - exclude all non loss-of-function variants from the test (less than missense in severity)."
    echo "     -f  - include only HC and LC loftee variants."
    echo "     -j  - include only HC loftee variants."
    echo "     -z  - config file for Burden_testing.pl."
    echo ""
    echo "Parameters to set up scores for variants:"
    echo "     -s  - turn weights on. Arguments: CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred,"
    echo "           Linsight, Mixed"
    echo "     -t  - the value with which the scores will be shifted"
    echo "     -k  - below the specified cutoff value, the variants will be excluded"
    echo ""
    echo "Gene list and chunking:"
    echo "     -L  - list file with the gene IDs (required, no default)."
    echo "     -d  - chunk count: how many chunks the gene list should be split into (default: 1)."
    echo "     -c  - chunk number or jobindex (default: 1)."
    echo ""
    echo "General options:"
    echo "     -w  - working directory (default: current working directory)"
    echo "     -b  - Keep temporary files."
    echo "     -i  - p-value threshold (default: 1e-5) below which the temporary files will be saved."
    echo ""
    echo "Monster parameters:"
    echo "     -p  - phenotype (required, no default)"
    echo "     -P  - phenotype file (required, no default)"
    echo "     -K  - Kinship file (required, no default)"
    echo "     -V  - VCF file (required, no default, Use % character for chromosome name eg 'chr%.vcf.gz')"
    echo ""
    echo "Other options:"
    echo "     -h  - print this help message and exit"
    echo ""
    echo ""

    exit 1
}

# --- If a gene fails at some point, generate a report, and save files in a failed folder.
function failed (){

# TODO: $gene ?
    
    message="${1}"

    echo "[Error] Gene has failed: ${gene}"
    echo "[Error] Gene has failed: ${gene}" >&2 # Reporting to the error log as well.
    echo "[Error] ${message}"
    echo "[Error] folder moved to ${folder}/failed"
    echo ""

    # Adding NA-s to the output file:
    echo -e "${gene}\tNA\tNA" >> ${workingDir}/gene_set.${chunkNo}/results

    # Generating failed directory and compress gene dir:
    mkdir -p ${workingDir}/failed
    tar -czvf ${workingDir}/failed/${gene}.tar.gz -C ${workingDir}/gene_set.${chunkNo}/ ${gene}

    # Once the folder is compressed, we add to the pooled tar archive:
    if [[ ! -e ${workingDir}/failed/failed_genes.tar ]]; then
        tar -c --remove-files -f ${workingDir}/failed/failed_genes.tar -C ${workingDir}/failed  ${gene}.tar.gz # Create archive
    else
        tar -r --remove-files -f ${workingDir}/failed/failed_genes.tar -C ${workingDir}/failed  ${gene}.tar.gz  # Add to archive
        rm ${workingDir}/failed/${gene}.tar.gz
    fi

    # Removing directory:
    cd ${workingDir}
    rm -rf ${workingDir}/gene_set.${chunkNo}/${gene}
    return 0
}

# --- If a p-value is really low, the whole run is backed up:
function savehit (){
    mkdir -p ${workingDir}/hits
    tar -czvf ${workingDir}/hits/${gene}.tar.gz -C ${workingDir}/gene_set.${chunkNo}/ ${gene}
    cd ${workingDir}
    return 0
}

# --- Capture command line options --------------------------------------------

# Help message is printed if no command line arguments has been passed:
if [ $# == 0 ]; then display_help; fi

# Looping through all command line options:
OPTIND=1
while getopts ":hL:c:d:p:P:K:V:bi:g:m:s:l:e:x:k:t:ofw:jC:" optname; do
    case "$optname" in
      # Gene list related parameters:
        "L") geneListFile=${OPTARG} ;;
        "c") export chunkNo=${OPTARG} ;;
        "d") chunkCount=${OPTARG} ;;

      # MONSTER input files:
        "p" ) phenotype=${OPTARG} ;;
        "P" ) phenotypeFile=${OPTARG} ;;
        "K" ) kinshipFile=${OPTARG} ;;
        "V" ) vcfFile="${OPTARG}" ;;

      # Wrapper option:
        "b") keep_temp="yes";;
        "i") signifThreshold=${OPTARG};;

      # variant filter parameters:
        "g") gencode=${OPTARG} ;;
        "m") MAF=${OPTARG} ;;
        "s") score=${OPTARG} ;;
        "l") overlap=${OPTARG} ;;
        "e") gtex=${OPTARG} ;;
        "x") xtend=${OPTARG} ;;
        "k") cutoff=${OPTARG} ;;
        "t") scoreshift=${OPTARG} ;;
        "o") lof=1 ;;
        "f") loftee=1 ;;
        "j") lofteeHC=1 ;;
        "C") configFile=${OPTARG} ;;

      # Other parameters:
        "w") rootDir=${OPTARG} ;;
        "h") display_help ;;
        "?") display_help "[Error] Unknown option $OPTARG" ;;
        ":") display_help "[Error] No argument value for option $OPTARG" ;;
        *) display_help "[Error] Unknown error while processing options" ;;
    esac
done

#--- checking input files - if any of the tests fails, the script exits.---------

if [[ -z "${configFile}" ]]; then
    echo "[Error] Config file was not specified";
    exit;
fi

if [[ ! -e "${configFile}" ]]; then
    echo "[Error] Config file does not exist: $geneListFile";
    exit;
fi

# Phenotype file (Is it set? Does it exists?):
if [[ -z "${phenotypeFile}" ]]; then
    display_help "[Error] Phenotype file has to be specified!";
elif [[ ! -e "${phenotypeFile}" ]]; then
    echo "[Error] Phenotype file could not be opened: $phenotypeFile";
    exit;
fi

# Kinship file (Is it set? Does it exists?):
if [[ -z "${kinshipFile}" ]]; then
    display_help "[Error] Kinship file has to be specified!";
elif [[ ! -e "${kinshipFile}" ]]; then
    echo "[Error] Kiship file could not be opened: $kinshipFile";
    exit;
fi

# TODO
# vcf file (Is it set? Does it exists?):
if [[ -z "${vcfFile}" ]]; then
    display_help "[Error] VCF file has to be specified!";
elif [[ ! -e $( echo "${vcfFile}" | sed -e 's/\%/12/') ]]; then
    echo "[Error] VCF file could not be opened: $vcfFile";
    exit;
else
    commandOptions=" --vcfFile ${vcfFile} "
fi

# Select genome build:
#commandOptions="${commandOptions} --build ${build}"

# Gene list file (Is it set? Does it exists?):
if [[ ! -e "${geneListFile}" ]]; then
    echo "[Error] Gene list file could not be opened: $geneListFile";
    exit;
fi

# Working directory (Is it set? Does it exists?):
if [[ -z "${rootDir}" ]]; then
    rootDir=$(pwd)
elif [[ ! -d "${rootDir}" ]]; then
    echo "[Error] The directory does not exists: $rootDir";
    exit;
fi

# GENCODE -expecting a list of feature names separated by a comma.
if [[ ! -z "${gencode}" ]]; then
    folder="${gencode}"
    commandOptions="${commandOptions} --GENCODE ${gencode}"
fi

# GTEx - expecting a list of feature names separeted by comma.
if [[ ! -z "$gtex" ]]; then
    folder="${folder}.GTEX.${gtex}"
    commandOptions="${commandOptions} --GTEx ${gtex}"
fi

# Overlap - expecting a list of features separeated by comma.
if [[ ! -z "${overlap}" ]]; then
    folder="${folder}.Overlap.${overlap}"
    commandOptions="${commandOptions} --overlap ${overlap}"
fi

# MAF - expecting a float between 0 and 0.5
if [[ ! -z "$MAF" ]]; then
    folder="${folder}.MAF.${MAF}"
    commandOptions="${commandOptions} --maf ${MAF}"
fi

# IF loftee is set, only those variants are included in the test that were predicted to be LC or HC lof variants by loftee:
if [[ ! -z "$loftee" ]]; then
    folder="${folder}.loftee"
    commandOptions="${commandOptions} --loftee "
fi

# IF lofteeHC is set, only those variants are included in the test that were predicted to be HC lof variants by loftee:
if [[ ! -z "$lofteeHC" ]]; then
    folder="${folder}.lofteeHC"
    commandOptions="${commandOptions} --lofteeHC "
fi

# If lof is set, only variants with severe consequences will be selected.
if [[ ! -z "$lof" ]]; then
    folder="${folder}.severe"
    commandOptions="${commandOptions} --lof "
fi

commandOptions="${commandOptions} --configFile ${configFile} "

# Score - If score is not given we apply no score. Otherwise we test the submitted value:
# Accepted scores:
if [[ ! -z "${score}" ]]; then
    score="${score^^}"
    case "${score}" in
        EIGEN )        score="Eigen";;
        EIGENPC )      score="EigenPC";;
        EIGENPHRED )   score="EigenPhred";;
        EIGENPCPHRED ) score="EigenPCPhred";;
        CADD )         score="CADD";;
        LINSIGHT )     score="Linsight";;
        MIXED )        score="Mixed";;
        * )            score="noweight";;
    esac
else
    echo "[Warning] Submitted score name is not recognized! Accepted scores: CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred, Linsight or Mixed."
    echo "[Warning] No scoreing will be applied."
    score="noweight"
fi

# Only adding score to command line if score is requested:
if [[ "${score}" != "noweight" ]]; then
    commandOptions="${commandOptions} --score ${score}";
fi

# If Eigen score is applied, we shift the scores by 1, if no other value is specified:
if [[ ("${score}" == "Eigen") && (-z "${scoreshift}" ) ]]; then scoreshift=1; fi

# Adjust dir name, noweight if no score is specified:
folder="${folder}.${score}"

# Exons might be extended with a given number of bps:
if [[ ! -z "${xtend}" ]]; then
    folder="${folder}.Xtend.${xtend}"
    commandOptions="${commandOptions} --extend ${xtend}"
fi

# Setting score cutoff, below which the variant will be removed from the test:
if [[ ! -z "${cutoff}" ]]; then
    folder="${folder}.cutoff.${cutoff}"
    commandOptions="${commandOptions} --cutoff ${cutoff}"
fi

# TODO: scoreshift, general vs Eigen ?
# Setting score shift, a number that will be added to every scores (MONSTER does not accept scores <= 0!!):
if [[ ! -z "${scoreshift}" ]]; then
    folder="${folder}.shift.${scoreshift}"
    commandOptions="${commandOptions} --shift ${scoreshift}"
fi

# Checking if phenotype is provided and if that phenotype file:
if [[ -z "${phenotype}" ]]; then
    echo "[Error] Phenotype was not set! Exiting."
    exit 1
fi

# Based on the size of the list and the chunk count we determine the size of each chunk:
chunkSize=$(wc -l "${geneListFile}" | perl -lane 'print int($F[0]/$ENV{"chunkCount"}+0.5)')

# Once everything is checked, print out the parameters:


# Updating working dir, and create folder:
folder=$( echo $folder | perl -lane '$_ =~ s/^\.//; print $_;')
workingDir="${rootDir}/${folder}/Pheno.${phenotype}"
mkdir -p "${rootDir}"

# --- Reporting parameters ------------------------------------------------------
echo "##"
echo "## Genome-wide Monster wrapper version ${version}"
echo "## Date: ${today}"
echo "##"
echo ""

echo "[Info] General options:"
echo -e "\tVariant selector: ${regionSelector}"
echo -e "\tScript dir: ${scriptDir}"
echo -e "\tSignificance threshold: ${signifThreshold}"
echo -e "\tWorking directory: ${workingDir}/gene_set.${chunkNo}"
echo -e "\tKeeping temporary files: ${keep_temp}"
echo ""

echo "[Info] Gene list options:"
echo -e "\tGene list file: ${geneListFile}"
echo -e "\tNumber of chunks the gene list is split into: ${chunkCount}"
echo -e "\tCurrent chunk: ${chunkNo}"
echo -e "\tNumber of genes in one chunk: ${chunkSize}"
echo ""

echo "[Info] Variant filtering options:"
echo -e "\tvcf file: ${vcfFile}"
echo -e "\tGENCODE feaures: ${gencode:--}"
echo -e "\tGTEx feaures: ${gtex:--}"
echo -e "\tOverlapping reg.features: ${overlap:-NA}"
echo -e "\tFeatures are extended by ${xtend:-0}bp"
echo -e "\tUpper minor allele frequency: ${MAF:-1}"
echo ""

echo "[Info] Weighting options:"
echo -e "\tWeighting: ${score}"
echo -e "\tScore cutoff: ${cutoff:-0}"
echo -e "\tScores shifted by: ${scoreshift:-0}"
echo -e "\tOutput folder: ${workingDir}"
echo ""

echo "[Info] command line options for burden get region: ${commandOptions}"
echo ""

echo "[Info] MONSTER options:"
echo -e "\tMONSTER executable: ${MONSTER}"
echo -e "\tMissingness: ${missing_cutoff}"
echo -e "\tImputation method: ${imputation_method:-BLUP}"
echo -e "\tKinship matrix: ${kinshipFile}"
echo -e "\tPhenotype file: ${phenotypeFile}"
echo -e "\tPhenotype: ${phenotype}"
echo ""

# --- Main loop executed for all genes --------------------------------------------

# Creating results file with header:
mkdir -p ${workingDir}/gene_set.${chunkNo}
if [[ ! -d ${workingDir}/gene_set.${chunkNo} ]]; then
    echo "[Error] Chunk directory (${workingDir}/gene_set.${chunkNo}) could not be created. Exiting."
fi

# Looping through all genes in the gene set:
awk -v cn="${chunkNo}" -v cs="${chunkSize}" 'NR > (cn-1)*cs && NR <= cn*cs' ${geneListFile} > ${workingDir}/gene_set.${chunkNo}/input_gene.list

# Entering working directory:
cd ${workingDir}/gene_set.${chunkNo};

# Reporting call:
echo "${scriptDir}/${regionSelector}  --build 38 --input input_gene.list --output gene_set_output ${commandOptions} --verbose > output.log"
${scriptDir}/${regionSelector}  --input input_gene.list --output gene_set_output ${commandOptions} --verbose > output.log

# We are expecting to get 2 files: gene_set_output_genotype_file.txt & gene_set_output_SNPinfo_file.txt
echo "[Info] Checking output..."
# We have to check if both files are generated AND they have enough lines.
gene_notenough=$(cat output.log | grep -c NOT_ENOUGH_VAR)
gene_toomany=$(cat output.log | grep -c TOO_MANY_VAR)
gene_noremain=$(cat output.log | grep -c NO_VAR_REMAIN)
gene_absent=$(cat output.log | grep -c NO_GENE)
region_absent=$(cat output.log | grep -c NO_REGION)

echo -e "[Warning] ERROR REPORTING FROM REGION SELECTOR"
echo -e "[Warning] ===================================="

if [[ "$gene_notenough" != 0 ]]; then
        echo -e "[Warning] Not enough variants [NOT_ENOUGH_VAR]:\t $(cat output.log | grep NOT_ENOUGH_VAR | sed 's/.*Gene.//;s/ .*//' | tr '\n' ' ')"
fi
if [[ "$gene_toomany" != 0 ]]; then
        echo -e "[Warning] Too many variants [TOO_MANY_VAR]:\t $(cat output.log | grep TOO_MANY_VAR | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')"
fi
if [[ "$gene_noremain" != 0 ]]; then
        echo -e "[Warning] All scoring failed [NO_VAR_REMAIN]:\t $(cat output.log | grep NO_VAR_REMAIN | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')"
fi
if [[ "$gene_absent" != 0 ]]; then
        echo -e "[Warning] Gene name unknown [NO_GENE]:\t $(cat output.log | grep NO_GENE | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')"
fi
if [[ "$region_absent" != 0 ]]; then
        echo -e "[Warning] No region in gene [NO_REGION]:\t $(cat output.log | grep NO_REGION | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')"
fi

if [[ ! -e gene_set_output_genotype_file.txt ]]; then
    echo "[Error] Gene set ${chunkNo} has failed. No genotype file has been generated. Exiting."
    exit
elif [[ $(cat gene_set_output_genotype_file.txt | wc -l ) -lt 2 ]]; then
    echo "[Error] Gene set ${chunkNo} has failed, genotype file is empty. Exiting."
    exit
elif [[ ! -e gene_set_output_variant_file.txt ]]; then
    echo "[Error] Gene set ${chunkNo} has failed, SNP file was not generated. Exiting."
    exit
elif [[ $( cat gene_set_output_variant_file.txt | wc -l ) -lt 1 ]]; then
    echo "[Error] Gene set ${chunkNo} has failed, SNP file is empty. Exiting."
    exit
fi

# At this point the we have to process the above created files to syncronize wi/nfs/team144/ds26/scripts/burden_testing/MONSTERgenome-wide_updated.sh -g exon -x 50 -e promoter,enhancer,TF_bind -l promoter,enhancer,TF_bind -s EigenPCPhred -c 3 -P /lustre/scratch115/projects/t144_helic_15x/analysis/HA/phenotypes/correct_names.andmissing/MANOLIS.HDL.txt -w /lustre/scratch115/realdata/mdt0/projects/t144_helic_15x/analysis/HA/burdentesting/arthur_rerun -V /lustre/scratch115/realdata/mdt0/projects/t144_helic_15x/analysis/HA/release/postrelease_missingnessfilter/chr%.missingfiltered-0.01_consequences.lof.HWE.vcf.gz -K /lustre/scratch115/realdata/mdt0/projects/t144_helic_15x/analysis/HA/relmat/final/burden/matrix.monster.txt -p HDL -d 50th the phenotype file and the kinship matrix.
# So Monster can process it:

# TODO: remove hardcoded filenames
# TODO: phenotype file format ? pheno file format ?
# Get the phenotype:
echo "[Info] Extracting phenotype."
cat ${phenotypeFile} | grep -v NA | awk 'NR != 1 {printf "1\t%s\t0\t0\t0\t%s\n", $1, $3}' > pheno.txt

# TODO: order means sorted IDs ?
# Order the sample IDs in the phenotype file:
echo "[Info] Re-ordering samples in the phenotype file."
head -n1 gene_set_output_genotype_file.txt | tr "\t" "\n" | tail -n+2 |sort| perl -lane 'BEGIN {open $pf, "< pheno.txt";
            while ($l = <$pf>){chomp $l;@a = split(/\s/, $l);$h{$a[1]} = $l;}}{print $h{$F[0]} if exists $h{$F[0]}}' > pheno.ordered.txt

# Get the list of samples that are missing from the phenotype file:
export samples=$(grep -v -w -f <(cut -f2 pheno.ordered.txt) <(head -n1 gene_set_output_genotype_file.txt | cut -f2- | tr "\t" "\n"))

# From the genotype file, extract only those samples that are present in the phenotype file:
echo "[info] Extracting un-used samples from the genotype file."
head -n1 gene_set_output_genotype_file.txt | tr "\t" "\n" | perl -lane 'BEGIN {foreach $s ( split /\s/, $ENV{"samples"}){$h{$s} = 1;}}{
    push @a, $. unless exists $h{$F[0]}} END{$s = sprintf("cut -f%s gene_set_output_genotype_file.txt > genotype.filtered.txt", join(",", @a));`$s`}'

# Generate a mapping file that helps to convert HELIC IDs to numbers:
echo "[Info] Generate sample mapping file."
cut -f2 pheno.ordered.txt | awk '{printf "s/%s/%s/g\n", $1, NR+2 }' > sample.map.sed

# Generate an inclusion list with the samples to be kept:
cut -f2 pheno.ordered.txt > samples.to.keep.txt

# Get the kinship matrix:
echo "[Info] Processing kinship file."
R --slave -e 'library(data.table); mlong=fread("'$kinshipFile'"); tokeep=fread("samples.to.keep.txt", header=F)$V1; direct=mlong[(mlong$V2 %in% tokeep) & (mlong$V3 %in% tokeep),];  mapping = fread("sample.map.sed", sep="/", header=FALSE);direct$V3 = mapping[match(direct$V3, mapping$V2),]$V3; direct$V2 = mapping[match(direct$V2, mapping$V2),]$V3;write.table(direct, file="filtered_kinship.txt", quote=FALSE, sep=" ", col.names = FALSE, row.names=FALSE)'

# Remap IDs and remove special characters from the snp, phenotype and genotype files:
echo "[Info] Changing IDs and variant names."
sed -i -f sample.map.sed pheno.ordered.txt
sed -i -f sample.map.sed genotype.filtered.txt
# TODO: keep underscores and capital letters too ?
cat genotype.filtered.txt | perl -lane '$_ =~ s/[^0-9a-z\-\t\.]//gi; print $_'  > genotype.filtered.mod.txt
cat gene_set_output_variant_file.txt | perl -lane '$_ =~ s/[^0-9a-z\t\.]//gi; $_ =~ s/Inf/0.0001/g; ;print $_'  > snpfile.mod.txt

# Filter out genes which have only monomorphic variants, as it might cause a crash:
echo "[Info] Looking for monomorphic variants..."
#tail -n+2 genotype.filtered.mod.txt | while read snp genotype ; do if [[ -z $( echo $genotype | awk '$0 ~ 1' ) ]]; then echo $snp; fi; done | awk '{printf "s/%s//g\n", $1}' > mono_remove.sed
tail -n+2 genotype.filtered.mod.txt | perl -lne '@f=split(/\s+/);$\="\n";$s=shift(@f);foreach (@f){$H{$_}=1;}if (scalar(keys(%H))==1){print $s;}' | awk '{printf "s/%s//g\n", $1}' > mono_remove.sed

echo "[Info] Removing monomorphic variants from the SNPs file."
sed -f mono_remove.sed snpfile.mod.txt > snpfile.nomono

echo "[Info] Get genes where only monomorphics remain. Exclude them."
grep -v -w -f <(cat snpfile.nomono | awk 'NF  == 2 {print $1; a+=1}END{if(a == 0){print "noting to remove"}}') snpfile.mod.txt > snpfile.mod.nomono.txt

# Calling MONSTER
echo "[info] MONSTER call: MONSTER -k filtered_kinship.txt -p pheno.ordered.txt -m 1 -g genotype.filtered.mod.txt  -s snpfile.mod.txt ${imputation_method}"
while true; do
    MONSTER -k filtered_kinship.txt -p pheno.ordered.txt -m 1 -g genotype.filtered.mod.txt  -s snpfile.mod.nomono.txt ${imputation_method}

    # We break the loop if the run was successful.
    if [[ $? == 0 ]]; then break; fi
    
    # If we have all the genes we are good:
    if [[ $(awk 'NF == 5' MONSTER.out | cut -f1 | tail -n+2 | wc -l ) == $(cut -f1 snpfile.mod.txt | sort -u | wc -l) ]]; then break; fi

    # Why did it fail?
    if [[ ! -e MONSTER.out ]]; then
        echo "[Error] MONSTER failed before creating the output file. Cannot be resolved. Exiting" >&2;
        break;
    elif [[ $( cat MONSTER.out | wc -l) -eq 1 ]]; then
        firstGene=$(cut -f1 snpfile.mod.nomono.txt | head -n1)
        echo "[Warning] It seems that the first gene (${firstGene}) has failed. Re-run MONSTER." >&2
        grep -vw $(firstGene) snpfile.mod.nomono.txt | sponge snpfile.mod.nomono.txt
    else
        lastGene=$(awk 'NF == 5' MONSTER.out | tail -n1 | cut -f1 )
        failedGene=$(grep -A1 -w ${lastGene} snpfile.mod.nomono.txt | cut -f1 | tail -n1 )
        echo "[Warning] Monster has failed after ${lastGene}, next gene (${failedGene}) is removed and re-run."
        grep -vw ${failedGene} snpfile.mod.nomono.txt | sponge snpfile.mod.nomono.txt
    fi
done

# Once MONSTER is finished, we remove the un-used temporary files:
rm genotype.filtered.txt

# Moving MONSTER.out to the root directory:
if [[ -e MONSTER.out ]]; then
    cp MONSTER.out ../MONSTER.${phenotype}.${chunkNo}.out
else
    echo "[Error] MONSTER.out file was not found. Something went wrong." >&2
fi

# Compress folder:
echo "[Info] Compressing and removing files."
tar -zcvf gene_set.${chunkNo}.tar.gz *
mv gene_set.${chunkNo}.tar.gz ..
cd .. && rm -rf gene_set.${chunkNo}
