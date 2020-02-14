#!/bin/bash

# phenofile should have 2 tab separated columns, no header
function checkPhenoFile {
    fname=$1
    code=$(cat $fname|awk 'BEGIN{FS="\t";c=0;}NF!=2{c=1;}END{print c;}')
    return $code
}

version="v12 Last modified: 2020.Feb.14"
today=$(date "+%Y.%b.%d")

# Folder with the variant selector script:
scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MONSTER=$(which MONSTER)
imputation_method='-A' # The default imputation method is BLUP, slowest, but the most accurate. For other options, see MONSTER documentation.

# --- print out help message and exit:
function display_help() {
    echo "$1"
    echo ""
    echo "Genome-wide Monster wrapper"
    echo "version: ${version}"
    echo ""
    echo "This script was written to run MONSTER genome wide"
    echo ""
    echo "Usage: $0 <parameters>"
    echo ""
    echo "Options:"
    echo "     -i  - input directory (required, no default)"
    echo "     -c  - chunk number (if not specified, all chunks will be analyzed)"
    echo "     -p  - phenotype name (required, no default)"
    echo "     -P  - phenotype file (two tab separated columns, no header; required, no default)"
    echo "     -K  - kinship matrix (required, no default)"
    echo ""
    echo "Other options:"
    echo "     -h  - print this message and exit"
    echo ""
    echo ""

    exit 0
}

# --- Capture command line options --------------------------------------------

if [ $# == 0 ]; then display_help; fi

OPTIND=1
inputDir=""
chunkNo=""
targetDirs=()
while getopts ":hp:P:K:i:c:" optname; do
    case "$optname" in
      # MONSTER input files:
        "p" ) phenotype=${OPTARG} ;;
        "P" ) phenotypeFile=${OPTARG} ;;
        "K" ) kinshipFile=${OPTARG} ;;
        "i" ) inputDir=${OPTARG} ;;
        "c" ) chunkNo=${OPTARG} ;;

      # Other parameters:
        "h") display_help ;;
        "?") display_help "[Error] Unknown option $OPTARG" ;;
        ":") display_help "[Error] No argument value for option $OPTARG" ;;
        *) display_help "[Error] Unknown error while processing options" ;;
    esac
done

if [[ -z "${inputDir}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Input directory not specified"
    exit 1
fi

if [[ ! -e "${inputDir}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Input directory does not exist"
    exit 1
fi

if [[ ! -z ${chunkNo} ]];then
    if [[ ${chunkNo} -le 0 ]];then
	echo `date "+%Y.%b.%d_%H:%M"` "[Error] Chunk number should be a positive integer; provided value: ${chunkNo}"
	exit 1
    fi
fi

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
	chunkNo=0 # default, all chunks
    fi
fi

if [[ ${chunkNo} -eq 0 ]];then
    for f in $(find ${inputDir} -maxdepth 1 -type d -name "gene_set*" | sort);do
	targetDirs+=( $f )
    done
else
    targetDirs+=${inputDir}"/gene_set.${chunkNo}"
fi

if [[ -z "${phenotypeFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Phenotype file has to be specified"
    exit 1
elif [[ ! -e "${phenotypeFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Phenotype file could not be opened: $phenotypeFile"
    exit 1
fi

if [[ -z "${kinshipFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Kinship file has to be specified!"
    exit 1
elif [[ ! -e "${kinshipFile}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Kiship file could not be opened: $kinshipFile"
    exit;
fi

# Checking if phenotype is provided and if that phenotype file:
if [[ -z "${phenotype}" ]]; then
    echo `date "+%Y.%b.%d_%H:%M"` "[Error] Phenotype was not set"
    exit 1
fi

for targetDir in ${targetDirs[@]}; do
    if [[ ! -e ${targetDir} ]];then
	echo "[Warning]: ${targetDir}; skipping"
	continue
    fi
    
    LOGFILE=${targetDir}/"MONSTER.log"

    dname=$(basename ${targetDir})
    cn=$(echo ${dname} | sed 's/.*\([0-9][0-9]*\)$/\1/')
    selectorLog=${targetDir}"/chunk_${cn}.output.log"
    if [[ ! -e ${selectorLog} ]];then
	echo `date "+%Y.%b.%d_%H:%M"` "[Warning] Selector log (${selectorLog}) does not exist; skipping" >> ${LOGFILE}
	continue
    fi

    # ------------------------------- Reporting parameters --------------------------
    echo `date "+%Y.%b.%d_%H:%M"` "##"  >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "## Genome-wide Monster wrapper version ${version}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "## Date: ${today}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "##" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Options:" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"`  "MONSTER executable: ${MONSTER}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"`  "Input directory: ${inputDir}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"`  "Working directory: ${targetDir}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"`  "Phenotype: ${phenotype}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"`  "Imputation method: ${imputation_method:-BLUP}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"`  "Kinship matrix: ${kinshipFile}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"`  "Phenotype file: ${phenotypeFile}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"`  "Phenotype: ${phenotype}" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` "" >> ${LOGFILE}

    cd ${targetDir}

    # We have to check if both files are generated AND they have enough lines.
    gene_notenough=$(cat ${selectorLog} | grep -c NOT_ENOUGH_VAR)
    gene_toomany=$(cat ${selectorLog} | grep -c TOO_MANY_VAR)
    gene_noremain=$(cat ${selectorLog} | grep -c NO_VAR_REMAIN)
    gene_absent=$(cat ${selectorLog} | grep -c NO_GENE)
    region_absent=$(cat ${selectorLog} | grep -c NO_REGION)

    echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] ERROR REPORTING FROM VARIANT SELECTOR" >> ${LOGFILE}
    echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] =====================================" >> ${LOGFILE}

    if [[ "$gene_notenough" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Not enough variants [NOT_ENOUGH_VAR]: $(cat ${selectorLog} | grep NOT_ENOUGH_VAR | sed 's/.*Gene.//;s/ .*//' | tr '\n' ' ')" >> ${LOGFILE}
    fi
    if [[ "$gene_toomany" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Too many variants [TOO_MANY_VAR]: $(cat ${selectorLog} | grep TOO_MANY_VAR | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
    fi
    if [[ "$gene_noremain" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] All scoring failed [NO_VAR_REMAIN]: $(cat ${selectorLog} | grep NO_VAR_REMAIN | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
    fi
    if [[ "$gene_absent" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] Gene name unknown [NO_GENE]: $(cat ${selectorLog} | grep NO_GENE | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
    fi
    if [[ "$region_absent" -ne 0 ]]; then
        echo `date "+%Y.%b.%d_%H:%M"` -e "[Warning] No region in gene [NO_REGION]: $(cat ${selectorLog} | grep NO_REGION | sed 's/.*Gene.//;s/ .*//'| tr '\n' ' ')" >> ${LOGFILE}
    fi

    if [[ ! -e gene_set_output_genotype_file.txt ]]; then
	echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set ${targetDir} has failed. No genotype file has been generated; skipping" >> ${LOGFILE}
	continue
    elif [[ $(cat gene_set_output_genotype_file.txt | wc -l ) -lt 2 ]]; then
	echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set in ${targetDir} has failed, genotype file is empty; skipping" >> ${LOGFILE}
	continue
    elif [[ ! -e gene_set_output_variant_file.txt ]]; then
	echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set in ${targetDir} has failed, SNP file was not generated; skipping" >> ${LOGFILE}
	continue
    elif [[ $( cat gene_set_output_variant_file.txt | wc -l ) -lt 1 ]]; then
	echo `date "+%Y.%b.%d_%H:%M"` "[Error] Gene set in ${targetDir} has failed, SNP file is empty; skipping" >> ${LOGFILE}
	continue
    fi

    # ASSUMING PHENOTYPE FILE HAS 2 COLUMNS: ID PHENOTYPE, TAB DELIMITED, NO HEADER
    checkPhenoFile ${phenotypeFile}
    if [[ $? -eq 1 ]];then
	echo `date "+%Y.%b.%d_%H:%M"` "[Error] wrong format of the phenotype file ($phenotypeFile); skipping"
	continue
    fi

    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Extracting samples with present phenotype values from the phenotype file and saving them in 01.pheno.txt." >> ${LOGFILE}
    cat ${phenotypeFile} | awk 'BEGIN{FS="\t";OFS="\t";}{if ($2!="NA") {printf "1\t%s\t0\t0\t0\t%s\n", $1, $2} }' > 01.pheno.txt # subset of samples without NA phenotype values
    n1=$(cat ${phenotypeFile} | wc -l)
    n2=$(cat 01.pheno.txt | wc -l)
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Original number of samples in the phenofile: ${n1}; samples in 01.pheno.txt: ${n2}" >> ${LOGFILE}
    echo  >> ${LOGFILE}

    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Selecting sample names from gene_set_output_genotype_file.txt that are already in 01.pheno.txt; saving the result in 02.pheno.ordered.txt" >> ${LOGFILE}
    head -n 1 gene_set_output_genotype_file.txt | tr "\t" "\n" | tail -n +2 |sort| perl -lne 'BEGIN {open $pf, "< 01.pheno.txt";while ($l = <$pf>){chomp $l;@a = split(/\s/, $l);$h{$a[1]} = $l;}close($pf);}{print $h{$_} if exists $h{$_}}' > 02.pheno.ordered.txt
    n3=$(cat 02.pheno.ordered.txt | wc -l)
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Samples in 02.pheno.ordered.txt: ${n3}" >> ${LOGFILE}
    echo  >> ${LOGFILE}

    # Get samples from gene_set_output_genotype_file.txt that are not in the 02.pheno.ordered.txt file:
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Getting samples from gene_set_output_genotype_file.txt that are not in the 02.pheno.ordered.txt file; saving results in 03.samples.to.exclude.txt" >> ${LOGFILE}
    cut -f 2 02.pheno.ordered.txt > temp.txt
    head -n 1 gene_set_output_genotype_file.txt | cut -f 2- | tr "\t" "\n" >> temp.txt
    sort temp.txt | uniq -u > 03.samples.to.exclude.txt
    rm temp.txt
    n4=$(cat  03.samples.to.exclude.txt | wc -l)
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Samples in 03.samples.to.exclude.txt: ${n4}" >> ${LOGFILE}
    echo  >> ${LOGFILE}

    # From the genotype file, extract only those samples that are present in the pheno file:
    echo `date "+%Y.%b.%d_%H:%M"` "[info] Removing samples in 03.samples.to.exclude.txt from gene_set_output_genotype_file.txt; saving results in 04.genotype.filtered.txt" >> ${LOGFILE}
    echo  >> ${LOGFILE}
    head -n 1 gene_set_output_genotype_file.txt | tr "\t" "\n" | perl -lne 'BEGIN {open $pf, "< 03.samples.to.exclude.txt";while ($l = <$pf>){chomp $l;$h{$l} = 1;}close($pf);}{push @a, $. unless exists $h{$_}} END{$s = sprintf("cut -f%s gene_set_output_genotype_file.txt > 04.genotype.filtered.txt", join(",", @a));`$s`}'

    # Generate a mapping file that helps to convert IDs to numbers:
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Generate SED sample ID to integer mapping file; saving result in 05.sample.map.sed" >> ${LOGFILE}
    echo  >> ${LOGFILE}
    cut -f 2 02.pheno.ordered.txt | awk '{printf "s/%s/%s/g\n", $1, NR+2 }' > 05.sample.map.sed

    # Generate an inclusion list with the samples to be kept:
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Generate an inclusion list from 02.pheno.ordered.txt with the samples to be kept; saving result in 06.samples.to.keep.txt" >> ${LOGFILE}
    echo  >> ${LOGFILE}
    cut -f 2 02.pheno.ordered.txt > 06.samples.to.keep.txt

    # Get the kinship matrix:
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Selecting samples in 06.samples.to.keep.txt from the kinship file;saving result in 07.kinship.filtered.txt" >> ${LOGFILE}
    echo  >> ${LOGFILE}
    R --slave -e 'library(data.table); mlong=fread("'$kinshipFile'"); tokeep=fread("06.samples.to.keep.txt", header=F)$V1; direct=mlong[(mlong$V2 %in% tokeep) & (mlong$V3 %in% tokeep),]; mapping = fread("05.sample.map.sed", sep="/", header=FALSE);direct$V3 = mapping[match(direct$V3, mapping$V2),]$V3; direct$V2 = mapping[match(direct$V2, mapping$V2),]$V3;write.table(direct, file="07.kinship.filtered.txt", quote=FALSE, sep=" ", col.names = FALSE, row.names=FALSE)'

    # Remap IDs and remove special characters from the snp, phenotype and genotype files:
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Changing sample names in 02.pheno.ordered.txt and 04.genotype.filtered.txt (using 05.sample.map.sed); saving results in 08.pheno.ordered.txt and 09.genotype.filtered.txt" >> ${LOGFILE}
    echo  >> ${LOGFILE}
    sed -f 05.sample.map.sed 02.pheno.ordered.txt > 08.pheno.ordered.txt
    sed -f 05.sample.map.sed 04.genotype.filtered.txt > 09.genotype.filtered.txt

    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Renaming variants in 09.genotype.filtered.txt and gene_set_output_variant_file.txt; saving results in 10.genotype.filtered.mod.txt and 11.snpfile.mod.txt" >> ${LOGFILE}
    echo  >> ${LOGFILE}
    cat 09.genotype.filtered.txt | perl -lane '$_ =~ s/[^0-9A-Za-z\-\t\._]//gi;$_ =~ s/_/x/g; print $_'  > 10.genotype.filtered.mod.txt
    cat gene_set_output_variant_file.txt | perl -lane '$_ =~ s/[^0-9A-Za-z\-\t\._]//gi;$_ =~ s/_/x/g; $_ =~ s/Inf/0.0001/g;print $_'  > 11.snpfile.mod.txt

    # Filter out genes which have only monomorphic variants, as it might cause a crash:
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Looking for monomorphic variants in 10.genotype.filtered.mod.txt; saving them in 12.mono.variants.txt" >> ${LOGFILE}
    echo  >> ${LOGFILE}
    tail -n +2 10.genotype.filtered.mod.txt | perl -lne '@f=split(/\s+/);$\="\n";$s=shift(@f);%H=();foreach $x (@f){$H{$x}=1;}if (scalar(keys(%H))==1){print $s;}' > 12.mono.variants.txt

    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Removing monomorphic variants and genes with less than two variants from 11.snpfile.mod.txt; saving result in 13.snpfile.nomono.txt" >> ${LOGFILE}
    echo  >> ${LOGFILE}
    exclude_mono.pl --input 11.snpfile.mod.txt --output 13.snpfile.final.txt --exclude 12.mono.variants.txt 2>mono.genes.txt

    # sorting kinship and genotype files
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Sorting 10.genotype.filtered.mod.txt and 07.kinship.filtered.txt; saving results in 10.genotype.filtered.mod.srt.txt and 07.kinship.filtered.srt.txt" >> ${LOGFILE}
    echo  >> ${LOGFILE}

    sort -k2,2n -k3,3n 07.kinship.filtered.txt > 07.kinship.filtered.srt.txt
    cat 10.genotype.filtered.mod.txt | transpose | sort -k1,1n | transpose > 10.genotype.filtered.mod.srt.txt
    cp 13.snpfile.final.txt 13.snpfile.final.original.txt

    # Calling MONSTER
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] MONSTER call: MONSTER -k 07.kinship.filtered.srt.txt -p 08.pheno.ordered.txt -m 1 -g 10.genotype.filtered.mod.srt.txt  -s 13.snpfile.final.txt ${imputation_method}" >> ${LOGFILE}

    flag=0
    while true; do
	MONSTER  -k 07.kinship.filtered.srt.txt -p 08.pheno.ordered.txt -m 1 -g 10.genotype.filtered.mod.srt.txt  -s 13.snpfile.final.txt ${imputation_method} 2>>${LOGFILE}

	# Break the loop if the run was successful
	if [[ $? -eq 0 ]]; then break; fi

	# No MONSTER.out
	if [[ ! -e MONSTER.out ]]; then
	    echo `date "+%Y.%b.%d_%H:%M"` "[Error] MONSTER failed before creating the output file" >> ${LOGFILE}
	    break
	fi

	# Empty MONSTER.out
	if [[ $( cat MONSTER.out | wc -l) -eq 0 ]]; then
            echo `date "+%Y.%b.%d_%H:%M"` "[Error] MONSTER.out is empty. Trying to run MONSTER gene by gene" >> ${LOGFILE}
	    flag=1
            break
	fi
	
	# Test if we've analyzed all genes
	if [[ $(awk 'BEGIN{FS="\t";}NF==5{print $0;}' MONSTER.out | tail -n +2 | wc -l ) == $(cat 13.snpfile.final.txt | wc -l) ]]; then break; fi

	# Only header is in MONSTER.out
	if [[ $( cat MONSTER.out | wc -l) -eq 1 ]]; then
            firstGene=$(cut -f 1 13.snpfile.final.txt | head -n 1)
            echo `date "+%Y.%b.%d_%H:%M"` "[Warning] It seems that the first gene (${firstGene}) has failed. Re-running MONSTER after excluding it." >> ${LOGFILE}
            awk 'BEGIN{FS="\t";}NR>1{print $0;}' 13.snpfile.final.txt | sponge 13.snpfile.final.txt
	else
            lastGene=$(awk 'BEGIN{FS="\t";}NF==5{print $1;}' MONSTER.out | tail -n 1 )
	    failedGene=$(awk -v g=${lastGene} 'BEGIN{FS="\t";f=0;}{if (f==1 && $1!=g){print $1;exit;} if ($1==g){f=1;}}' 13.snpfile.final.txt)
            echo `date "+%Y.%b.%d_%H:%M"` "[Warning] Monster has failed after ${lastGene}, next gene (${failedGene}) is removed and re-run." >> ${LOGFILE}
	    
            awk -v g=${failedGene} 'BEGIN{FS="\t";}$1!=g{print $0;}' 13.snpfile.final.txt | sponge 13.snpfile.final.txt
	fi
    done

    # gene by gene
    if [[ $flag -eq 1 ]];then
	k=$(cat 13.snpfile.final.txt | wc -l)
	for i in $(seq 1 $k);do
	    head -n $i 13.snpfile.final.txt | tail -n 1 > temp.snpfile.txt
	    gene=$(head -n $i 13.snpfile.final.txt | tail -n 1 | cut -f 1)
	    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Trying only one gene ${gene}" >> ${LOGFILE}
	    MONSTER  -k 07.kinship.filtered.srt.txt -p 08.pheno.ordered.txt -m 1 -g 10.genotype.filtered.mod.srt.txt -s temp.snpfile.txt ${imputation_method} 2>>${LOGFILE}
	    if [[ ! -e MONSTER.out ]]; then
		echo `date "+%Y.%b.%d_%H:%M"` "[Error] MONSTER failed before creating the output file for gene ${gene}" >> ${LOGFILE}
	    elif [[ $( cat MONSTER.out | wc -l) -eq 0 ]]; then
		echo `date "+%Y.%b.%d_%H:%M"` "[Error] MONSTER.out for gene ${gene} is empty" >> ${LOGFILE}
	    else
		mv MONSTER.out MONSTER.out."$i"	    
	    fi
	done
	
	echo -e 'SNP_set_ID\tn_individual\tn_SNP\trho_MONSTER\tp_MONSTER' > MONSTER.out

	shopt -s nullglob
	mfiles=(MONSTER.out.*)
	if [[ ${#mfiles[@]} == 0 ]];then
	    echo `date "+%Y.%b.%d_%H:%M"` "[Error] MONSTER.out files were not found"  >> ${LOGFILE}
	else
	    for f in "${mfiles[@]}"
	    do
		cat $f | grep -v "SNP_set_ID" >> MONSTER.out
	    done
	    rm MONSTER.out.* temp.snpfile.txt
	fi
	shopt -u nullglob	
    fi
    
    # Copying MONSTER.out to the root directory:
    if [[ -e MONSTER.out ]]; then
	cp MONSTER.out ../MONSTER.${phenotype}.${cn}.out
    else
	echo `date "+%Y.%b.%d_%H:%M"` "[Error] MONSTER.out file was not found"  >> ${LOGFILE}
    fi

    cp ${selectorLog} ..

    # Compress folder:
    echo `date "+%Y.%b.%d_%H:%M"` "[Info] Compressing and removing files" >> ${LOGFILE}
    tar -zcvf gene_set.${cn}.tar.gz *
    mv gene_set.${cn}.tar.gz ..
    cd .. && rm -rf ${targetDir}    
done







