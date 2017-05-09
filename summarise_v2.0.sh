#!/usr/local/bin/bash


# Setting up environment:
singlePointDir=/lustre/scratch115/projects/t144_helic_15x/analysis/HA/single_point/output/missing_chunks

# Working dir is read from the command line parameters:
workingDir=$(pwd)/$1

# We exit if the directory does not exist:
if [[ ! -d ${workingDir} ]]; then
    echo "[Error] The submitted working dir does not exists. Exiting."
    exit
fi

cd $workingDir

# Initialize output file:
echo -e "Trait\tGene\tP_burden\tSNP_cnt\tSNP\tP_sp\tP_cond" > ${workingDir}/conditioned_pvalues.txt

# Looping through all phenotypes:
for trait in $(ls -d Pheno.* | sed -e 's/Pheno\.//'); do
    
    # Go to the next phenotype if there is no hits:
    if [[ ! -d ${workingDir}/Pheno.${trait}/hits ]]; then continue; fi
    
    # Get single point phenotype (the names of the phenotypes are different):
    sTrait=${trait}
    case $trait in
        'adiponectin' ) sTrait="ADIPONECTIN" ;;
        "BGP" ) sTrait="OSTEOCALCIN" ;;
        "Bilirubin" ) sTrait="BILIRUBIN" ;;
        "Fe_iron" ) sTrait="FEIRON" ;;
        "Ferritin" ) sTrait="FERRITIN" ;;
        "gammaGT" ) sTrait="GAMMAGT" ;;
        "leptin" ) sTrait="LEPTIN" ;;
    esac
    
    # OK enter directory:
    cd ${workingDir}/Pheno.${trait}/hits
    
    # Looping through all the genes:
    for gene in $(ls *tar.gz | sed -e 's/\.tar\.gz//'); do
        
        # Enter directory:
        cd ${workingDir}/Pheno.${trait}/hits/
        
        # Extract p-value:
        read oldpval SNPcnt <<< $(find ${workingDir}/Pheno.${trait} -name results | xargs cat | awk -v gene=${gene} '$1 == gene { print $2, $3 }')

        # decompress file:
        tar zxvf ${gene}.tar.gz
        
        # Does the folder exists?
        if [[ ! -d ${workingDir}/Pheno.${trait}/hits/${gene} ]]; then echo "[Error] $gene has failed!"; continue; fi        
        
        # Enter directory:
        cd ${workingDir}/Pheno.${trait}/hits/${gene}
        SinglePointFile=${singlePointDir}/${sTrait}/MANOLIS.${sTrait}.assoc.missingfilled.txt.bgz

        # Checking if the assoc file is tabixed or not:
        if [[ ! -e ${SinglePointFile} ]]; then
            echo "[Warning] Assoc file for ${trait} is not bgzipped! Re-compressing and indexing!"
            zcat ${SinglePointFile} | bgzip > ${SinglePointFile}
            tabix -s 1 -b 3 -e 3 -S 1 ${SinglePointFile}
        fi
        
        # If index file is not created create:
        if [[ ! -e ${SinglePointFile}.tbi ]]; then
            echo "[Warning] Assoc file for ${trait} was not indexed!"
            tabix -f -s 1 -b 3 -e 3 -S 1 ${SinglePointFile}
        fi
        
        # Get the most significant variant and the corresponding sigle point p-value:
        read highvar Psp <<< $(tabix ${SinglePointFile} $( cat ${gene}_output_variants |  head -n1 | 
            cut -f3- | tr "\t" "\n"  | perl -F"_" -lane '$F[0]=~s/chr//; printf "%s:%s-%s ", $F[0], $F[1], $F[1]' ) | cut -f2,14 | perl -lane 'if ($. == 1){$v = $F[0]; $m = $F[1]; next};
            if ($F[1] < $m){$v = $F[0]; $m = $F[1]; }; END {print "$v $m"}')
        
        # Format variable name:
        varMod=$(echo $highvar | perl -F":" -lane '$F[1]=~s/\[b38\]//g; print $F[0].$F[1]')

        # If the tabix was not successful, we save what we have and go to the next gene:
        if [[ -z  $Psp ]]; then
            echo "[Warning] No single point p-values for $gene in $trait. Gene skipped!"
            echo -e "${trait}\t${gene}\t${oldpval}\t${SNPcnt}\t-\t-\t-" >> ${workingDir}/conditioned_pvalues.txt
            continue;
        fi

        # Adding the genotype as a covariate to the phenotype file:
        #paste <( zgrep -m1 "#CHROM" ${vcfDir}/chr11.missingfiltered-0.01_consequences.vcf.gz | cut -f10- | tr "\t" "\n" )\
        #    <( tabix ${vcfDir}/chr${chr}.missingfiltered-0.01_consequences.vcf.gz $(echo $highvar | perl -F":" -lane '$F[1]=~ s/\[b38\]//; printf "%s:%s-%s", $F[0], $F[1], $F[1]') \
        #        | cut -f10- | tr "\t" "\n" | cut -f1 -d":" | perl -F"/" -lane '$F[0] eq "." ? print "-9" : print $F[0] + $F[1]') \
        #        | sed -f ${scriptDir}/HELIC.to.Num.sed | tr " " "\t" | cut -f1,3 | perl -lane '
        #    $h{$F[0]} = $F[1];
        #    END{
        #        open($ph, "< pheno.mod.ordered.txt");
        #        while ( $l = <$ph> ){
        #            chomp $l;
        #            @a = split("\t", $l);
        #            next unless exists $h{$a[1]};
        #            next if $h{$a[1]} == -9;
        #            print "$l\t$h{$a[1]}"
        #        }
        #    }' > pheno.mod.signifadded.txt        
        
        
        paste pheno.mod.ordered.txt <(grep $varMod genotype.mod.txt | cut -f $(head -n1 genotype.mod.txt | cut -f2- | tr "\t" "\n" \
                | perl -lane 'BEGIN { open $P, "< pheno.mod.ordered.txt"; while ($l = <$P>){ chomp $l; @a = split("\t", $l); $h{$a[1]} = 1}}; push(@b, $. + 1) if exists $h{$F[0]}; END {print join ",", @b} ') | tr "\t" "\n" ) \
                > pheno.mod.signifadded.txt
        
        # Calling monster:
        mv MONSTER.out MONSTER.backup
        /nfs/team144/software/MONSTER_v1.3/MONSTER -k kinship.mod.filtered.txt -p pheno.mod.signifadded.txt -m 1 -g genotype.mod.filtered.txt  -s snpfile.mod.txt
        
        # Process output:
        mv MONSTER.out MONSTER.test.out
        pval=$(tail -n1 MONSTER.test.out | cut -f5)
        
        # Saving data:
        echo -e "${trait:-NA}\t${gene:-NA}\t${oldpval:-NA}\t${SNPcnt:-NA}\t${highvar:-NA}\t${Psp:-NA}\t${pval:-NA}" >> ${workingDir}/conditioned_pvalues.txt
    done
done

cd $workingDir

