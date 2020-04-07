#!/usr/bin/env Rscript
library(GMMAT)
library(data.table)
## Expects 2 args: a directory and a variant set file
argv = commandArgs(trailingOnly = TRUE)
dir=argv[1]
group_file=argv[2]

singlecohort_files=list.files(dir, pattern="*.singlecohort.out")
cohorts=gsub(".singlecohort.out", "", singlecohort_files, fixed=T)
cat(paste("[INFO] detected", length(cohorts), "cohorts :", paste(cohorts, collapse=" "), "\n"))
nfiles=sapply(cohorts, function(x){
    scorefiles=list.files(dir, pattern=paste0(unlist(x),".score.*"))
    return(length(scorefiles))
})

cat("[INFO] Meta-analysing...\n")
flush.console()
meta=SMMAT.meta(meta.files.prefix=paste0(dir, "/", cohorts), n.files=nfiles,
           group.file=group_file,
        MAF.range=c(1e-10, 0.05),
        MAF.weights.beta = c(1,1),
        miss.cutoff=0.01,
        tests=c("O", "E"),
        rho=(0:10)/10,
        use.minor.allele=T
          )

cat("[INFO] Done. Writing results.\n")
fwrite(meta, paste(c(cohorts, "meta", "txt"), collapse="."), sep="\t", na=NA, quote=F)
