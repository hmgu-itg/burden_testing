#!/usr/bin/env Rscript
list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
for(lib in list.of.packages){
  library(lib, character.only=T)
}
argv=commandArgs(T)
## expected arguments: OUTFILE, followed by list of files following Step1 naming convention
OUTFILE=argv[1]

argv=argv[-1]
cat(paste("\t[INFO] Outputting to", OUTFILE, "\n"))

MAFfromAF=function(x){
    out=x
    out[x>0.5]=1-out[x>0.5]
    return(out)
}

merged=NULL
for(f in argv){
    cohort=fread(f)
    setnames(cohort, c("chr", "pos", "ref", "alt", "ac", "an"))
    cohort[,c("id", "af") := list(paste(chr, pos, ref, alt), ac/an)]
    cohort[,maf := MAFfromAF(af)]
    cohort_name=sub("\\..*", "", f)
    cat(paste("\t[INFO] Read cohort", cohort_name, "\n"))
    stopifnot(!(is.null(cohort_name) | cohort_name==""))
    cols=colnames(cohort)
    cols[cols!="id"]=paste0(cohort_name, ".", cols[cols!="id"])
    setnames(cohort, cols)
    if(is.null(merged)){
      merged=cohort
    }else{
      merged=merge(merged, cohort, all=T, by="id")
      cat(paste("\t[INFO] Successfully merged cohort", cohort_name, ", df is now of dim", paste(dim(merged), collapse=" x "), "\n"))
    }
}

#print(colnames(merged))
mafs=grep("\\.maf", colnames(merged))
mafs=colnames(merged)[mafs]
ans=grep("\\.an", colnames(merged))
ans=colnames(merged)[ans]
merged[,"avgmaf" := weighted.mean(unlist(mget(mafs)), w=unlist(mget(ans)), na.rm=T), by = seq_len(nrow(merged))]
merged=merged[avgmaf<=0.05,]
cat(paste("\t[INFO] Filtering out avgMAF leaves",nrow(merged),"rows.\n"))
flush.console()
maxall=merged[,lapply(.SD, max, na.rm=T), .SDcols=ans]
maxall=unlist(maxall)
merged[,an:=sum(maxall[!is.na(unlist(.SD))]),by=seq_len(nrow(merged)),.SDcols=ans]
merged[,miss := 1-rowSums(.SD, na.rm=T)/(an), .SDcols = ans]
merged=merged[miss<0.01,]
cat(paste("\t[INFO] Filtering out missingness leaves",nrow(merged),"rows.\n"))
merged=data.table(id=merged[,id])
merged[, c("chr", "pos", "ref", "alt") := tstrsplit(id, " ", fixed=TRUE)]
#add_column(merged, id="-", .after = "pos")
merged[,id := "-"]
#print(head(merged))
#merged=as.data.table(merged)
merged[,chr := sub("chr", "", chr)]
setcolorder(merged, c("chr", "pos", "id", "ref", "alt"))
merged[, c("chr", "pos") := lapply(.SD, as.numeric), .SDcols=c("chr", "pos")]
setorder(merged, chr, pos)
#head(merged)
fwrite(merged, OUTFILE, quote=F, col.names=F, sep="\t")
cat(paste("\t[INFO] File",OUTFILE,"has been written. Bgzipping and tabixing.\n"))
system(paste("bgzip", OUTFILE))
system(paste0("tabix -s 1 -b 2 -e 2 ", OUTFILE, ".gz"))
cat("Done.\n")
