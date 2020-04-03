#!/usr/bin/env Rscript
library(data.table)
argv=commandArgs(T)
## expected arguments: list of files following Step1 naming convention

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
    cohort_name=sub("\\..*", "")
    cat(paste("Read cohort", cohort_name))
    stopif(is.null(cohort_name) | cohort_name=="")
    cols=colnames(cohort)
    cols[cols!="id"]=paste0(cohort_name, ".", cols[cols!="id"])
    setnames(cohort, cols)
    if(is.null(merged)){
      merged=cohort
    }else{
      merged=merge(merged, cohort, all=T, by="id")
      cat(paste("Successfully merged cohort", cohort_name, "df is now of dim", paste(nrow(merged), collapse=" x ")))
    }
}
