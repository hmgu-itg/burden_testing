#!/usr/bin/env Rscript

inform=function(...){
  cat(paste("INFO :", ... , "\n"))
}

err=function(...){
  cat(paste("ERROR :", ... , "\n"))
}


validateSMMATOutput = function(group.file, meta.files.prefix){
  dname=dirname(meta.files.prefix)
  fname=basename(meta.files.prefix)
  n.cohort=length(meta.files.prefix)
  n.files=sapply(1:length(meta.files.prefix), function(i){length(list.files(dname[i], pattern=paste0(fname[i], ".var.*")))})
  inform("Provided", n.cohort, "cohorts.")
  for(i in 1:n.cohort){
    inform("Fileset", fname[i], "\t\t", n.files[i], "files")
  }

  library(data.table)
  group.info=tryCatch(fread(group.file, header = FALSE,
    col.names = c("group", "chr", "pos", "ref", "alt", "weight"),
    colClasses = c("character","character","integer","character","character","numeric")),
  error=function(e){err("Failed to read file", group.file,". The error message was:");print(e);stop()})
  n.groups=length(unique(group.info$group))
  inform("Read", n.groups, "groups from file", group.file)
  group.info[,c("vid", "group.idx"):=list(paste(group, chr, pos, ref, alt, sep=":"), as.numeric(factor(group.info$group)))]
  group.info[,c("group", "chr", "ref", "alt", "pos"):=NULL]
  scores=list()
  numvar=list()
  for(i in 1:n.cohort) {

    scores[[i]]=tryCatch(
      rbindlist(lapply(1:n.files[i], function(j) {
        fn=paste0(meta.files.prefix[i], ".score.", j)
        x=fread(fn);
        x[,file:=j];
        inform("Read", nrow(x), "variants in", length(unique(x$group)),"groups from file", fn);
        flush.console()
        return(x)})),
    error=function(e){
      err("Failed to read",n.files[i], "files with prefix", meta.files.prefix[i],". FYI, files matched were:\n", paste(list.files(dname[i], pattern=paste0(fname[i], ".score.*")), sep="\n", collapse="\n"));
      err("\n\nThe error message was:");print(e);stop()})

    scores[[i]][,vid:=paste(group, chr, pos, ref,alt, sep=":")];
    scores[[i]]=merge(scores[[i]], group.info, by="vid");
    scores[[i]][,vid:=NULL];
    scores[[i]]=unique(scores[[i]])
    if(nrow((scores[[i]][,length(unique(file)), by="group"])[V1>1])>1){
      err("Some groups have their scores distributed across two files:")
      err("Prefix\t\t", meta.files.prefix[i])
      print(scores[[i]][,length(unique(file)), by="group"][V1>1])
      stop()
    }
    numvar[[i]]=scores[[i]][,list(.N, unique(group.idx), unique(file)),by=group]
    setnames(numvar[[i]], c("group", "numvar", "group.idx", "file"))
  }

  # first opening
  cons=lapply(1:n.cohort, function(j) file(paste0(meta.files.prefix[j], ".var.1"), "rb"))
  file.idx=rep(1, n.cohort)
  cat("Validating groups :\n")
  pb <- txtProgressBar(min = 0, max = n.groups, style = 3)
  grp.read=rep(0, n.cohort)

  for(i in 1:n.groups){
    if(i %% ceiling(n.groups/100) == 0) {
        setTxtProgressBar(pb, i)
    }
    for(j in 1:n.cohort){
      if(i %in% numvar[[j]]$group.idx){
        ## test whether we need to change file
        if(numvar[[j]][group.idx==i]$file != file.idx[j]){
          cat("\n")
          inform("Closing file", paste0(meta.files.prefix[j], ".var.", file.idx[j]), "after reading", i-1, "groups.")
          close(cons[[j]])
          file.idx[j]=file.idx[j]+1
          inform("Opening file", paste0(meta.files.prefix[j], ".var.", file.idx[j]))
          cons[[j]]=file(paste0(meta.files.prefix[j], ".var.", file.idx[j]), "rb")
        }
        tmp.n.p=numvar[[j]][group.idx == i]$numvar
        tmp.V <- matrix(0, tmp.n.p, tmp.n.p)
        lotri=tryCatch(
          {val = readBin(cons[[j]], what = "numeric", n = (1+tmp.n.p)*tmp.n.p/2, size = 4)
          tmp.V[lower.tri(tmp.V, diag = TRUE)] = val}
          ,
          error=function(e){
          err("\n\nError: failed to read group ", i, " (",numvar[[j]][group.idx==i]$group,").\n
          Score file ", paste0(meta.files.prefix[j], ".var.", file.idx[j]), " may be truncated.\n")
          stop(e)},
        warning=function(w){
          err("\nSucceeded in reading group but vector is the wrong size.")
          err("\nExpected ", (1+tmp.n.p)*tmp.n.p/2, "numbers but received", length(val));
          err("group ", i, " (",numvar[[j]][group.idx==i]$group,").\n
          Score file ", paste0(meta.files.prefix[j], ".var.", file.idx[j]))
          stop()})
        grp.read[j]=grp.read[j]+1

      }
    }

  }
  setTxtProgressBar(pb, n.groups)
  sapply(cons, close)
  cat("\n")
  inform("Finished checking", i ,"groups, no error encountered. All done.")
}
