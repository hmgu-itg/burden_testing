#!/usr/bin/bash

if [[ "$#" -lt 4 ]];then
    echo "ERROR: $# arguments are provided; need at least 4"
    exit 1
fi

# group file
gfile=$1
# pheno name
pname=$2

if [[ ! -f "$gfile" ]];then
    echo "ERROR: group file ($gfile) is not a regular file"
    exit 1
fi

declare -a array
array+=($gfile)
declare -a tmpdirs

argc=$#
argv=("$@")

# the rest of the arguments are input dirs
for (( i=2; i<argc; i++ )); do
    if [[ ! -d "${argv[i]}" ]];then
	echo "ERROR: ${argv[i]} is not a directory"
	exit 1
    fi
    tmp_dir=$(mktemp -d -t run.meta.pheno-XXXXXXX)
    tmpdirs+=($tmp_dir)
    echo "INFO: created temp dir $tmp_dir"
    cp -v "${argv[i]}"/*."$pname".* "$tmp_dir"
    array+=($tmp_dir)
    echo
done

echo "${array[@]}"
./run.meta "${array[@]}"
echo

echo "INFO: deleting temp dirs:"
for d in "${tmpdirs[@]}";do
    echo "INFO: $d"
    rm -rf "$d"
done
echo

exit 0
