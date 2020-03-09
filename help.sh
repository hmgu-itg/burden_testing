#/bin/bash

echo "Burden testing pipeline using MONSTER"
echo ""
echo "Analysis consists of two steps"
echo "1) prepare BED file with GENCODE, GTEx and Regulatory Elements features"
echo "   Available options: singularity exec <container name> prepare-regions -h"
echo ""
echo "2) select variants for MONSTER"
echo "   Available options: singularity exec <container name> select-variants -h"
echo ""
echo "3) run genome-wide MONSTER burden testing"
echo "   Available options: singularity exec <container name> run-monster -h"

