#/bin/bash

echo "Burden testing pipeline using MONSTER"
echo ""
echo "Analysis consists of two steps"
echo "1) prepare BED file with GENCODE, GTEx and Regulatory Elements features"
echo "   Available options: singularity exec burden.1.1 prepare-regions -h"
echo ""
echo "2) run genome-wide MONSTER burden testing"
echo "   Available options: singularity exec burden.1.1 monster-burden -h"
echo ""



