#! /bin/bash

set -euo pipefail

im=/nfs/cellgeni/singularity/images/cellbender0.2.0-pytorch1.11-cuda11.3.1-commit2507742.sif

dir=$1 #sampleID which will be used as directory name
epochs=150 #epochs value
cells=5000 #expected cells value
droplets=15000 #total droplets value
fpr=0.01 #fpr value
learn=0.0001 #learning value

cd ../data/$dir

h5=raw_gene_bc_matrices_h5.h5
outdir=cellbender

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

[[ ! -e "$h5" ]] && echo "No $h5 in $dir" && false

echo $dir ok
mkdir -p $outdir

singularity run --nv --bind /nfs,/lustre $im cellbender remove-background \
    --input $h5 \
    --output "${outdir}/cellbender_out.h5" \
    --cuda \
    --expected-cells $cells \
    --epochs $epochs \
    --total-droplets-included $droplets \
    --fpr $fpr \
    --learning-rate $learn
