#! /bin/bash

set -euo pipefail

TAG=$1 #sampleID which will be used as directory name
WDIR=`pwd` #setting current working path to a variable
EPOCHS=150 #epochs value
EXPECTED_CELLS=5000 #expected cells value
DROPLETS=150000 #total droplets value
FPR=0.01 #fpr value
LEARN=0.0001 #learning value

# add singularity to the path
PATH="/software/singularity-v3.6.4/bin:$PATH"

# set the path to the image we want to use
CELLBENDER_IMAGE=/nfs/cellgeni/singularity/images/cellbender0.2.0-pytorch1.11-cuda11.3.1-commit2507742.sif

# path to output folder (samples will have a folder inside this one)
OUTPUT_FOLDER="$WDIR/results/$TAG/outs"

# create output folder if it does not exist
if [[ ! -d "${OUTPUT_FOLDER}" ]]; then
    mkdir -p "${OUTPUT_FOLDER}"
fi

# setting input path
INPUT_FILE="$WDIR/cellbender_data/$TAG/raw_gene_bc_matrices_h5.h5"

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

singularity run --nv --bind /nfs,/lustre $CELLBENDER_IMAGE cellbender remove-background \
    --input "$INPUT_FILE" \
    --output "${OUTPUT_FOLDER}/cellbender_out.h5" \
    --cuda \
    --expected-cells $EXPECTED_CELLS \
    --epochs $EPOCHS \
    --total-droplets-included $DROPLETS \
    --fpr $FPR \
    --learning-rate $LEARN
