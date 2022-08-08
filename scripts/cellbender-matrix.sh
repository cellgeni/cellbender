#! /bin/bash

#SIF=/nfs/cellgeni/singularity/images/cellbender-0.2.0.sif
SIF=/nfs/cellgeni/singularity/images/cellbender0.2.1-pytorch11.1-cuda11.3.1.sif

## this version uses the input dir rather than h5 file 
TAG=$1
TYPE=Gene ## change to GeneFull for snRNA-seq!
IN_DIR=$TAG/output/$TYPE/raw/
OUTPUT=$TAG.cellbender.out

echo "WARNING: Processing sample $TAG using input folder $TYPE. If this is single nucleus, please change input directory to GeneFull!"

NEXP=`zcat $TAG/output/$TYPE/filtered/barcodes.tsv.gz | wc -l`
TOTAL=""
if (($NEXP > 20000)) 
then
  NTOT=$((NEXP+10000))
  echo "Modifying presets: expecting more than 20k cells ($NEXP), total number of droplets is $NTOT.."
  TOTAL="--total-droplets-included $NTOT"
else 
  echo "Standard presets: expected number of cells is $NEXP.."
fi

# create output folder if it does not exist
if [[ ! -d $OUTPUT ]]
then
    mkdir -p $OUTPUT
fi

singularity run --nv --bind /nfs,/lustre $SIF cellbender remove-background \
  --input $IN_DIR --output "$OUTPUT/cellbender_out.h5" \
  --cuda --expected-cells $NEXP $TOTAL
