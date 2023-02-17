#! /bin/bash

#SIF=/nfs/cellgeni/singularity/images/cellbender-0.2.0.sif
SIF=/nfs/cellgeni/singularity/images/cellbender0.2.1-pytorch11.1-cuda11.3.1.sif

## this version uses the input dir rather than h5 file 
## it should be produced using multiome.R via "bsub < actions/cellbender/scripts/bsub.extract.gex.from.multiome.sh"
TAG=$1
IN_DIR=../data/$TAG/gex
OUTPUT=$TAG

LOWCOUNTHR=5
NEXP=`cut -d ',' -f4  ../data/$TAG/summary.csv | tail -n1`
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
  --cuda --expected-cells $NEXP $TOTAL \
  --low-count-threshold $LOWCOUNTHR
