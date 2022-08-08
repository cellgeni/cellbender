#! /bin/bash

set -euo pipefail

#run this from work directory

#script=../actions/cellbender_matrix.sh #selecting cellbender script to run
script=../actions/cellbender-0-2-0.sh # script changes whether using h5 files or matrix format
sf=../actions/sample-list  #selecting sample file

CPU=8 #selecting cpus
MEM=40000 #selecting memory 
GROUP="cellgeni" #selecting group to submit with 
#QUE="gpu-cellgeni-a100" #selecting queue to submit to
QUE="gpu-normal" #if GMEM > 30GB then can use gpu-normal or gpu-huge or gpu-basement
GMEM=6000 #selecting gpu conditions: memory, gpu cores, whether to shared gpu
MODEL="" #selecting gpu model
#MODEL="dgx-b11" #selecting gpu model

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

if true; then
  mkdir -p logs
  cat $sf | while read name; do
    if [ "$MODEL" != "" ]; then
      bsub -n $CPU -Rspan[hosts=1] -M $MEM -R"select[mem>${MEM}] rusage[mem=${MEM}]" -G $GROUP -q $QUE -gpu "mode=shared:j_exclusive=no:gmem=${GMEM}:num=1" -m $MODEL -o logs/ooo.$name.%J.txt -e logs/eee.$name.%J.txt $script $name
    else
      bsub -n $CPU -Rspan[hosts=1] -M $MEM -R"select[mem>${MEM}] rusage[mem=${MEM}]" -G $GROUP -q $QUE -gpu "mode=shared:j_exclusive=no:gmem=${GMEM}:num=1" -o logs/ooo.$name.%J.txt -e logs/eee.$name.%J.txt $script $name
    fi
  done
fi
