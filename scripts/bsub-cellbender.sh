#! /bin/bash

set -euo pipefail

#run this from work directory

script=../actions/cellbender-0-2-0.sh #selecting cellbender script to run
sf=../actions/samples.txt #selecting sample file

CPU=8 #selecting cpus
MEM=40000 #selecting memory 
GROUP="cellgeni" #selecting group to submit with 
QUE="gpu-cellgeni" #selecting queue to submit to
GPU_MODE="mode=shared:j_exclusive=no:gmem=6000:num=1" #selecting gpu conditions: memory, gpu cores, whether to shared gpu
GPU_MODEL="dgx-b11" #selecting gpu model

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

if true; then
  mkdir -p logs
  cat $sf | while read name; do
    bsub -n $CPU -Rspan[hosts=1] -M $MEM -R"select[mem>${MEM}] rusage[mem=${MEM}]" -G $GROUP -q $QUE -gpu $GPU_MODE -m $GPU_MODEL -o logs/ooo.$name.%J.txt -e logs/eee.$name.%J.txt $script $name
  done
fi
