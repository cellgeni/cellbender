#! /bin/bash

set -euo pipefail

#run this from work directory

script=../actions/cellbender-0-2-0.sh #selecting cellbender script to run
sf=../actions/samples.txt #selecting sample file

cpu=8 #selecting cpus
mem=40000 #selecting memory 
group="cellgeni" #selecting group to submit with 
que="gpu-cellgeni" #selecting queue to submit to
gpu_mode="mode=shared:j_exclusive=no:gmem=60000:num=1" #selecting gpu conditions: memory, gpu cores, whether to shared gpu
gpu_model="dgx-b11" #selecting gpu model

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

if true; then
  mkdir -p logs
  cat $sf | while read name; do
    bsub -n $cpu -Rspan[hosts=1] -M $mem -a "memlimit=True" -G $group -q $que -gpu $gpu_mode -m $gpu_model -o logs/ooo.$name.%J.txt -e logs/eee.$name.%J.txt $script $name
  done
fi
