#! /bin/bash

set -euo pipefail

#run this from work directory

#script=../actions/cellbender-matrix.sh #selecting cellbender script to run
script=../actions/cellbender-0-2-0.sh # script changes whether using h5 files or matrix format
sf=../actions/sample-list  #selecting sample file

CPU=8 #selecting cpus
MEM=40000 #selecting memory 
GROUP="cellgeni" #selecting group to submit with 
GMEM=6000 #selecting gpu conditions: memory, gpu cores, whether to shared gpu
QUE="gpu-normal" #selecting queue to submit to
# if GMEM < 30GB you can use "gpu-normal", "gpu-huge" or "gpu-basement"
# if you need GMEM >30GB use QUE="gpu-cellgeni-a100"

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

mkdir -p logs
cat $sf | while read name; do
   bsub -n $CPU -Rspan[hosts=1] -M $MEM -R"select[mem>${MEM}] rusage[mem=${MEM}]" -G $GROUP -q $QUE -gpu "mode=shared:j_exclusive=no:gmem=${GMEM}:num=1" -o logs/ooo.$name.%J.txt -e logs/eee.$name.%J.txt $script $name
done
