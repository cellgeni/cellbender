#! /bin/bash
#BSUB -G cellgeni
#BSUB -J sr.visium
#BSUB -o %J.cb2mtx.out
#BSUB -e %J.cb2mtx.err
#BSUB -q normal
#BSUB -n 2
#BSUB -M24000
#BSUB -R "span[hosts=1] select[mem>24000] rusage[mem=24000]"
#BSUB -q normal

im=/nfs/cellgeni/singularity/images/r4.1.0-seurat4.3.0-full.sif
singularity run --bind /nfs,/lustre $im Rscript actions/cellbender/scripts/cellbender_output_to_mtx.R
