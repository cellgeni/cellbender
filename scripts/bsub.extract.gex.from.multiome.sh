#! /bin/bash
#BSUB -G cellgeni
#BSUB -J sr.visium
#BSUB -o %J.extract.gex.out
#BSUB -e %J.extract.gex.err
#BSUB -q normal
#BSUB -n 2
#BSUB -M24000
#BSUB -R "span[hosts=1] select[mem>24000] rusage[mem=24000]"
#BSUB -q normal

im=/nfs/cellgeni/singularity/images/r4.1.0-seurat4.3.0-full.sif
singularity run --bind /nfs,/lustre $im Rscript actions/cellbender/scripts/multiome.R
