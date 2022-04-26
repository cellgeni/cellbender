#!/bin/bash 

set -euo pipefail

sf=../actions/irods.txt
h5=raw_feature_bc_matrix.h5

cat $sf | while read name path; do 
  mkdir -p $name
  ( cd $name
    echo "-- $path"
    iget -f -v -N 4 -K $path/$h5 raw_gene_bc_matrices_h5.h5; 
    echo "✓"
  )
done
