#!/bin/bash 

set -euo pipefail

# Run this for the data directory

sf=../actions/irods.txt #tab separated file containing sampleIDs and corresponding irods path to directory
h5=raw_feature_bc_matrix.h5 #name of h5 file within irods directory

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

cat $sf | while read name path; do 
  mkdir -p $name
  ( cd $name
    echo "-- $path"
    iget -f -v -N 4 -K $path/$h5 raw_gene_bc_matrices_h5.h5; 
    echo "✓"
  )
done
