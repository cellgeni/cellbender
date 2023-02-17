library(Seurat)  
library(Matrix)

# this scripts goes through subfolders in data, extracts gene expression part of multiome 
# and saves it into gex subfolder as matrix.mtx.gz accompained by features.tsv.gz and barcodes.tsv.gz
# it should works with both h5 and mtx.gz input (but I didn't test it with the second)
# Pasha M

subsetGEX = function(d,out,type='Gene Expression'){
  print(d)
  if(file.exists(paste0(d,'/matrix.mtx.gz'))){
    m = readMM(paste0(d,'/matrix.mtx.gz'))
    features = read.table(paste0(d,'/features.tsv.gz'),sep='\t')
    barcodes = readLines(paste0(d,'/barcodes.tsv.gz'))
    f = features$V3 == type
    m = m[f,]
    features = features[f,]
  }else{
    m=Seurat::Read10X_h5(d,unique.features = FALSE,use.names = FALSE)$`Gene Expression`
    n=Seurat::Read10X_h5(d,unique.features = FALSE,use.names = TRUE)$`Gene Expression`
    features = data.frame(rownames(m),rownames(n),'Gene Expression')
    barcodes = colnames(m)
  }
  
  if(!dir.exists(out)) dir.create(out)
  
  gz1 = gzfile(paste0(out,'/features.tsv.gz'), "w")
  write.table(features,gz1,quote = FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  close(gz1)
  
  gz1 = gzfile(paste0(out,'/barcodes.tsv.gz'), "w")
  writeLines(barcodes,gz1)
  close(gz1)
  
  if(file.exists(paste0(out,'/matrix.mtx.gz'))) file.remove(paste0(out,'/matrix.mtx.gz'))
  writeMM(m,paste0(out,'/matrix.mtx')) # it doesn't work (or was extremely slow) with gzfile probably because of files sizes 
  system({paste0("gzip ",out,'/matrix.mtx')})
}

paths = list.files('data/',pattern = 'raw_feature_bc_matrix.h5',recursive = TRUE,full.names = TRUE)
if(length(paths) == 0){
  paths = list.files('data/',pattern = 'matrix.mtx.gz',recursive = TRUE,full.names = TRUE)
  paths = dirname(paths)
}
for(path in paths)
  subsetGEX(path, paste0(dirname(path),'/gex'))
