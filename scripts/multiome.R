library(Matrix)

subsetGEX = function(d,out,type='Gene Expression'){
  m = readMM(paste0(d,'/matrix.mtx.gz'))
  features = read.table(paste0(d,'/features.tsv.gz'),sep='\t')
  if(!dir.exists(out)) dir.create(out)
  file.copy(paste0(d,'/barcodes.tsv.gz'),paste0(out,'/barcodes.tsv.gz'))
  f = features$V3 == type
  gz1 = gzfile(paste0(out,'/features.tsv.gz'), "w")
  write.table(features[f,],gz1,quote = FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  close(gz1)
  
  if(file.exists(paste0(out,'/matrix.mtx.gz'))) file.remove(paste0(out,'/matrix.mtx.gz'))
  m = m[f,]  
  writeMM(m,paste0(out,'/matrix.mtx')) # it doesn't work (or was extremely slow) with gzfile probably because of files sizes 
  system({paste0("gzip ",out,'/matrix.mtx')})
}

theargs <- R.utils::commandArgs(asValues=TRUE)

path <- theargs$path

subsetGEX(paste0(path,"/raw_feature_bc_matrix/"), paste0(path,"/raw_feature_bc_matrix_gex/"))
