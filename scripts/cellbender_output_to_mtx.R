library(hdf5r)
library(Matrix)
# convert to mtx #########

sids = readLines('actions/sample-list')
for(sid in sids){
  if(!file.exists(paste0('work/',sid,'/cellbender_out.h5'))) next
  infile = hdf5r::H5File$new(filename = paste0('work/',sid,'/cellbender_out.h5'), mode = "r")
  n=names(infile)
  cnts = sparseMatrix(i=infile[[n]][['indices']]$read()+1, p=infile[[n]][['indptr']]$read(),x = as.numeric(infile[[n]][['data']]$read()),dims=infile[[n]][['shape']]$read())
  barcodes = infile[[n]][['barcodes']]$read()
  fnames = c('id','name','feature_type')
  features = as.data.frame(lapply(fnames,function(nn)infile[[n]][['features']][[nn]]$read()))
  colnames(features) = fnames
  infile$close_all()
  
  dir.create(paste0('work/',sid,'/cellbender_out'),showWarnings = FALSE)
  gz1 = gzfile(paste0('work/',sid,'/cellbender_out/features.tsv.gz'), "w")
  write.table(features,gz1,quote = FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  close(gz1)
  
  gz1 = gzfile(paste0('work/',sid,'/cellbender_out/barcodes.tsv.gz'), "w")
  writeLines(barcodes,gz1)
  close(gz1)
  
  if(file.exists(paste0('work/',sid,'/cellbender_out/matrix.mtx.gz'))) file.remove(paste0('work/',sid,'/cellbender_out/matrix.mtx.gz'))
  if(file.exists(paste0('work/',sid,'/cellbender_out/matrix.mtx'))) file.remove(paste0('work/',sid,'/cellbender_out/matrix.mtx'))
  writeMM(cnts,paste0('work/',sid,'/cellbender_out/matrix.mtx')) 
  system({paste0("gzip 'work/",sid,"/cellbender_out/matrix.mtx'")})
  
}

# QC ###############
# qc = read.table('/warehouse/cellgeni/tic-2055/work/qc.txt',header = TRUE)
# qc$cellsf = qc$common.detected.cells/(qc$probable.cells + qc$additional.cells)
# qc[qc$cellsf>0.99,]
# qc[order(qc$cellsf,decreasing = TRUE)[1:15],]
# 
# 
# plot(qc$cellsf,qc$median.soup.fraq)
# fisher.test(table(qc$cellsf>0.99,qc$median.soup.fraq<0.01))
