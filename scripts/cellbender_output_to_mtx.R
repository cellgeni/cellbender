library(hdf5r)
library(Matrix)
# cellbender v2 outputs invalid h5 if it starts from mtx (no genome in header) so scanpy cannot read it
# this script converts h5 into mtx format so anyone can read it
# find all cellbender outputs in current folder recursively
# convert to mtx #########

# function takes cellbender outdir as input and creates subfolder there with unfiltered counts in mtx format
h5_to_mtx = function(cb_out,outdir=paste0(cb_out,'/cellbender_out')){
  h5file = paste0(cb_out,'/cellbender_out.h5')
  infile = hdf5r::H5File$new(filename = h5file, mode = "r")
  n=names(infile)
  # counts
  cnts = sparseMatrix(i=infile[[n]][['indices']]$read()+1, p=infile[[n]][['indptr']]$read(),x = as.numeric(infile[[n]][['data']]$read()),dims=infile[[n]][['shape']]$read())
  # barcodes
  barcodes = infile[[n]][['barcodes']]$read()
  # features
  fnames = c('id','name','feature_type')
  features = as.data.frame(lapply(fnames,function(nn)infile[[n]][['features']][[nn]]$read()))
  
  # save; I'll not delete folder if exists, just overwrite
  dir.create(outdir,showWarnings = FALSE)
  
  gz1 = gzfile(paste0(outdir,'/features.tsv.gz'), "w")
  write.table(features,gz1,quote = FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
  close(gz1)
  
  gz1 = gzfile(paste0(outdir,'/barcodes.tsv.gz'), "w")
  writeLines(barcodes,gz1)
  close(gz1)
  
  writeMM(cnts,paste0(outdir,'/matrix.mtx')) 
  system({paste0("gzip -f '",outdir,"/matrix.mtx'")})
}


cb_outs = list.files('.',recursive = TRUE,pattern = 'cellbender_out.h5')
cb_outs = sub('cellbender_out.h5$','',cb_outs)
for(o in cb_outs){
  print(o)
  h5_to_mtx(o)
}
