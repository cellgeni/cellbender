library(Matrix)
#library(optparse)
library(argparse)

# TODO 2) add thresholds to identify bad samples
# TODO 3) plot PCA/umap? on expression or on latent_gene_encoding?
# TODO 4) report cell with zero UMIs?



# option_list = list(
#   make_option(c("-i", "--input"), type="character", default='.',dest='dir', 
#               help="mandatory, path to folder that contains per-sample folders with CellBender results", metavar="character"),
#   make_option(c("-t", "--train-history-pdf"), type="character",default = NULL,dest = 'train.pdf', 
#               help="saves plots of training history into specified file (no plots are saved by default)", metavar="character")
# ) 
# 
# opt_parser = OptionParser(option_list=option_list,
#                           usage = "usage: %prog [options] dir",
#                           description='The programm collects summary statistics from CellBender outputs and prints it in stdout. The programm attempts to place most problematic samples first.')
# opt = parse_args2(opt_parser)


parser = ArgumentParser(description='The programm collects summary statistics from CellBender outputs and prints it in stdout. The programm attempts to place most problematic samples first.')
parser$add_argument("dir",nargs=1,help='path to folder that contains per-sample folders with CellBender results')
parser$add_argument("-m","--mode",help="run mode: 1 (default) is fast, it doesn't read count matrices, produce minimala output; 2 extracts cell probabilities from output h5, takes 2-10 sec per sample; 3 calculates soup fraction, takes 10-60 sec per sample",dest='mode',default=1,type='integer')
parser$add_argument("-o","--out",help='path to save results. -t, -c, -s and -r are defined relative to output folder. Default is currect folder.',default='.',dest='out')

parser$add_argument("-t","--train-history-pdf",help='saves plots of training history into specified file (train.pdf by default)',dest='train.pdf',default='train.pdf')
parser$add_argument("-c","--cell-probability-pdf",help='saves plots with cell probabilities into specified file, make sense only if MODE > 1 (prob.pdf by default)',dest='prob.pdf',default='prob.pdf')
parser$add_argument("-s","--soup-fraction-pdf",help='saves plots with soup-fraction vs total UMI plot into specified file. Only cells (that is barcodes with cell probability above 0.5) are considered. , make sense only if MODE > 2 (soup.pdf by default)',dest='soup.pdf',default='soup.pdf')
parser$add_argument("-r",'--rds',help='folder to save per-sample rds with summary information. Sample will not be reprocessed if corresponding rds is found in the folder (rds by default)',dest='rds',default='rds')
parser$add_argument("-n",help='process only first N folders (default is -1, process all)',dest='N',default=-1,type='integer')
parser$add_argument("-v","--verbose",help='Specifies whether programm should pring progress to stderr (silent by default). It makes sense to use this flag when mode is greater than 1 since reading of h5 files takes time.',action='store_true')
parser$add_argument("-w","--warnings",help='Specifies whether programm should pring warnings to stderr (silent by default).',action='store_true')


#opt = parser$parse_args(c("-p",'t.pdf','processed/multiome/gex/starsolo/cellbender.v03'))
opt = parser$parse_args()
opt$train.pdf = paste0(opt$out,'/',opt$train.pdf)
opt$prob.pdf = paste0(opt$out,'/',opt$prob.pdf)
opt$soup.pdf = paste0(opt$out,'/',opt$soup.pdf)
opt$rds = paste0(opt$out,'/',opt$rds)

parseCBlog = function(f){
  l = readLines(f)
  nloss = l[-grep('cellbender:remove-background: [epoch',l,fixed = T)]
  # parse loss
  lloss = l[grep('cellbender:remove-background: [epoch',l,fixed = T)]
  lloss=sub(' +\\(.+ seconds per epoch\\)$','',lloss)
  ls = as.numeric(stringr::str_match(lloss,'\\d+.\\d+$'))
  epoch = as.numeric(sub(']','',stringr::str_match(lloss,'\\d+\\]'),fixed = T))
  istest = grepl('average test loss',lloss)
  training = data.frame(epoch=sort(unique(epoch)),test=NA,train=NA)
  
  training$test[match(epoch[istest],training$epoch)] = -ls[istest]
  training$train[match(epoch[!istest],training$epoch)] = -ls[!istest]

    # parse other things
  patterns = c(prior.empty.count = "Prior on counts in empty droplets is ",
               prior.cell.count = "Prior on counts for cells is ",
               barcode.low.thr = "Excluding barcodes with counts below ",
               posterior.reg.factor = "Optimal posterior regularization factor = ",
               largest.empty.count = "Largest surely-empty droplet has ",
               genes.used = "Including ",
               probable.cells = "Using ")
  info = as.data.frame(lapply(patterns,function(p){
    t = nloss[grep(p,nloss)]
    r = stringr::str_match(t,paste0(p,"\\d+.?\\d*"))
    r = as.numeric(sub(p,'',r))
    if(length(r)==0)
      r = NA
    r
  }))  
  t = nloss[grep('probable cell barcodes, plus an additional ',nloss)]
  z = 'probable cell barcodes, plus an additional '
  info$additional.cells = as.numeric(sub(z,'',stringr::str_match(t,paste0(z,"\\d+.?\\d*"))))
  z=' barcodes, and '
  info$probable.empty.cells = as.numeric(sub(z,'',stringr::str_match(t,paste0(z,"\\d+.?\\d*"))))
  # find input and output
  io = strsplit(l[grep("cellbender remove-background --input",l)],' ')[[1]]
  input = io[which(io=='--input')+1]
  output = io[which(io=='--output')+1]
  
  
  info$input = sub(output,input,gsub('log$','h5',normalizePath(f)))
  return(list(training=training,info=info))
}

plotCBtraining = function(f1,...){
  plot(f1$training$epoch,f1$training$train,t='l',col='blue',ylim=range(f1$training$test,f1$training$train,na.rm=T),lwd=3,xlab='Epoch',ylab='ELBO',...)
  f = !is.na(f1$training$test)
  lines(f1$training$epoch[f],f1$training$test[f],col='orange',lwd=3)
  legend('bottomright',lwd=3,col=c('blue','orange'),c('Train','Test'),bty='n')
}

plotCellProbability = function(l,...){
  if(is.null(l$cell.probability))
    plot.new()
  else{
    plot(l$cell.probability,type='l',xlab='Barcode index, sorted by UMI count',ylab='Cell probability',ylim=c(0,1),...)
  }
}

parseCB.h5 = function(f){
  infile = hdf5r::H5File$new(filename = f, mode = "r")
  n=names(infile)
  #str(sapply(names(infile[[n]]),function(m){infile[[n]][[m]]$maxdims}))
  r = list()
  r$cell.probability = infile[[n]][['latent_cell_probability']]$read()
  names(r$cell.probability) = infile[[n]][['barcodes']]$read()[infile[[n]][['barcode_indices_for_latents']]$read()+1]
  
  training = infile[[n]][['training_elbo_per_epoch']]$read()
  training = data.frame(epoch=1:length(training),test=NA,train=training)
  
  training$test[match(infile[[n]][['test_epoch']]$read(),training$epoch)] = infile[[n]][['test_elbo']]$read()
  
  z = sparseMatrix(i=infile[[n]][['indices']]$read()+1, p=infile[[n]][['indptr']]$read(),x = as.numeric(infile[[n]][['data']]$read()),dims=infile[[n]][['shape']]$read())
  r$filtered.tUMI = setNames(colSums(z),infile[[n]][['barcodes']]$read())
  infile$close_all()
  r$training = training
  r
}

plot.as.hm = function(x,y,xbins=100,ybins=100,cols=c('white','gray','blue','orange','red'),zfun=identity,leg.title='',num.leg.tic=NULL,legend=TRUE,trimZq=0,xlim=NULL,ylim=NULL,new=TRUE,xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),...){
  xlab;ylab;
  f = !is.na(x) & !is.infinite(x) & !is.na(y) & !is.infinite(y) 
  if(!is.null(xlim)) f = f & x >= xlim[1] & x <= xlim[2]
  if(!is.null(ylim)) f = f & y >= ylim[1] & y <= ylim[2]
  x = x[f]
  y = y[f]
  if(length(xbins)==1) xbins = seq(min(x),max(x),length.out = xbins+1)
  if(length(ybins)==1) ybins = seq(min(y),max(y),length.out = ybins+1)
  
  f = x >= xbins[1] & x <= xbins[length(xbins)] & y >= ybins[1] & y <= ybins[length(ybins)]
  x = x[f]
  y = y[f]
  
  xb = findInterval(x,xbins,rightmost.closed = TRUE)
  yb = findInterval(y,ybins,rightmost.closed = TRUE)
  
  
  z = as.matrix(table(factor(yb,levels=1:(length(ybins)-1)),factor(xb,levels=1:(length(xbins)-1))))
  z = trimQ(z,trimZq)
  if(is.null(xlim)) xlim = range(x)
  if(is.null(ylim)) ylim = range(y)
  if(new)
    plot(1,t='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  z2col=function(x)num2col(zfun(x),cols)
  rect(rep(xbins[-length(xbins)],each =length(ybins)-1),
       rep(ybins[-length(ybins)],times=length(xbins)-1),
       rep(xbins[-1]            ,each =length(ybins)-1),
       rep(ybins[-1]            ,times=length(xbins)-1),
       col=z2col(z),border = NA)
  if(legend)
    plotColorLegend2(grconvertX(1,'npc','nfc'),1,grconvertY(0,'npc','nfc'),grconvertY(1,'npc','nfc'),
                     range(z),range(z),zfun,z2col=z2col,leg=num.leg.tic,title=leg.title)
  invisible(list(xbins=xbins,ybins=ybins,z=z))
}

trimQ = function(x,q){
  if(q==0) return(x)
  qq=quantile(x,sort(c(q,1-q)))
  x[x<=qq[1]] = qq[1]
  x[x>=qq[2]] = qq[2]
  x
}

num2col = function(d,col=c('blue','cyan','gray','orange','red'),minx=min(d),maxx=max(d)){
  if(sd(d)==0)
    return(rep(col[1],length(d)))
  d[d<minx] = minx
  d[d>maxx] = maxx
  bwr = colorRamp(col,alpha=TRUE)
  apply(bwr(scaleTo(d,0,1,minx=minx,maxx=maxx))/255,1,function(x)rgb(x[1],x[2],x[3],x[4]))
}

plotColorLegend2 = function(x0,x1,y0,y1,fullzlim,zlim,zfun,z2col,N=100,ntic=5,leg=NULL,title=NULL){
  # make tics
  if(is.null(leg)){
    ztic = seq(zlim[1],zlim[2],length.out = 1e5)
    ztict = zfun(ztic)
    ztticat = seq(zfun(zlim[1]),zfun(zlim[2]),length.out = ntic)
    leg = ztic[findInterval(ztticat,ztict)]
    leg[1] = zlim[1]
    leg[ntic] = zlim[2]
    
    # adjust tic step
    d = (10^floor(log10(leg[2]-leg[1])))/2
    leg = unique(d*round(leg / d))
    leg = leg[leg>=zlim[1]]
  }
  # make colors
  ztat = zfun(leg)
  ztcol = sort(unique(c(ztat,seq(zfun(zlim[1]),zfun(zlim[2]),length.out = N))))
  col = z2col(c(zfun(fullzlim),ztcol))[-(1:length(fullzlim))]
  at = match(ztat,ztcol)
  plotColorLegend(x0,x1,y0,y1,col,at=at,legend=leg,title=title)
}

plotColorLegend = function(x0,x1,y0,y1,col,at,legend,title=NULL){
  xpd = par(xpd=TRUE)
  y = seq(grconvertY(y0,'nfc','user'),grconvertY(y1,'nfc','user'),length.out = length(col)+1)
  rect(grconvertX(x0,'nfc','user'),y[-length(y)],grconvertX(x0+(x1-x0)*0.25,'nfc','user'),y[-1],col=col,border = NA)
  at = y[at]+(y[2]-y[1])/2
  text(grconvertX(x0+(x1-x0)*0.3,'nfc','user'),at,legend,adj=c(0,0.5))
  if(!is.null(title)){
    #text(grconvertX(x1,'nfc','user'),y[length(y)],title,adj=c(1,-0.5))
    text(grconvertX(x0,'nfc','user'),y[length(y)],title,adj=c(0,-0.5))
  }
  par(xpd=xpd)
}

log10p1 = function(x)log1p(x)/log(10)

scaleTo = function(x,from=0,to=1,minx=min(x,na.rm=TRUE),maxx=max(x,na.rm=TRUE),fraction=1){
  x = (x-minx)/(maxx-minx)
  x*(to-from)*fraction + from + (to-from)*(1-fraction)/2
}

plotSoupFraq = function(l,...){
  cbc = names(l$cell.probability)[l$cell.probability>0.5]
  x = log10p1(l$raw.tUMI[cbc])
  y = (l$raw.tUMI[cbc]-l$filtered.tUMI[cbc])/l$raw.tUMI[cbc]*100
  plot.as.hm(x,y,zfun = log1p,xlab='log10(UMI + 1)',ylab='soup %',leg.title = '#cells',ylim=c(0,100),...)
}

myRead10X = function(f){
  require(Matrix)
  if(file.exists(paste0(f,'/matrix.mtx'))){
    m = Matrix::readMM(paste0(f,'/matrix.mtx'))
    rownames(m) = readLines(paste0(f,'/features.tsv'))
    colnames(m) = readLines(paste0(f,'/barcodes.tsv'))
  }else{
    m = Matrix::readMM(paste0(f,'/matrix.mtx.gz'))
    rownames(m) = readLines(paste0(f,'/features.tsv.gz'))
    colnames(m) = readLines(paste0(f,'/barcodes.tsv.gz'))
  }
  m
}

if(!dir.exists(opt$dir))
  stop('Folder "',opt$dir,'" does not exist.')

sids = list.dirs(opt$dir,recursive = FALSE,full.names = FALSE)

if(!dir.exists(opt$rds)) dir.create(opt$rds)
# parse logs 
logs = list()
n = 0
start.time = Sys.time() 
for(i in 1:length(sids)){
  if(opt$N != -1 & n == opt$N)
    break
  s = sids[i]
  sample.rds = paste0(opt$rds,'/',s,opt$mode,'.rds')
  if(file.exists(sample.rds)){
    logs[[s]] = readRDS(sample.rds)
    if(opt$verbose)
      cat('INFO: ',s,', pre-saved rds found; ',n,' folders processed, ',i-1,' folders (out of ',length(sids),') checked. Time elapsed: ',format(Sys.time()  - start.time),'\n',file =  stderr(),sep='')
    n = n+1
    next
  }
  
  if(opt$verbose)
    cat('INFO: ',s,'; ',n,' folders processed, ',i-1,' folders (out of ',length(sids),') checked. Time elapsed: ',format(Sys.time()  - start.time),'\n',file =  stderr(),sep='')
  
  log.file = list.files(paste0(opt$dir,'/',s),pattern = '.log',recursive = TRUE,full.names = TRUE)
  if(length(log.file) != 1){
    if(opt$warn)
      cat("WARN: there are ",length(log.file),' log files in folder "',s,'", while one is expected.\n',file =  stderr(),sep='')
    next
  }
  
  file.h5 = sub('.log','.h5',log.file)
  if(!file.exists(file.h5)){
    cat("ERROR: file ",file.h5,' does not exist. Sample "',s,'" skipped\n',file =  stderr(),sep='')
    next
  }
  
  
  l = parseCBlog(log.file)
  last.notna = max(which(!is.na(l$training$test)))
  
  l$info$final.test.elbo.minus.max = l$training$test[last.notna] - max(l$training$test,na.rm=TRUE)
  l$info$test.minus.train.elbo = l$training$test[last.notna ]-l$training$train[last.notna]
  l$info$epochs = max(l$training$epoch)
  
  # check called cells
  cbbcs = readLines(sub('.log','_cell_barcodes.csv',log.file))
  l$info$cellbender.detected.cells = length(cbbcs)
  l$info$input.detected.cells = NA
  l$info$common.detected.cells = NA
  try({
    inbcs = readLines(list.files(gsub('raw','filtered/',l$info$input),pattern = '*barcodes*',full.names = T,recursive = TRUE))
    l$info$input.detected.cells = length(inbcs)
    l$info$common.detected.cells = length(intersect(inbcs,cbbcs))
  },silent = TRUE)
  
  
  if(opt$mode>1){ # reading of h5 file take more time, so I'll read it only if requested
    h5 = parseCB.h5(file.h5)
    l$training.h5 = h5$training # it is identical to training curve readed from log, but it is more explicit
    l$cell.probability = h5$cell.probability
    l$filtered.tUMI = h5$filtered.tUMI
    #l$info$detected.cells2 = sum(l$cell.probability>0.5)
    l$info$top10umi.min.cell.prob = min(l$cell.probability[1:10])
  }
  if(opt$mode>2){
    r = myRead10X(l$info$input)
    l$raw.tUMI = colSums(r)[names(l$filtered.tUMI)]
    cbc = names(l$cell.probability)[l$cell.probability>0.5]
    soup.fraq = (l$raw.tUMI-l$filtered.tUMI)/l$raw.tUMI
    soup.fraq = soup.fraq[cbc]
    l$info$mean.soup.fraq = mean(soup.fraq)
    l$info$median.soup.fraq = median(soup.fraq)
  }
  
  if(!file.exists(sample.rds)){
    saveRDS(l,sample.rds)
  }
  logs[[s]] = l
  
  # preliminary stop
  n = n + 1
}

# extract summary
info = do.call(rbind,lapply(logs,function(l)l$info))
info = cbind(sample.id=rownames(info),info)

if(nrow(info) == 0)
  stop("No CellBender output were found.")

# order samples
o = order(info$final.test.elbo.minus.max)
if(!is.null(info$top10umi.min.cell.prob)){
  o = order(info$top10umi.min.cell.prob>0.5,info$final.test.elbo.minus.max)
}
info = info[o,]
logs = logs[o]

# write summary
write.table(info,quote = FALSE,sep = '\t',row.names = FALSE)

# plot training history

pdf(opt$train.pdf,w=6*3,h=4*3)
par(mfrow=c(4,6),tcl=-0.2,mgp=c(1.1,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,1),bty='n')
for(i in names(logs)){
  plotCBtraining(logs[[i]],main=i)
}
t=dev.off()


# plot cell probabilities
if(opt$mode>1){
  pdf(opt$prob.pdf,w=6*3,h=4*3)
  par(mfrow=c(4,6),tcl=-0.2,mgp=c(1.1,0.2,0),mar=c(2.5,2.5,1.5,0),oma=c(0,0,1,1),bty='n')
  for(i in names(logs)){
    plotCellProbability(logs[[i]],main=i)
  }
  t=dev.off()
}

# plot soup fraq
if(opt$mode>2){
  pdf(opt$soup.pdf,w=6*3.5,h=4*3)
  par(mfrow=c(4,6),tcl=-0.2,mgp=c(1.1,0.2,0),mar=c(2.5,2.5,1.5,4),oma=c(0,0,1,1),bty='n')
  for(i in names(logs)){
    plotSoupFraq(logs[[i]],main=i)
  }
  t=dev.off()
}



