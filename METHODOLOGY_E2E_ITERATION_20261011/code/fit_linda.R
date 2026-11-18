# Execute the pinned current upstream function unchanged. All 100 taxa retained.
args <- commandArgs(trailingOnly=TRUE); scale <- args[[1]]
stopifnot(scale %in% c('pilot','formal'))
.libPaths(c(file.path(getwd(),'runtime','R-library'),.libPaths()))
suppressPackageStartupMessages(library(modeest))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(matrixStats))
source('upstream/MicrobiomeStat_R_linda.R')
manifest <- jsonlite::fromJSON(file.path('data',scale,'MANIFEST.json'))
dest <- file.path('estimates',scale);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
for (i in seq_len(nrow(manifest$datasets))) {
  row <- manifest$datasets[i,];p <- row$path
  con <- gzfile(p,'rt');dat <- jsonlite::fromJSON(paste(readLines(con,warn=FALSE),collapse=''));close(con)
  g <- dat$parameters$pairs;s <- dat$subsamples;m <- dat$taxa
  meta <- data.frame(arm=factor(rep(c('0','1'),g),levels=c('0','1')),pair=factor(rep(seq_len(g),each=2)))
  result <- list(cell=dat$cell,replicate=dat$replicate,input_sha256=row$sha256,engine='MicrobiomeStat 1.4 unmodified LinDA source',years=list())
  for (year in c('A1','A2','B')) {
    counts <- dat[[year]]$counts
    W <- log(counts+0.5);W <- W-rowMeans(W)
    U <- rowsum(W,rep(seq_len(2*g),each=s),reorder=FALSE)/s
    P <- exp(U-apply(U,1,max));P <- t(P/rowSums(P))
    rownames(P) <- sprintf('T%03d',seq_len(m));colnames(P) <- sprintf('U%03d',seq_len(2*g));rownames(meta)<-colnames(P)
    warnings <- character()
    fit <- withCallingHandlers(linda(feature.dat=P,meta.dat=meta,formula='~arm+pair',feature.dat.type='proportion',prev.filter=0,mean.abund.filter=0,max.abund.filter=0,is.winsor=FALSE,adaptive=FALSE,zero.handling='pseudo-count',pseudo.cnt=.5,p.adj.method='BH',alpha=.05,n.cores=1,verbose=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    stopifnot(identical(rownames(fit$output$arm1),rownames(P)),length(fit$output$arm1$log2FoldChange)==100)
    treatment_row <- which(rownames(fit$wald$beta)=='arm1')
    nclr <- t(log(fit$feature.dat.use))-colMeans(log(fit$feature.dat.use))
    result$years[[year]] <- list(native=fit$output$arm1,raw_coefficient_log2=unname(fit$wald$beta[treatment_row,]),bias_log2=unname(fit$bias[which(fit$variables=='arm1')]),unit_clr=U,unit_profile=P,native_clr=nclr,profile_CLR_maxabs=max(abs(U-nclr)),warnings=warnings,df=2*g-ncol(fit$wald$X))
  }
  target <- file.path(dest,sub('.json.gz','.json',basename(p),fixed=TRUE));stopifnot(!file.exists(target))
  jsonlite::write_json(result,target,digits=17,auto_unbox=TRUE,pretty=FALSE,na='null')
  if(i %% 10 == 0) {cat(sprintf('completed %d / %d\n',i,nrow(manifest$datasets)));flush.console();gc(FALSE)}
}
jsonlite::write_json(list(scale=scale,datasets=nrow(manifest$datasets),complete=TRUE,source='upstream/CURRENT_LINDA_MANIFEST.json',modeest=as.character(packageVersion('modeest')),matrixStats=as.character(packageVersion('matrixStats')),R=R.version.string),file.path(dest,'COMPLETION.json'),pretty=TRUE,auto_unbox=TRUE)
