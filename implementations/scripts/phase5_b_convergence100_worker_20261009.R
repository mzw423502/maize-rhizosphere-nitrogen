args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L,args[1] %in% c("primary","prefix_sensitivity"),
          args[2] %in% c("synthetic","actual"),args[3] %in% c("original","observed"))
stopifnot(identical(args,c("primary","actual","observed")))
kind <- args[2]; variant <- args[3]
stopifnot(kind=="synthetic" || variant=="observed")
id <- args[1]
root <- "."
directory <- file.path(root,"15_FINAL_AUDIT/PHASE5_B_PRIMARY_CONVERGENCE100_20261009_V1",kind,paste(id,variant,sep="_"))
stopifnot(dir.exists(directory),!file.exists(file.path(directory,"COMPLETION.json")),
          !file.exists(file.path(directory,"FAILURE.json")))
source(file.path(root,"scripts/phase5_b_random_model_common_20261009.R"))
options(contrasts=c("contr.treatment","contr.poly"),digits=17,warn=1)
suppressPackageStartupMessages(library(ANCOMBC))
stopifnot(as.character(getRversion())=="4.4.3",as.character(packageVersion("ANCOMBC"))=="2.8.0",
          as.character(packageVersion("lme4"))=="2.0.6",as.character(packageVersion("lmerTest"))=="3.2.1")
contract <- jsonlite::fromJSON(file.path(root,"15_FINAL_AUDIT/PHASE5_B_CONVERGENCE100_CONTRACT_20261009.json"),simplifyVector=FALSE)
entry <- Filter(function(x)x$id==id,contract$models_in_execution_order)[[1]]
count <- readRDS(file.path(root,"15_FINAL_AUDIT/PHASE5_NATIVE_COUNT_INPUT_20261009_V1/B_COUNTS.rds"))
meta <- jsonlite::fromJSON(file.path(root,"15_FINAL_AUDIT/PHASE5_B_RANDOM_MODEL_INPUTS_20261009_V1/B_MODEL_METADATA.json"))
stopifnot(identical(dim(count),c(119L,2123L)),all(is.finite(count)),all(count>=0),all(count==floor(count)),
          identical(meta$sample_id,colnames(count)),!anyDuplicated(rownames(count)),!anyDuplicated(colnames(count)))
rownames(meta) <- meta$sample_id
meta$N <- factor(meta$N,levels=c("HN","LN"))
meta$spatial_pair <- factor(meta$spatial_pair,levels=c("1","2"))
meta$allocation_unit <- factor(meta$allocation_unit,levels=c("1","2","3","4"))
meta$genotype <- factor(meta$genotype,levels=sort(unique(meta$genotype)))
meta$library_prefix <- factor(meta$library_prefix,levels=paste0("D",14:21))
stopifnot(!anyNA(meta),nlevels(meta$genotype)==528L)
original_zeros <- count==0
if (kind=="synthetic") {
  set.seed(2026100903)
  count[] <- sample.int(9970L,length(count),replace=TRUE)+30L
  count[original_zeros] <- 0
  stopifnot(identical(count==0,original_zeros),all(count[!original_zeros]>0))
}
write_new_json(list(kind=kind,variant=variant,dimensions=dim(count),families=rownames(count),
  sample_ids=colnames(count),zero_cells=sum(original_zeros),zero_layout_preserved=identical(count==0,original_zeros),
  generated_positive_seed=if(kind=="synthetic")2026100903 else NULL,
  synthetic_input_digest=digest::digest(count,algo="sha256"),
  pid=Sys.getpid(),utc=format(Sys.time(),tz="UTC",usetz=TRUE)),file.path(directory,"INPUT_IDENTITY.json"))
rm(original_zeros)
observer <- if(variant=="observed") make_observed_ancom(directory) else NULL
call_model <- if(variant=="observed") observer$run else ANCOMBC::ancombc2
set.seed(entry$seed)
warning_log <- character()
cat("B_LOW_MEMORY_WORKER_BEGIN model=",id," utc=",format(Sys.time(),tz="UTC",usetz=TRUE),"\n",sep="")
fit <- withCallingHandlers(call_model(data=count,meta_data=meta,taxa_are_rows=TRUE,
  assay_name="counts",tax_level=NULL,aggregate_data=NULL,
  fix_formula=entry$fix_formula,rand_formula=contract$rand_formula,
  p_adj_method="BH",pseudo=0,pseudo_sens=TRUE,prv_cut=0,lib_cut=0,s0_perc=0.05,
  group="N",struc_zero=FALSE,neg_lb=FALSE,alpha=0.05,n_cl=1,verbose=TRUE,
  global=FALSE,pairwise=FALSE,dunnet=FALSE,trend=FALSE,
  iter_control=list(tol=0.01,max_iter=100,verbose=TRUE),
  em_control=list(tol=1e-5,max_iter=100),lme_control=lme4::lmerControl(),
  mdfdr_control=list(fwer_ctrl_method="holm",B=100),
  trend_control=list(contrast=NULL,node=NULL,solver="ECOS",B=100)),
  warning=function(w)warning_log <<- c(warning_log,conditionMessage(w)))
stopifnot(if(variant=="observed") observer$verify() else TRUE,identical(dimnames(fit$feature_table),dimnames(count)),
          isTRUE(all.equal(unname(fit$feature_table),unname(count),tolerance=0,check.attributes=FALSE)))
saveRDS(fit,file.path(directory,"NATIVE_RESULT.rds"),compress=FALSE,version=3)
write_new_json(list(status="B_WORKER_EXECUTION_COMPLETE_NOT_RELEASE",model=id,kind=kind,variant=variant,
  utc=format(Sys.time(),tz="UTC",usetz=TRUE),samples=2123L,families=119L,
  fixed_formula=entry$fix_formula,random_formula=contract$rand_formula,seed=entry$seed,
  effective_iter_control=list(tol=0.01,max_iter=100,verbose=TRUE),
  warnings=warning_log,installed_function_bodies_unchanged=TRUE,
  worker_count_table_exactly_retained=TRUE,original_positive_counts_used=(kind=="actual"),inference="EXPLORATORY_ASSOCIATION_ONLY",
  library_proxy_not_certified_batch=TRUE,independent_N_allocation_units=4L,
  Bonn_access=FALSE,scientific_release=FALSE),file.path(directory,"COMPLETION.json"))
cat("B_LOW_MEMORY_WORKER_COMPLETED model=",id,"\n",sep="")
