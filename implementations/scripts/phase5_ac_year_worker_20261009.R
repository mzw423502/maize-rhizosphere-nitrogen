args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L,args[1]%in%c("A_Y2018","A_Y2019","C_2017","C_2018"),
          args[2]%in%c("synthetic","actual"),args[3]%in%c("original","observed"),
          args[2]=="synthetic"||args[3]=="observed")
root <- "."
source(file.path(root,"scripts/phase5_ac_year_common_20261009.R"))
id <- args[1];kind <- args[2];variant <- args[3]
directory <- file.path(ac_base,kind,paste(id,variant,sep="_"))
stopifnot(dir.exists(directory),!file.exists(file.path(directory,"COMPLETION.json")))
options(contrasts=c("contr.treatment","contr.poly"),digits=17,warn=1)
suppressPackageStartupMessages(library(ANCOMBC))
stopifnot(as.character(getRversion())=="4.4.3",as.character(packageVersion("ANCOMBC"))=="2.8.0",
  as.character(packageVersion("lme4"))=="2.0.6",as.character(packageVersion("lmerTest"))=="3.2.1")
loaded <- ac_load(id);entry <- loaded$entry;count <- loaded$count;meta <- loaded$metadata
original_selected_digest <- digest::digest(count,algo="sha256")
original_zeros <- count==0
if(kind=="synthetic") {
  set.seed(entry$synthetic_positive_count_seed)
  count[] <- sample.int(9970L,length(count),replace=TRUE)+30L
  count[original_zeros] <- 0
}
write_new_json(list(kind=kind,variant=variant,model_id=id,dimensions=dim(count),families=rownames(count),
  sample_ids=colnames(count),zero_cells=sum(original_zeros),zero_layout_preserved=identical(count==0,original_zeros),
  input_digest=digest::digest(count,algo="sha256"),original_selected_digest=original_selected_digest,
  original_RDS_digest=loaded$original_digest,original_dimensions=loaded$original_dimensions,
  synthetic_positive_seed=if(kind=="synthetic")entry$synthetic_positive_count_seed else NULL,
  full_frozen_family_universe=loaded$input$families_in_order,
  not_estimable_families=entry$not_estimable_families,allocation_support_warning_only=entry$allocation_support_warning_only,
  pid=Sys.getpid(),utc=format(Sys.time(),tz="UTC",usetz=TRUE)),file.path(directory,"INPUT_IDENTITY.json"))
rm(loaded,original_zeros);gc()
observer <- if(variant=="observed") make_observed_ancom(directory) else NULL
call_model <- if(variant=="observed")observer$run else ANCOMBC::ancombc2
set.seed(entry$seed);warning_log <- character()
cat("AC_YEAR_WORKER_BEGIN model=",id," kind=",kind," utc=",format(Sys.time(),tz="UTC",usetz=TRUE),"\n",sep="")
fit <- withCallingHandlers(call_model(data=count,meta_data=meta,taxa_are_rows=TRUE,
  assay_name="counts",tax_level=NULL,aggregate_data=NULL,
  fix_formula=entry$fix_formula,rand_formula=entry$rand_formula,
  p_adj_method="BH",pseudo=0,pseudo_sens=TRUE,prv_cut=0,lib_cut=0,s0_perc=0.05,
  group="N",struc_zero=FALSE,neg_lb=FALSE,alpha=0.05,n_cl=1,verbose=TRUE,
  global=FALSE,pairwise=FALSE,dunnet=FALSE,trend=FALSE,
  iter_control=list(tol=0.01,max_iter=100,verbose=TRUE),
  em_control=list(tol=1e-5,max_iter=100),lme_control=lme4::lmerControl(),
  mdfdr_control=list(fwer_ctrl_method="holm",B=100),
  trend_control=list(contrast=NULL,node=NULL,solver="ECOS",B=100)),
  warning=function(w)warning_log <<- c(warning_log,conditionMessage(w)))
stopifnot(if(variant=="observed")observer$verify() else TRUE,
  identical(dimnames(fit$feature_table),dimnames(count)),
  isTRUE(all.equal(unname(fit$feature_table),unname(count),tolerance=0,check.attributes=FALSE)))
saveRDS(fit,file.path(directory,"NATIVE_RESULT.rds"),compress=FALSE,version=3)
write_new_json(list(status="AC_YEAR_EXECUTION_COMPLETE_NOT_SCIENTIFIC_RELEASE",model_id=id,kind=kind,variant=variant,
  utc=format(Sys.time(),tz="UTC",usetz=TRUE),samples=ncol(count),model_families=nrow(count),
  full_frozen_family_universe=entry$original_family_universe,fixed_formula=entry$fix_formula,random_formula=entry$rand_formula,
  seed=entry$seed,effective_iter_control=list(tol=0.01,max_iter=100,verbose=TRUE),warnings=warning_log,
  model_input_exactly_retained=TRUE,installed_functions_unchanged=TRUE,source_counts_modified=FALSE,
  original_positive_counts_used=kind=="actual",not_estimable_rows_not_reinterpreted_as_null_results=TRUE,
  inference="POST_QC_DESIGN_EXPLORATORY_ASSOCIATION",independent_N_units_not_sample_count=entry$unit_count,
  common_bias_feature_set_changed_by_explicit_human_support_rule=TRUE,
  cross_year_pooling=FALSE,Bonn_access=FALSE,scientific_release=FALSE),file.path(directory,"COMPLETION.json"))
cat("AC_YEAR_WORKER_COMPLETE\n")
