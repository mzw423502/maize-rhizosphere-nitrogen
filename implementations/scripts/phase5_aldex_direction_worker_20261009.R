args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L,args[1] %in% c("synthetic","actual"))
kind<-args[1];id<-args[2]
source("./scripts/phase5_aldex_direction_common_20261009.R")
suppressPackageStartupMessages(library(ALDEx2))
al_native_guard()
input<-al_load(id)
folder<-file.path(al_base,kind,id);stopifnot(dir.exists(folder),!file.exists(file.path(folder,"COMPLETION.json")))
if(kind=="synthetic") {
  # Complete bounded synthetic tables, never a partial feature/denominator call.
  # Several zero layouts exercise both cohort-total-zero and sample-zero behavior.
  checks<-list()
  for(case in seq_len(3L)) {
    nf<-input$entry$universe;nc<-min(64L,input$entry$samples)
    v<-seq_len(nf*nc);x<-matrix(as.double((v*37+case*11)%%1001+1),nf,nc)
    x[(v+case)%%7==0]<-0
    x[,1L]<-seq_len(nf) # Every intended nonzero feature is in the complete table.
    if(input$entry$cohort_total_zero>0L)x[seq.int(nf-input$entry$cohort_total_zero+1L,nf),]<-0
    if(case==2L)x[,nc]<-0 # Prior-only sample retained; no column-local feature deletion.
    if(case==3L)x[seq_len(nf- input$entry$cohort_total_zero),2L]<-1
    dimnames(x)<-list(paste0("fixture_f",seq_len(nf)),paste0("fixture_s",seq_len(nc)))
    condition<-rep(c("HN","LN"),length.out=nc)
    w<-ifelse(condition=="LN",1/sum(condition=="LN"),-1/sum(condition=="HN"))
    al_rng(input$entry$seed+case)
    native<-aldex.clr(x,condition,mc.samples=128L,denom="all",gamma=NULL,useMC=FALSE,verbose=FALSE)
    end_state<-.Random.seed
    native_p<-getDirichletInstances(native);native_z<-getMonteCarloInstances(native)
    checker<-function(i,p,z)stopifnot(identical(p,native_p[[i]]),identical(z,native_z[[i]]))
    testpath<-file.path(folder,paste0("NATIVE_EQUIVALENCE_CASE_",case,".bin"))
    stream<-al_stream(x,condition,w,rep("fixture_pair",nc),"fixture_pair",input$entry$seed+case,testpath,checker)
    stopifnot(identical(stream$final_rng,end_state))
    direct<-Reduce(`+`,Map(function(z,w)z*w,native_z,as.list(w)))
    stopifnot(identical(unname(direct),unname(stream$aggregate)))
    checks[[case]]<-list(case=case,features=length(stream$features),samples=nc,mc=128L,
       exact_dirichlet_and_clr=TRUE,exact_aggregate=TRUE,exact_final_rng=TRUE,
       compared_clr_cells=length(stream$features)*nc*128L)
    rm(native,native_p,native_z,stream,direct,x);invisible(gc())
  }
  al_json(list(status="NATIVE_FULL_BOUNDED_TABLE_EXACT_EQUIVALENCE_PASS",checks=checks,
     actual_effect_calls=0L,native_largest_full_cohort_RSS_measured=FALSE,
     proof_scope="Complete synthetic table for each case; full cohort scaling is tested separately below, not claimed as a full cohort native run."),file.path(folder,"NATIVE_EQUIVALENCE.json"))
  full<-input$full;pos<-which(full>0)
  full[pos]<-as.double((seq_along(pos)*37+11)%%10001+1)
  stopifnot(identical(full==0,input$full==0))
} else full<-input$full
result<-al_stream(full,input$meta$nitrogen_treatment,input$coefficient,input$meta$spatial_pair_year,
                 input$pairs,input$entry$seed,file.path(folder,"CLR_STREAM.bin"))
al_save_result(result,full,input,folder,kind)
al_json(list(status="ALDEX_DIRECTION_STREAM_COMPLETED_NUMERIC_AUDIT_PENDING",utc=format(Sys.time(),tz="UTC",usetz=TRUE),
  kind=kind,cohort=id,features=length(result$features),samples=length(result$samples),mc=128L,
  R_version=as.character(getRversion()),ALDEx2_version=as.character(packageVersion("ALDEx2")),
  biological_inference_computed=FALSE,scientific_release=FALSE),file.path(folder,"COMPLETION.json"))
cat("ALDEX_DIRECTION_STREAM_COMPLETE",kind,id,"\n")
