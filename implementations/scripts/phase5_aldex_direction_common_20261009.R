# New directional adapter only. Frozen ALDEx2 functions are called, never altered.
al_root <- "."
al_base <- file.path(al_root,"15_FINAL_AUDIT/PHASE5_ALDEX2_UNIT_DIRECTION_20261009_V1")
al_contract_path <- file.path(al_root,"15_FINAL_AUDIT/PHASE5_ALDEX2_DIRECTION_CONTRACT_20261009.json")
al_json <- function(x,path) {
  stopifnot(!file.exists(path))
  jsonlite::write_json(x,path,pretty=TRUE,auto_unbox=TRUE,digits=NA,null="null",na="null")
}
al_entry <- function(id) {
  c <- jsonlite::fromJSON(al_contract_path,simplifyVector=FALSE)
  e <- Filter(function(x)identical(x$id,id),c$methods$cohort_plan)
  stopifnot(length(e)==1L,c$methods$mc_samples==128L,identical(c$methods$denom,"all"),is.null(c$methods$gamma))
  e[[1]]
}
al_native_guard <- function() {
  ns <- loadNamespace("ALDEx2")
  stopifnot(as.character(packageVersion("ALDEx2"))=="1.38.0")
  stopifnot(as.character(getRversion())=="4.4.3",Sys.getenv("OMP_NUM_THREADS")=="1",Sys.getenv("OPENBLAS_NUM_THREADS")=="1")
  ref <- parse(file.path(al_root,"15_FINAL_AUDIT/PHASE5_ALDEX2_SOURCE_REFERENCE_20261009/clr_function.r"),keep.source=FALSE)
  ex <- Filter(function(x)is.call(x)&&identical(x[[1]],as.name("<-"))&&identical(x[[2]],as.name("aldex.clr.function")),as.list(ref))
  original <- eval(ex[[1]][[3]],envir=baseenv())
  installed <- get("aldex.clr.function",ns,inherits=FALSE)
  stopifnot(identical(body(original),body(installed),ignore.srcref=TRUE),identical(formals(original),formals(installed)))
  # Keep the package's RNG primitive; stream scheduling does not substitute one.
  get("rdirichlet",ns,inherits=FALSE)
}
al_load <- function(id) {
  entry <- al_entry(id)
  dataset <- substr(id,1L,1L);year <- substring(id,3L)
  counts <- readRDS(file.path(al_root,"15_FINAL_AUDIT/PHASE5_NATIVE_COUNT_INPUT_20261009_V1",paste0(dataset,"_COUNTS.rds")))
  mapping <- read.delim(file.path(al_root,"15_FINAL_AUDIT/PHASE5_UNIT_MAP_20261009_V1/SAMPLE_UNIT_MAP.tsv"),check.names=FALSE,stringsAsFactors=FALSE,colClasses="character")
  weights <- read.delim(file.path(al_root,"15_FINAL_AUDIT/PHASE5_UNIT_WEIGHTS_20261009_V1/SAMPLE_CONTRAST_WEIGHTS.tsv"),check.names=FALSE,stringsAsFactors=FALSE,colClasses="character")
  meta <- mapping[mapping$dataset==dataset & mapping$year==year,,drop=FALSE]
  stopifnot(nrow(meta)==entry$samples,!anyDuplicated(meta$run),all(meta$run %in% colnames(counts)))
  samples <- colnames(counts)[colnames(counts) %in% meta$run]
  meta <- meta[match(samples,meta$run),,drop=FALSE]
  w <- weights[match(samples,weights$run),,drop=FALSE]
  stopifnot(identical(w$run,samples),all(w$dataset==dataset),all(w$year==year))
  full <- counts[,samples,drop=FALSE]
  stopifnot(nrow(full)==entry$universe,ncol(full)==entry$samples,all(is.finite(full)),all(full>=0),all(full==floor(full)),all(full<=2^53))
  keep <- rowSums(full)>0
  stopifnot(sum(keep)==entry$clr_features,sum(!keep)==entry$cohort_total_zero)
  coefficient <- as.numeric(w$numerator)/as.numeric(w$denominator)
  stopifnot(all(is.finite(coefficient)),all(abs(coefficient)>0),abs(sum(coefficient[coefficient>0])-1)<1e-12,abs(sum(coefficient[coefficient<0])+1)<1e-12)
  pairs <- unique(meta$spatial_pair_year)
  stopifnot(length(pairs)==entry$pairs,length(unique(meta$allocation_unit_year))==entry$N_units)
  # Link sign to actual frozen metadata, not a newly inferred sample-name code.
  low <- if(dataset=="C")"Low_N" else "LN"
  high <- if(dataset=="C")"Std_N" else "HN"
  stopifnot(all(meta$nitrogen_treatment %in% c(low,high)),identical(coefficient>0,meta$nitrogen_treatment==low))
  stopifnot(!anyDuplicated(rownames(full)),identical(meta$run,samples))
  list(entry=entry,full=full,keep=keep,meta=meta,weight=w,coefficient=coefficient,pairs=pairs)
}
al_rng <- function(seed) {
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed)
}
al_stream <- function(full,conditions,weights,pair_labels,pairs,seed,path,checker=NULL) {
  stopifnot(!file.exists(path))
  keep <- rowSums(full)>0
  reads <- as.data.frame(full[keep,,drop=FALSE],check.names=FALSE)
  rn <- rownames(reads);stopifnot(length(rn)>=2L)
  subset <- get("aldex.set.mode",asNamespace("ALDEx2"))(reads,conditions,"all")
  stopifnot(is.numeric(subset),identical(as.integer(subset),seq_len(nrow(reads))))
  draw <- get("rdirichlet",asNamespace("ALDEx2"),inherits=FALSE)
  al_rng(seed);initial_rng <- .Random.seed
  aggregate <- matrix(0,nrow(reads),128L,dimnames=list(rn,NULL))
  pair_values <- lapply(pairs,function(x)aggregate);names(pair_values)<-pairs
  con <- file(path,"wb");on.exit(close(con),add=TRUE)
  for(i in seq_len(ncol(reads))) {
    # Identical native per-sample call, complete vector and call order.
    p <- t(draw(128L,reads[[i]]+0.5));rownames(p)<-rn
    z <- apply(log2(p),2L,function(col)col-mean(col[subset]))
    stopifnot(all(is.finite(p)),all(is.finite(z)),identical(dim(z),c(nrow(reads),128L)))
    if(!is.null(checker))checker(i,p,z)
    writeBin(as.double(z),con,size=8L,endian="little")
    aggregate <- aggregate + weights[i]*z
    name <- pair_labels[i]
    pair_values[[name]] <- pair_values[[name]] + length(pairs)*weights[i]*z
    if(i%%500L==0L)cat("STREAM_SAMPLE_COMPLETE",i,"OF",ncol(reads),"\n")
  }
  list(aggregate=aggregate,pairs=pair_values,initial_rng=initial_rng,final_rng=.Random.seed,
       features=rn,samples=colnames(reads),zero_rows=rownames(full)[!keep])
}
al_summary <- function(result,full,input) {
  entry <- input$entry
  summarize <- function(x)list(median_log2=unname(median(x)),conditional_dirichlet_low95_log2=unname(quantile(x,0.025,type=7)),
     conditional_dirichlet_high95_log2=unname(quantile(x,0.975,type=7)),conditional_dirichlet_sd_log2=sd(x))
  rows <- lapply(rownames(full),function(f) {
    observed <- f %in% result$features
    values <- if(observed)summarize(result$aggregate[f,]) else list(median_log2=NULL,conditional_dirichlet_low95_log2=NULL,conditional_dirichlet_high95_log2=NULL,conditional_dirichlet_sd_log2=NULL)
    supported_units <- unique(input$meta$allocation_unit_year[full[f,]>0])
    c(list(family=f,status=if(observed)"CONDITIONAL_DIRECTION_ONLY_NOT_INFERENTIAL" else "NOT_OBSERVED_IN_COHORT",
      cohort_observed_nonzero_samples=sum(full[f,]>0),cohort_observed_sample_count=ncol(full),
      units_with_observed_reads=as.list(supported_units),
      missing_unit_read_support=length(supported_units)<entry$N_units,
      direction=if(observed)sign(values$median_log2) else NULL,biological_SE=NULL,biological_CI=NULL,P=NULL,FDR=NULL,
      inference_missing_reason="NOT_ESTIMATED_IN_APPROVED_DIRECTIONAL_SENSITIVITY",
      pair_summaries=if(observed)lapply(result$pairs,function(x)summarize(x[f,])) else NULL),values)
  })
  list(cohort=entry$id,mc_samples=128L,units="log2 CLR difference",interval_scope="Conditional Dirichlet count uncertainty; not biological CI",
       N_units=entry$N_units,spatial_pairs=entry$pairs,samples=entry$samples,universe=entry$universe,rows=rows,
       ANCOM_annotation_join_status="PENDING_SEPARATE_ALL_FAMILY_REPORT_NOT_A_SELECTION_RULE",
       biological_inference_computed=FALSE,scientific_release=FALSE)
}
al_save_result <- function(result,full,input,folder,kind) {
  writeBin(as.double(full),file.path(folder,"INPUT_COUNTS.bin"),size=8L,endian="little")
  writeBin(as.double(result$aggregate),file.path(folder,"CONTRAST_DRAWS.bin"),size=8L,endian="little")
  writeBin(unlist(lapply(result$pairs,as.double),use.names=FALSE),file.path(folder,"PAIR_DRAWS.bin"),size=8L,endian="little")
  summary <- al_summary(result,full,input)
  summary$pair_unit_sample_counts <- lapply(names(result$pairs),function(p) {
    meta<-input$meta[input$meta$spatial_pair_year==p,,drop=FALSE]
    list(pair=p,units=lapply(unique(meta$allocation_unit_year),function(u) {
      m<-meta[meta$allocation_unit_year==u,,drop=FALSE]
      list(unit=u,nitrogen=unique(m$nitrogen_treatment),samples=nrow(m),observed_cells=length(unique(m$sampling_cell)))
    }))
  })
  al_json(summary,file.path(folder,"SUMMARY.json"))
  al_json(list(cohort=input$entry$id,kind=kind,features=result$features,universe=rownames(full),samples=colnames(full),
    pair_order=names(result$pairs),metadata=input$meta,coefficients=input$weight,
    shape=list(features=length(result$features),samples=ncol(full),mc=128L),
    binary_layout="CLR_STREAM.bin: sample-major blocks; each block R column-major feature by MC little-endian float64. CONTRAST same matrix; PAIR pair-major same matrices.",
    canonical_full_count_digest=digest::digest(full,algo="sha256"),zero_count=sum(full==0),
    original_full_count_digest=digest::digest(input$full,algo="sha256"),
    original_zero_count=sum(input$full==0),sample_or_family_filter_by_outcome=FALSE),file.path(folder,"INPUT_IDENTITY.json"))
  saveRDS(list(initial=result$initial_rng,final=result$final_rng),file.path(folder,"RNG_STATES.rds"),compress=FALSE)
  al_json(list(initial_sha256=digest::digest(result$initial_rng,algo="sha256"),final_sha256=digest::digest(result$final_rng,algo="sha256"),
               rng_kind=RNGkind()),file.path(folder,"RNG_IDENTITY.json"))
  writeLines(capture.output(sessionInfo()),file.path(folder,"SESSION_INFO.txt"),useBytes=TRUE)
}
