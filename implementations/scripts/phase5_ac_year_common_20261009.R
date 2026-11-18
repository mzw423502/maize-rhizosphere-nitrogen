# New A/C input adapter; installed ANCOM/lmer code remains untouched.
ac_root <- "."
ac_base <- file.path(ac_root,"15_FINAL_AUDIT/PHASE5_AC_YEAR_MODELS_20261009_V1")
source(file.path(ac_root,"scripts/phase5_b_random_model_common_20261009.R"))
ac_contract <- function() jsonlite::fromJSON(file.path(ac_root,
  "15_FINAL_AUDIT/PHASE5_AC_YEAR_MODEL_CONTRACT_20261009.json"),simplifyVector=FALSE)
ac_load <- function(id) {
  contract <- ac_contract()
  entries <- Filter(function(e)identical(e$id,id),contract$models)
  stopifnot(length(entries)==1L)
  entry <- entries[[1]]
  input <- jsonlite::fromJSON(file.path(ac_root,entry$model_input_file))
  meta <- input$metadata_in_source_count_order
  original <- readRDS(file.path(ac_root,input$native_count_source))
  stopifnot(nrow(original)==entry$original_family_universe,
    identical(rownames(original),input$families_in_order),
    identical(colnames(original)[meta$source_count_row_index_1based],meta$sample_id),
    !anyDuplicated(meta$sample_id),nrow(meta)==entry$samples)
  families <- unlist(entry$families_in_model_order,use.names=FALSE)
  excluded <- unlist(entry$not_estimable_families,use.names=FALSE)
  stopifnot(length(families)==entry$model_families,
    setequal(c(families,excluded),rownames(original)),length(intersect(families,excluded))==0L,
    identical(rownames(original)[rownames(original)%in%families],families))
  count <- original[families,meta$sample_id,drop=FALSE]
  stopifnot(all(is.finite(count)),all(count>=0),all(count==floor(count)),
            all(count<=2^53),all(rowSums(count)>0))
  rownames(meta) <- meta$sample_id
  meta$N <- factor(meta$N,levels=c("HN","LN"))
  for (key in intersect(c("spatial_pair","genotype","allocation_unit","split_plot_uid",
                         "split_plot_block_uid","subplot_uid","period","subplot_period_uid"),names(meta)))
    meta[[key]] <- factor(meta[[key]],levels=sort(unique(meta[[key]])))
  stopifnot(!anyNA(meta),identical(levels(meta$N),c("HN","LN")))
  list(entry=entry,input=input,count=count,metadata=meta,
       original_dimensions=dim(original),original_digest=digest::digest(original,algo="sha256"))
}
ac_canonical_digest <- function(count) {
  lines <- vapply(seq_len(ncol(count)),function(i)
    paste(c(colnames(count)[i],format(count[,i],scientific=FALSE,trim=TRUE)),collapse="\t"),character(1))
  digest::digest(paste0(paste(lines,collapse="\n"),"\n"),algo="sha256",serialize=FALSE)
}
