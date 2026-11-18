# New observational adapter. Installed scientific functions are never modified.
instrument_observer <- function(original, callback_name, private_env) {
  old <- as.list(body(original))
  stopifnot(identical(old[[1]], as.name("{")),
            is.call(old[[length(old)]]),
            identical(old[[length(old)]][[1]], as.name("return")))
  observer <- substitute(CALL(environment()), list(CALL=as.name(callback_name)))
  revised <- append(old, list(observer), after=length(old)-1L)
  clone <- original
  body(clone) <- as.call(revised)
  environment(clone) <- private_env
  stopifnot(identical(as.list(body(clone))[-(length(old))], old),
            identical(formals(clone), formals(original)))
  clone
}

write_new_json <- function(x,path) {
  stopifnot(!file.exists(path))
  jsonlite::write_json(x,path,auto_unbox=TRUE,pretty=TRUE,na="null",null="null",digits=NA)
}
write_new_table <- function(x,path) {
  stopifnot(!file.exists(path))
  write.table(x,path,sep="\t",row.names=FALSE,quote=TRUE,na="NA",fileEncoding="UTF-8")
}

make_observed_ancom <- function(directory) {
  ns <- asNamespace("ANCOMBC")
  private <- new.env(parent=ns)
  seen <- new.env(parent=emptyenv()); seen$remle <- 0L; seen$top <- 0L
  installed_top <- get("ancombc2",ns); installed_remle <- get(".iter_remle",ns)
  originals <- list(top_body=body(installed_top),remle_body=body(installed_remle))
  private$phase5_observe_remle <- function(e) {
    seen$remle <- seen$remle+1L
    # Initial call computes theta; final call receives it. The frozen helper
    # does not expose its final convergence state in the public return object.
    stage <- if (exists("epsilon",envir=e,inherits=FALSE)) "initial" else "final"
    z <- get("output",e); expected <- get("fix_eff",e)
    fits <- z$fits
    diagnostics <- lapply(seq_along(fits),function(i) {
      fit <- fits[[i]]; ok <- inherits(fit,"lmerModLmerTest")
      if (!ok) return(list(family=rownames(z$beta)[i],fit_returned=FALSE,
                           qualified_numerical_fit=FALSE,reason="FITTING_RETURNED_NA"))
      sm <- summary(fit); vc <- as.data.frame(lme4::VarCorr(fit))
      fixed <- sm$coefficients
      drop <- attr(lme4::getME(fit,"X"),"col.dropped")
      msgs <- fit@optinfo$conv$lme4$messages
      opt <- fit@optinfo$conv$opt
      singular <- lme4::isSingular(fit,tol=1e-4)
      complete <- identical(rownames(fixed),expected)
      target <- if ("NLN" %in% rownames(fixed)) unname(fixed["NLN",c("Estimate","Std. Error","df")]) else rep(NA_real_,3)
      valid <- complete && length(drop)==0 && all(is.finite(target)) && target[2]>0 && target[3]>0 &&
        !singular && length(msgs)==0 && (length(opt)==0 || all(opt==0))
      list(family=rownames(z$beta)[i],fit_returned=TRUE,nobs=stats::nobs(fit),
           coefficient_names=rownames(fixed),N_raw_estimate=target[1],N_raw_standard_error=target[2],N_fit_df=target[3],
           fixed_basis_exact=complete,dropped_columns=drop,singular=singular,
           optimizer_code=opt,convergence_messages=msgs,variance_components=vc,
           qualified_numerical_fit=valid)
    })
    write_new_json(list(stage=stage,call=seen$remle,diagnostics=diagnostics,
                        iteration_count=if(stage=="initial") get("iterNum",e) else NULL,
                        final_epsilon=if(stage=="initial") get("epsilon",e) else NULL,
                        tolerance=get("tol",e),maximum_iterations=get("max_iter",e)),
                   file.path(directory,paste0(stage,"_fit_diagnostics.json")))
    invisible(NULL)
  }
  private$phase5_observe_top <- function(e) {
    seen$top <- seen$top+1L
    out <- get("out",e); para <- get("para2",e)
    expected <- get("fix_eff",e); taxa <- rownames(para$beta)
    raw_v <- para$var_hat; d <- get("var_delta",e); stabilizer <- get("s02",e)
    res <- out$res
    stopifnot(identical(as.character(res$taxon),taxa))
    rows <- lapply(seq_along(expected),function(k) {
      nm <- expected[k]
      data.frame(family=taxa,coefficient=nm,raw_beta=para$beta[,k],raw_variance=raw_v[,k],
        bias_variance=d[k],stabilizer=if(length(stabilizer)==1) stabilizer else stabilizer[k],
        dof=para$dof[,k],lfc=res[[paste0("lfc_",nm)]],se=res[[paste0("se_",nm)]],
        W=res[[paste0("W_",nm)]],p=res[[paste0("p_",nm)]],q=res[[paste0("q_",nm)]],
        passed_native_zero_sensitivity=res[[paste0("passed_ss_",nm)]],check.names=FALSE)
    })
    tab <- do.call(rbind,rows)
    tab$CI95_low <- tab$lfc-stats::qt(0.975,tab$dof)*tab$se
    tab$CI95_high <- tab$lfc+stats::qt(0.975,tab$dof)*tab$se
    write_new_table(tab,file.path(directory,"ALL_COEFFICIENT_COMPONENTS.tsv"))
    write_new_table(tab[tab$coefficient=="NLN",,drop=FALSE],file.path(directory,"N_ASSOCIATIONS_RAW.tsv"))
    write_new_json(list(shared_families=rownames(out$feature_table),samples=colnames(out$feature_table),
                        fixed_coefficients=expected,sampling_fractions=out$samp_frac,
                        bias_variance=d,stabilizers=stabilizer,full_global_table_retained=TRUE,
                        iteration_observer_calls=seen$remle,top_observer_calls=seen$top),
                   file.path(directory,"GLOBAL_COMPONENTS.json"))
    invisible(NULL)
  }
  private$.iter_remle <- instrument_observer(installed_remle,"phase5_observe_remle",private)
  fn <- instrument_observer(installed_top,"phase5_observe_top",private)
  verify <- function() {
    stopifnot(identical(body(get("ancombc2",ns)),originals$top_body),
              identical(body(get(".iter_remle",ns)),originals$remle_body),seen$remle==2L,seen$top==1L)
    TRUE
  }
  list(run=fn,verify=verify)
}
