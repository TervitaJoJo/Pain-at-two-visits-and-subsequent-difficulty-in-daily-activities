suppressPackageStartupMessages(library(metafor))
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1],winslash='/',mustWork=TRUE)
newroot <- dirname(root)
oldroot <- file.path(dirname(newroot),'pain_transition_diabetes_adl_20261029')
counts <- read.csv(file.path(newroot,'results/four_cohort_primary_comparison.csv'))
cohorts <- c('HRS','ELSA','SHARE','CHARLS')
inputs <- lapply(cohorts,function(cc) {
 fn <- if(cc %in% c('HRS','ELSA')) file.path(oldroot,'results',paste0(cc,'_contrasts.csv')) else file.path(newroot,'results',paste0(cc,'_estimates.csv'))
 a <- read.csv(fn)
 a <- a[a$variant=='primary_MI' & a$contrast=='11_vs_00',]
 stopifnot(nrow(a)==1, a$imputations==50, a$se>0)
 z <- counts[counts$cohort==cc,]
 stopifnot(abs(a$estimate-z$aHR)<1e-10,a$persons==z$with_followup)
 data.frame(cohort=cc,logHR=a$log_estimate,se=a$se,HR=a$estimate,lower=a$lower,upper=a$upper,cohort_p=a$p,MI_df=a$df,eligible=z$eligible,persons=z$with_followup,events=z$events,source=sub(paste0(normalizePath(getwd(),winslash='/'),'/'),'',fn,fixed=TRUE))
})
d <- do.call(rbind,inputs)
stopifnot(!anyDuplicated(d$cohort),all(abs(exp(d$logHR)-d$HR)<1e-10))
run_model <- function(dat,label,method='REML',test='adhoc') {
 fit <- rma.uni(yi=dat$logHR,sei=dat$se,method=method,test=test,slab=dat$cohort,control=list(threshold=1e-12,maxiter=1000))
 dfmean <- if(test %in% c('knha','adhoc','t')) fit$k-1 else Inf
 tq <- if(is.finite(dfmean)) qt(.975,dfmean) else qnorm(.975)
 pi_df <- fit$k-2
 pi_margin <- if(pi_df>0) qt(.975,pi_df)*sqrt(fit$tau2+as.numeric(fit$se)^2) else NA_real_
 out <- data.frame(analysis=label,k=fit$k,persons=sum(dat$persons),events=sum(dat$events),method=method,test=test,logHR=as.numeric(coef(fit)),se=as.numeric(fit$se),HR=exp(as.numeric(coef(fit))),lower=exp(fit$ci.lb),upper=exp(fit$ci.ub),p=fit$pval,df=dfmean,tau2=fit$tau2,I2=fit$I2,Q=fit$QE,Q_df=fit$k-1,Q_p=fit$QEp,HK_scale=if(is.null(fit$s2w)) NA_real_ else fit$s2w,pred_lower=exp(as.numeric(coef(fit))-pi_margin),pred_upper=exp(as.numeric(coef(fit))+pi_margin),pred_df=pi_df)
 list(fit=fit,row=out)
}
models <- list()
models[['all_four_modified_HK']] <- run_model(d,'all_four_modified_HK')
models[['all_four_HK']] <- run_model(d,'all_four_HK',test='knha')
models[['all_four_Wald']] <- run_model(d,'all_four_Wald',test='z')
models[['all_four_common_effect']] <- run_model(d,'all_four_common_effect',method='FE',test='z')
for(cc in cohorts) {
 lab <- paste0('omit_',cc)
 models[[lab]] <- run_model(d[d$cohort!=cc,],lab)
}
out <- do.call(rbind,lapply(models,function(m) m$row))
d$weight_percent <- as.numeric(weights(models[['all_four_modified_HK']]$fit))
write.csv(d,file.path(root,'results/meta_inputs.csv'),row.names=FALSE)
write.csv(out,file.path(root,'results/meta_estimates.csv'),row.names=FALSE)
write.csv(data.frame(cohort=d$cohort,weight_percent=d$weight_percent),file.path(root,'results/cohort_weights.csv'),row.names=FALSE)
saveRDS(lapply(models,function(m)m$fit),file.path(root,'results/meta_models.rds'))
writeLines(capture.output(sessionInfo()),file.path(root,'audit/R_sessionInfo.txt'))
ci <- tryCatch(confint(models[['all_four_modified_HK']]$fit),error=function(e)conditionMessage(e))
writeLines(capture.output(print(ci)),file.path(root,'results/heterogeneity_uncertainty.txt'))
writeLines(capture.output(print(models[['all_four_modified_HK']]$fit)),file.path(root,'results/primary_model.txt'))
print(out[,c('analysis','k','HR','lower','upper','p','tau2','I2','Q_p','pred_lower','pred_upper')],row.names=FALSE,digits=7)
