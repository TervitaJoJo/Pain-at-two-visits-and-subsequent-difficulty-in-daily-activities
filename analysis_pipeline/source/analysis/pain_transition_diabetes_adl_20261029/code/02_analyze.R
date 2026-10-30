Sys.setenv(LANGUAGE='en')
suppressPackageStartupMessages({library(mice); library(sandwich); library(splines)})
args <- commandArgs(trailingOnly=TRUE)
cohort <- if(length(args)) args[1] else 'HRS'
stopifnot(cohort %in% c('HRS','ELSA'))
root <- normalizePath('analysis/pain_transition_diabetes_adl_20261029',winslash='/',mustWork=TRUE)
out <- file.path(root,'results'); audit <- file.path(root,'audit')
seed <- if(cohort=='HRS') 202610291 else 202610292
covars <- c('age','female','education','smoking','arthritis','cvd','cesd')
rhs <- 'exposure + age + female + education + smoking + arthritis + cvd + cesd + entry + interval'
results <- list(); diagnostics <- list(); weightdiag <- list(); impdiag <- list(); modelwarnings <- list(); mainfits <- list()
load_data <- function(sample) {
 d <- read.csv(file.path(root,'data_private',paste0(cohort,'_',sample,'.csv')),colClasses=if(sample=='longitudinal') c(exposure='character') else NA)
 d$education <- factor(d$education,levels=1:3); d$entry <- factor(d$entry,levels=c(2010,2012))
 for(v in c('female','smoking','arthritis','cvd')) d[[v]] <- factor(d[[v]],levels=c(0,1))
 if(sample=='longitudinal') d$exposure <- factor(d$exposure,levels=c('00','01','10','11'))
 return(d)
}
intervals <- function(d,composite=FALSE,risk=FALSE) {
 z <- d
 if(composite) {
  for(j in 1:2) z[[paste0('y',j)]] <- ifelse(z[[paste0('dead',j)]]==1,1,z[[paste0('y',j)]])
 }
 a <- z; a$interval <- factor(1,levels=1:2); a$y <- z$y1; a$years <- if(composite) 2 else z$t1
 keep2 <- !is.na(z$y1) & z$y1==0
 b <- z[keep2,,drop=FALSE]; b$interval <- factor(2,levels=1:2); b$y <- b$y2; b$years <- if(composite) 2 else b$t2-b$t1
 q <- rbind(a,b); q$observed <- as.integer(!is.na(q$y)); q <- q[order(q$pid,q$interval),]
 if(!risk) q <- q[q$observed==1,]
 stopifnot(all(is.finite(q$years[q$observed==1])),all(q$years[q$observed==1]>0))
 rownames(q) <- NULL
 q
}
impute_data <- function(d,sample) {
 if(sample=='longitudinal') {
  obs <- intervals(d)
  d$observed_event <- as.integer(d$pid %in% obs$pid[obs$y==1])
  cnt <- table(obs$pid); yrs <- tapply(obs$years,obs$pid,sum)
  d$n_intervals <- as.numeric(cnt[as.character(d$pid)]); d$n_intervals[is.na(d$n_intervals)] <- 0
  d$followup_years <- as.numeric(yrs[as.character(d$pid)]); d$followup_years[is.na(d$followup_years)] <- 0
  cols <- c('exposure',covars,'entry','function3','observed_event','n_intervals','followup_years')
 } else cols <- c('pain',covars,'entry','function3','adl')
 dat <- d[,cols]; met <- setNames(rep('',length(cols)),cols)
 for(v in cols) if(anyNA(dat[[v]])) met[v] <- if(v=='education') 'polyreg' else if(v %in% c('female','smoking','arthritis','cvd')) 'logreg' else 'pmm'
 pred <- matrix(1,length(cols),length(cols),dimnames=list(cols,cols)); diag(pred)<-0
 imp <- mice(dat,m=50,maxit=20,method=met,predictorMatrix=pred,seed=seed+if(sample=='cross_sectional') 100L else 0L,printFlag=FALSE,donors=5)
 saveRDS(imp,file.path(root,'data_private',paste0(cohort,'_',sample,'_mice50.rds')))
 if(!is.null(imp$loggedEvents)) write.csv(imp$loggedEvents,file.path(audit,paste0(cohort,'_',sample,'_mice_events.csv')),row.names=FALSE)
 cv <- mice::convergence(imp)
 write.csv(cv,file.path(audit,paste0(cohort,'_',sample,'_mice_convergence.csv')),row.names=FALSE)
 pdf(file.path(root,'figures',paste0(cohort,'_',sample,'_imputation_chains.pdf')),width=9,height=6)
 print(plot(imp)); dev.off()
 for(v in names(met)[met!='']) {
  vals <- unlist(lapply(1:50,function(i) as.numeric(as.character(complete(imp,i)[[v]][is.na(dat[[v]])]))))
  impdiag[[length(impdiag)+1]] <<- data.frame(cohort=cohort,sample=sample,variable=v,method=met[v],missing=sum(is.na(dat[[v]])),imputed_mean=mean(vals),imputed_sd=sd(vals),imputed_min=min(vals),imputed_max=max(vals),logged_events=if(is.null(imp$loggedEvents)) 0 else nrow(imp$loggedEvents))
 }
 completed <- lapply(1:50,function(i) {
  dd <- d; cc <- complete(imp,i)
  for(v in cols) {
   known <- !is.na(dat[[v]])
   stopifnot(all(as.character(cc[[v]][known])==as.character(dat[[v]][known])))
   dd[[v]] <- cc[[v]]
  }
  stopifnot(!anyNA(dd[,covars])); dd
 })
 completed
}
fit_glm <- function(d,formula,family,variant,i,weights=NULL) {
 if(is.null(weights)) weights <- rep(1,nrow(d))
 d$.weight <- weights
 warns <- character()
 fit <- withCallingHandlers(glm(formula,data=d,family=family,weights=.weight,control=glm.control(maxit=100)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
 if(length(warns)) modelwarnings[[length(modelwarnings)+1]] <<- data.frame(cohort=cohort,variant=variant,imputation=i,message=unique(warns))
 stopifnot(fit$converged,all(is.finite(coef(fit))),fit$rank==length(coef(fit)))
 V <- vcovCL(fit,cluster=d$pid,type='HC1',cadjust=TRUE)
 stopifnot(all(is.finite(V)),all(diag(V)>0))
 X <- model.matrix(fit)
 diagnostics[[length(diagnostics)+1]] <<- data.frame(cohort=cohort,variant=variant,imputation=i,n=nrow(d),persons=length(unique(d$pid)),events=sum(d$y),parameters=length(coef(fit)),converged=fit$converged,iterations=fit$iter,min_fitted=min(fitted(fit)),max_fitted=max(fitted(fit)),max_leverage=max(hatvalues(fit)),max_cook=max(cooks.distance(fit)),min_cov_eigen=min(eigen(V,symmetric=TRUE,only.values=TRUE)$values),kappa_scaled=if(ncol(X)>1) kappa(scale(X[,-1,drop=FALSE])) else NA)
 list(b=coef(fit),V=V,n=length(unique(d$pid)),k=length(coef(fit)),fit=if(i==1) fit else NULL)
}
contrast_vector <- function(nms,type,interaction=FALSE) {
 v<-setNames(rep(0,length(nms)),nms)
 if(type=='pain_vs_no_pain') v['pain']<-1
 else if(type=='10_vs_11_severity') v['exposure10']<-1
 else {
  groups<-strsplit(type,'_vs_')[[1]]
  for(j in 1:2) if(groups[j]!='00') {
   key<-paste0('exposure',groups[j]);v[key]<-if(j==1) 1 else -1
   if(interaction) {
    key2<-paste0(key,':interval2'); if(!key2 %in% nms) key2<-paste0('interval2:',key)
    stopifnot(key2 %in% nms);v[key2]<-if(j==1) 1 else -1
   }
  }
 }
 stopifnot(identical(names(v),nms));v
}
pool_contrasts <- function(fits,variant,types=c('11_vs_00','01_vs_00','10_vs_00','10_vs_11'),second_interval=FALSE) {
 rows<-lapply(types,function(type) {
  v<-contrast_vector(names(fits[[1]]$b),type,second_interval)
  Q<-vapply(fits,function(f) sum(v*f$b),numeric(1)); U<-vapply(fits,function(f) as.numeric(t(v)%*%f$V%*%v),numeric(1))
  n<-min(vapply(fits,function(f) f$n,numeric(1))); k<-max(vapply(fits,function(f) f$k,numeric(1)))
  if(length(fits)>1) {ps<-pool.scalar(Q,U,n=n,k=k);est<-ps$qbar;se<-sqrt(ps$t);df<-ps$df;fmi<-ps$fmi;mcse<-sqrt(ps$b/length(fits))}
  else {est<-Q;se<-sqrt(U);df<-n-k;fmi<-0;mcse<-0}
  p<-2*pt(abs(est/se),df=df,lower.tail=FALSE);crit<-qt(.975,df=df)
  data.frame(cohort=cohort,variant=variant,contrast=type,imputations=length(fits),persons=n,log_estimate=est,se=se,estimate=exp(est),lower=exp(est-crit*se),upper=exp(est+crit*se),p=p,df=df,fmi=fmi,mcse_log=mcse,p_holm_secondary=NA_real_)
 })
 ans<-do.call(rbind,rows)
 if(variant=='primary_MI') {ix<-ans$contrast!='11_vs_00';ans$p_holm_secondary[ix]<-p.adjust(ans$p[ix],method='holm')}
 results[[length(results)+1]] <<- ans
 ans
}
observation_weights <- function(d,i) {
 r<-intervals(d,risk=TRUE)
 den<-glm(as.formula(paste('observed ~',rhs)),data=r,family=binomial())
 num<-glm(observed~exposure+entry+interval,data=r,family=binomial())
 stopifnot(den$converged,num$converged)
 pd<-predict(den,type='response');pn<-predict(num,type='response');stopifnot(all(pd>0),all(pd<1))
 w<-ave(pn/pd,r$pid,FUN=cumprod);obs<-r$observed==1
 q<-quantile(w[obs],c(.01,.99));trunc<-pmin(pmax(w[obs],q[1]),q[2])
 weightdiag[[length(weightdiag)+1]] <<- data.frame(cohort=cohort,imputation=i,risk_intervals=nrow(r),observed_intervals=sum(obs),min_den=min(pd),max_den=max(pd),raw_min=min(w[obs]),raw_max=max(w[obs]),q01=q[1],q99=q[2],weight_mean=mean(trunc),weight_sd=sd(trunc),ess=sum(trunc)^2/sum(trunc^2),n_truncated=sum(w[obs]<q[1]|w[obs]>q[2]))
 list(d=r[obs,],w=trunc)
}
cat(cohort,'building imputations\n');flush.console()
d<-load_data('longitudinal'); completed<-impute_data(d,'longitudinal')
form<-as.formula(paste('y ~',rhs,'+ offset(log(years))'))
fits<-lapply(seq_along(completed),function(i) fit_glm(intervals(completed[[i]]),form,binomial('cloglog'),'primary_MI',i))
mainfits<-fits
print(pool_contrasts(fits,'primary_MI')[,c('cohort','contrast','estimate','lower','upper','p')]);flush.console()
cc<-d[complete.cases(d[,covars]),]
pool_contrasts(list(fit_glm(intervals(cc),form,binomial('cloglog'),'complete_case',1)),'complete_case')
cat(cohort,'running preplanned sensitivity analyses\n');flush.console()
fits<-lapply(seq_along(completed),function(i) fit_glm(intervals(completed[[i]]),update(form,.~.+function3),binomial('cloglog'),'function_adjusted_MI',i))
pool_contrasts(fits,'function_adjusted_MI')
fits<-lapply(seq_along(completed),function(i) {
 dd<-completed[[i]];dd<-dd[dd$exposure %in% c('10','11') & dd$pain_severity %in% 1:3,];dd$exposure<-factor(dd$exposure,levels=c('11','10'));dd$pain_severity<-factor(dd$pain_severity,levels=1:3)
 fit_glm(intervals(dd),update(form,.~.+pain_severity),binomial('cloglog'),'severity_adjusted_MI',i)
})
pool_contrasts(fits,'severity_adjusted_MI','10_vs_11_severity')
fits<-lapply(seq_along(completed),function(i) {
 ww<-observation_weights(completed[[i]],i)
 fit_glm(ww$d,form,quasibinomial('cloglog'),'observation_weighted_MI',i,ww$w)
})
pool_contrasts(fits,'observation_weighted_MI')
if(cohort=='HRS') {
 fits<-lapply(seq_along(completed),function(i) fit_glm(intervals(completed[[i]],composite=TRUE),form,binomial('cloglog'),'ADL_or_death_nominal_MI',i))
 pool_contrasts(fits,'ADL_or_death_nominal_MI')
}
fits<-lapply(seq_along(completed),function(i) fit_glm(intervals(completed[[i]]),update(form,.~.+exposure:interval),binomial('cloglog'),'interval_interaction_MI',i))
pool_contrasts(fits,'interval1_diagnostic_MI')
pool_contrasts(fits,'interval2_diagnostic_MI',second_interval=TRUE)
fits<-lapply(seq_along(completed),function(i) fit_glm(intervals(completed[[i]]),update(form,.~.-age-cesd+ns(age,df=3)+ns(cesd,df=3)),binomial('cloglog'),'nonlinear_covariates_MI',i))
pool_contrasts(fits,'nonlinear_covariates_MI')
cat(cohort,'cross-sectional supplement\n');flush.console()
cs<-load_data('cross_sectional');cscomp<-impute_data(cs,'cross_sectional')
csform<-y ~ pain + age + female + education + smoking + arthritis + cvd + cesd + entry
fits<-lapply(seq_along(cscomp),function(i) {dd<-cscomp[[i]];dd$y<-dd$adl;fit_glm(dd,csform,poisson('log'),'cross_sectional_MI',i)})
pool_contrasts(fits,'cross_sectional_MI','pain_vs_no_pain')
cc<-cs[complete.cases(cs[,covars]),];cc$y<-cc$adl
pool_contrasts(list(fit_glm(cc,csform,poisson('log'),'cross_sectional_CC',1)),'cross_sectional_CC','pain_vs_no_pain')
write.csv(do.call(rbind,results),file.path(out,paste0(cohort,'_contrasts.csv')),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(audit,paste0(cohort,'_model_diagnostics.csv')),row.names=FALSE)
write.csv(do.call(rbind,weightdiag),file.path(audit,paste0(cohort,'_weight_diagnostics.csv')),row.names=FALSE)
write.csv(do.call(rbind,impdiag),file.path(audit,paste0(cohort,'_imputation_summary.csv')),row.names=FALSE)
if(length(modelwarnings)) write.csv(do.call(rbind,modelwarnings),file.path(audit,paste0(cohort,'_model_warnings.csv')),row.names=FALSE)
saveRDS(list(primary=mainfits,results=results),file.path(root,'data_private',paste0(cohort,'_models.rds')))
writeLines(capture.output(sessionInfo()),file.path(audit,paste0(cohort,'_R_session.txt')))
writeLines(paste('Completed',cohort,format(Sys.time(),'%Y-%m-%d %H:%M:%S %z')),file.path(audit,paste0(cohort,'_complete.txt')))
cat(cohort,'ALL MODELS COMPLETED\n')
