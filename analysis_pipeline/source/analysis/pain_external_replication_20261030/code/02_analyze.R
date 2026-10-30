Sys.setenv(LANGUAGE='en')
suppressPackageStartupMessages({library(mice);library(sandwich)})
cohort<-commandArgs(trailingOnly=TRUE)[1];stopifnot(cohort %in% c('CHARLS','SHARE'))
root<-normalizePath('analysis/pain_external_replication_20261030',winslash='/',mustWork=TRUE)
filename<-if(cohort=='SHARE')'SHARE_primary_longitudinal.csv' else 'CHARLS_longitudinal.csv'
d<-read.csv(file.path(root,'data_private',filename),colClasses=c(exposure='character'))
d$exposure<-factor(d$exposure,levels=c('00','01','10','11'))
d$education<-factor(d$education,levels=1:3)
for(v in c('female','smoking','arthritis','cvd'))d[[v]]<-factor(d[[v]],levels=0:1)
if(cohort=='SHARE')d$country<-factor(d$country)
covars<-c('age','female','education','smoking','arthritis','cvd','cesd')
intervals<-function(z){
 a<-z;a$interval<-factor(1,levels=1:2);a$y<-a$y1;a$years<-a$t1
 b<-z[!is.na(z$y1)&z$y1==0,,drop=FALSE];b$interval<-factor(2,levels=1:2);b$y<-b$y2;b$years<-b$t2-b$t1
 q<-rbind(a,b);q<-q[!is.na(q$y),];q<-q[order(q$pid,q$interval),];rownames(q)<-NULL
 stopifnot(all(is.finite(q$years)),all(q$years>0),!anyDuplicated(paste(q$pid,q$interval)))
 q
}
obs<-intervals(d)
d$observed_event<-as.integer(d$pid %in% obs$pid[obs$y==1])
cnt<-table(obs$pid);yrs<-tapply(obs$years,obs$pid,sum)
d$n_intervals<-as.numeric(cnt[as.character(d$pid)]);d$n_intervals[is.na(d$n_intervals)]<-0
d$followup_years<-as.numeric(yrs[as.character(d$pid)]);d$followup_years[is.na(d$followup_years)]<-0
cols<-c('exposure',covars,if(cohort=='SHARE')c('country','function3'),'observed_event','n_intervals','followup_years')
dat<-d[,cols];met<-setNames(rep('',length(cols)),cols)
for(v in cols)if(anyNA(dat[[v]]))met[v]<-if(v=='education')'polyreg' else if(v %in% c('female','smoking','arthritis','cvd'))'logreg' else 'pmm'
pred<-matrix(1,length(cols),length(cols),dimnames=list(cols,cols));diag(pred)<-0
seed<-if(cohort=='CHARLS')202610301L else 202610302L
cat(cohort,'starting MI 50 x 20 at',format(Sys.time()),'\n');flush.console()
imp<-mice(dat,m=50,maxit=20,method=met,predictorMatrix=pred,seed=seed,donors=5,printFlag=FALSE)
saveRDS(imp,file.path(root,'data_private',paste0(cohort,'_mice50.rds')))
if(!is.null(imp$loggedEvents))write.csv(imp$loggedEvents,file.path(root,'audit',paste0(cohort,'_mice_logged_events.csv')),row.names=FALSE)
write.csv(mice::convergence(imp),file.path(root,'audit',paste0(cohort,'_mice_convergence.csv')),row.names=FALSE)
completed<-lapply(1:50,function(i){dd<-d;cc<-complete(imp,i);for(v in cols){known<-!is.na(dat[[v]]);stopifnot(all(as.character(dat[[v]][known])==as.character(cc[[v]][known])));dd[[v]]<-cc[[v]]};stopifnot(!anyNA(dd[,covars]));dd})
rhs<-paste(c('exposure',covars,'interval',if(cohort=='SHARE')'country'),collapse=' + ')
form<-as.formula(paste('y ~',rhs,'+ offset(log(years))'))
diags<-list();warns<-list();allresults<-list()
fitone<-function(dd,ff,variant,i){
 q<-intervals(dd);ww<-character()
 f<-withCallingHandlers(glm(ff,data=q,family=binomial('cloglog'),control=glm.control(maxit=100)),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')})
 V<-vcovCL(f,cluster=q$pid,type='HC1',cadjust=TRUE)
 stopifnot(f$converged,all(is.finite(coef(f))),f$rank==length(coef(f)),all(is.finite(V)),all(diag(V)>0))
 n<-length(unique(q$pid));k<-length(coef(f))
 diags[[length(diags)+1]]<<-data.frame(cohort=cohort,variant=variant,imputation=i,persons=n,intervals=nrow(q),events=sum(q$y),parameters=k,converged=f$converged,iterations=f$iter,min_fitted=min(fitted(f)),max_fitted=max(fitted(f)),max_leverage=max(hatvalues(f)),min_cov_eigen=min(eigen(V,symmetric=TRUE,only.values=TRUE)$values),warnings=length(ww))
 if(length(ww))warns[[length(warns)+1]]<<-data.frame(variant=variant,imputation=i,message=unique(ww))
 if(variant=='primary_MI'&&i==1){write.csv(q,file.path(root,'data_private',paste0(cohort,'_first_imputation_intervals.csv')),row.names=FALSE);write.csv(data.frame(term=names(coef(f)),estimate=coef(f),se=sqrt(diag(V))),file.path(root,'results',paste0(cohort,'_first_imputation_coefficients.csv')),row.names=FALSE)}
 list(b=coef(f),V=V,n=n,k=k,events=sum(q$y),intervals=nrow(q))
}
pool<-function(fits,variant){
 rr<-lapply(c('11_vs_00','01_vs_00','10_vs_00','10_vs_11'),function(type){
  v<-setNames(rep(0,length(fits[[1]]$b)),names(fits[[1]]$b));g<-strsplit(type,'_vs_')[[1]]
  for(j in 1:2)if(g[j]!='00')v[paste0('exposure',g[j])]<-if(j==1)1 else -1
  Q<-vapply(fits,function(f)sum(v*f$b),numeric(1));U<-vapply(fits,function(f)as.numeric(t(v)%*%f$V%*%v),numeric(1))
  n<-min(vapply(fits,function(f)f$n,numeric(1)));k<-max(vapply(fits,function(f)f$k,numeric(1)))
  if(length(fits)>1){ps<-pool.scalar(Q,U,n=n,k=k);est<-ps$qbar;se<-sqrt(ps$t);df<-ps$df;fmi<-ps$fmi;mcse<-sqrt(ps$b/length(fits))}else{est<-Q;se<-sqrt(U);df<-n-k;fmi<-0;mcse<-0}
  data.frame(cohort=cohort,variant=variant,contrast=type,imputations=length(fits),persons=n,events=fits[[1]]$events,intervals=fits[[1]]$intervals,estimate=exp(est),lower=exp(est-qt(.975,df)*se),upper=exp(est+qt(.975,df)*se),p=2*pt(abs(est/se),df,lower.tail=FALSE),log_estimate=est,se=se,df=df,fmi=fmi,mcse_log=mcse,p_holm_secondary=NA_real_)
 })
 rr<-do.call(rbind,rr);ix<-rr$contrast!='11_vs_00';rr$p_holm_secondary[ix]<-p.adjust(rr$p[ix],method='holm')
 allresults[[length(allresults)+1]]<<-rr;print(rr[,c('variant','contrast','persons','events','estimate','lower','upper','p')]);flush.console()
}
cat(cohort,'fitting primary models at',format(Sys.time()),'\n');flush.console()
fits<-lapply(seq_along(completed),function(i)fitone(completed[[i]],form,'primary_MI',i));saveRDS(fits,file.path(root,'data_private',paste0(cohort,'_primary_fits.rds')));pool(fits,'primary_MI')
cc<-d[complete.cases(d[,covars]),];pool(list(fitone(cc,form,'complete_case',1)),'complete_case')
if(cohort=='SHARE'){
 pool(lapply(seq_along(completed),function(i)fitone(completed[[i]],update(form,.~.+function3),'function_adjusted_MI',i)),'function_adjusted_MI')
 pool(lapply(seq_along(completed),function(i)fitone(completed[[i]][completed[[i]]$self_only==1,],form,'respondent_only_MI',i)),'respondent_only_MI')
}
write.csv(do.call(rbind,allresults),file.path(root,'results',paste0(cohort,'_estimates.csv')),row.names=FALSE)
write.csv(do.call(rbind,diags),file.path(root,'audit',paste0(cohort,'_model_diagnostics.csv')),row.names=FALSE)
if(length(warns))write.csv(do.call(rbind,warns),file.path(root,'audit',paste0(cohort,'_model_warnings.csv')),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'audit',paste0(cohort,'_R_session.txt')))
cat(cohort,'DONE at',format(Sys.time()),'\n')
