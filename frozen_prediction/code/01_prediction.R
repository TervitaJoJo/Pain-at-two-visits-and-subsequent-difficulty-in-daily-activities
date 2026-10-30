# Exploratory two-visit pain prediction. Never source the association MICE objects.
suppressPackageStartupMessages({library(splines);library(jsonlite)})
options(warn=1)
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2L)stop('Usage: Rscript 01_prediction.R INPUT_ANALYSIS_ROOT OUTPUT_ROOT [full|bootstrap]')
input_root<-normalizePath(args[1],winslash='/',mustWork=TRUE)
root<-normalizePath(args[2],winslash='/',mustWork=TRUE)
mode<-if(length(args)>2)args[3] else 'full'
seed<-20261030L
paths<-c(HRS='pain_transition_diabetes_adl_20261029/data_private/HRS_longitudinal.csv',ELSA='pain_transition_diabetes_adl_20261029/data_private/ELSA_longitudinal.csv',SHARE='pain_external_replication_20261030/data_private/SHARE_primary_longitudinal.csv')
features<-c('age','female','education','smoking','arthritis','cvd','depressed','function3')
allwarnings<-character();warning_records<-list();fit_count<-0L
quiet_glm<-function(form,data,family,weights=NULL){
 if(is.null(weights))weights<-rep(1,nrow(data));data$.w<-weights
 f<-withCallingHandlers(glm(form,data=data,family=family,weights=.w,control=glm.control(maxit=100)),warning=function(w){allwarnings<<-c(allwarnings,conditionMessage(w));warning_records[[length(warning_records)+1]]<<-data.frame(fit_number=fit_count+1L,formula=paste(deparse(form),collapse=' '),family=family$family,link=family$link,n=nrow(data),message=conditionMessage(w));invokeRestart('muffleWarning')})
 fit_count<<-fit_count+1L
 if(!f$converged || any(!is.finite(coef(f))))stop('Non-converged or rank-deficient GLM')
 f
}
load_data<-function(k){
 d<-read.csv(file.path(input_root,paths[k]),colClasses=c(exposure='character'))
 stopifnot(!anyDuplicated(d$pid),all(d$exposure %in% c('00','01','10','11')))
 d$gid<-seq_len(nrow(d));d$p1<-as.integer(substr(d$exposure,1,1));d$p2<-as.integer(substr(d$exposure,2,2));d$depressed<-as.integer(d$cesd>=if(k=='SHARE')4 else 3)
 d$education<-factor(d$education,levels=1:3)
 if(k=='SHARE')d$country<-factor(d$country) else d$entry<-factor(d$entry,levels=c(2010,2012))
 d$cohort<-k
 stopifnot(!anyNA(d[,c('p1','p2')]),all(na.omit(d$y1)%in%0:1),all(na.omit(d$y2)%in%0:1))
 d
}
fit_imputer<-function(d,raw=FALSE){
 fs<-if(raw)c(setdiff(features,'depressed'),'cesd') else features
 setNames(lapply(fs,function(v){x<-d[[v]];if(is.factor(x)||v%in%c('female','smoking','arthritis','cvd','depressed')){tb<-table(x);names(tb)[which.max(tb)]}else median(x,na.rm=TRUE)}),fs)
}
apply_imputer<-function(d,imp){
 for(v in names(imp)){miss<-is.na(d[[v]]);if(any(miss))d[[v]][miss]<-if(is.factor(d[[v]]))as.character(imp[[v]]) else as.numeric(imp[[v]])}
 stopifnot(!anyNA(d[,names(imp),drop=FALSE]));d
}
common_rhs<-"ns(age, knots=c(60,75), Boundary.knots=c(50,100)) + female + education + smoking + arthritis + cvd + depressed + function3"
model_formula<-function(model,raw=FALSE,cohort=NULL){
 rhs<-common_rhs
 if(raw)rhs<-paste(sub('depressed','cesd',rhs,fixed=TRUE),if(cohort=='SHARE')'+ country' else '+ entry')
 as.formula(paste('y ~ factor(interval) +',rhs,switch(model,C='',M0='+ p2',M1='+ p1 * p2')))
}
interval_data<-function(d){
 a<-d;a$y<-a$y1;a$interval<-1
 b<-d[!is.na(d$y1)&d$y1==0,,drop=FALSE];b$y<-b$y2;b$interval<-2
 z<-rbind(a,b);z[!is.na(z$y),,drop=FALSE]
}
fit_pipeline<-function(d,raw=FALSE,models=c('C','M0','M1')){
 imp<-fit_imputer(d,raw);dd<-apply_imputer(d,imp);q<-interval_data(dd)
 fits<-setNames(lapply(models,function(m)quiet_glm(model_formula(m,raw,unique(d$cohort)),q,binomial('cloglog'))),models)
 list(imputer=imp,fits=fits,raw=raw)
}
predict_pipeline<-function(fit,d){
 dd<-apply_imputer(d,fit$imputer);res<-array(NA_real_,c(nrow(d),length(fit$fits),2),dimnames=list(NULL,names(fit$fits),c('visit1','visit2')))
 for(m in names(fit$fits)){
  dd$interval<-1;h1<-predict(fit$fits[[m]],newdata=dd,type='response');dd$interval<-2;h2<-predict(fit$fits[[m]],newdata=dd,type='response')
  res[,m,1]<-h1;res[,m,2]<-1-(1-h1)*(1-h2)
 }
 stopifnot(all(is.finite(res)),all(res>=0&res<=1),all(res[,,2]>=res[,,1]-1e-12));res
}
make_folds<-function(d,seed,stratified=FALSE){
 set.seed(seed);ids<-unique(d$gid)
 if(stratified){
  z<-d[match(ids,d$gid),];s<-ifelse(is.na(z$y1),'missing1',ifelse(z$y1==1,'event1',ifelse(is.na(z$y2),'missing2',ifelse(z$y2==1,'event2','noevent'))))
  v<-integer(length(ids));for(g in unique(s)){ix<-which(s==g);v[ix]<-sample(rep(1:5,length.out=length(ix)))}
 }else v<-sample(rep(1:5,length.out=length(ids)))
 v[match(d$gid,ids)]
}
# Evaluation observation nuisance models are distinct from the frozen risk model.
observation_probs<-function(d,seed){
 folds<-make_folds(d,seed,FALSE);out<-matrix(NA_real_,nrow(d),2,dimnames=list(NULL,c('q1','q2')))
 form<-as.formula(paste('response ~',common_rhs,'+ p1 * p2'))
 for(k in 1:5){
  tr<-d[folds!=k,];te<-d[folds==k,];imp<-fit_imputer(tr);tr<-apply_imputer(tr,imp);te<-apply_imputer(te,imp)
  tr$response<-as.integer(!is.na(tr$y1));f1<-quiet_glm(form,tr,binomial());out[folds==k,1]<-predict(f1,te,type='response')
  tr<-tr[!is.na(tr$y1)&tr$y1==0,];tr$response<-as.integer(!is.na(tr$y2));f2<-quiet_glm(form,tr,binomial());out[folds==k,2]<-predict(f2,te,type='response')
 }
 stopifnot(all(is.finite(out)));out
}
endpoint<-function(d,h){
 y<-rep(NA_real_,nrow(d));i1<-!is.na(d$y1)
 if(h==1)y[i1]<-d$y1[i1] else {
  y[i1&d$y1==1]<-1
  ii<-i1&d$y1==0&!is.na(d$y2);y[ii]<-d$y2[ii]
 }
 y
}
get_weights<-function(d,pr,h,type='conditional'){
 if(type=='marginal'){
  q1<-mean(!is.na(d$y1));q2<-mean(!is.na(d$y2[!is.na(d$y1)&d$y1==0]));pr<-cbind(rep(q1,nrow(d)),rep(q2,nrow(d)))
 }
 q<-if(type=='untruncated')pmax(pr,1e-6) else pmin(pmax(pr,.02),.98)
 w<-1/q[,1]
 if(h==2){ii<-!is.na(d$y1)&d$y1==0;w[ii]<-1/(q[ii,1]*q[ii,2])}
 if(type=='conditional')w<-pmin(w,20)
 if(type=='unweighted')w[]<-1
 w[is.na(endpoint(d,h))]<-0;w
}
weighted_auc<-function(y,p,w){
 ix<-order(p);y<-y[ix];p<-p[ix];w<-w[ix]
 grp<-cumsum(c(TRUE,diff(p)!=0));cases<-as.vector(rowsum(w*(y==1),grp,reorder=FALSE));controls<-as.vector(rowsum(w*(y==0),grp,reorder=FALSE))
 sum(cases*(cumsum(controls)-.5*controls))/(sum(cases)*sum(controls))
}
metric<-function(d,p,pr,h,type='conditional',calibrate=TRUE){
 y<-endpoint(d,h);w<-get_weights(d,pr,h,type);ok<-!is.na(y)&w>0;n<-nrow(d)
 yy<-y[ok];pp<-p[ok];ww<-w[ok];bw<-sum(ww*(yy-pp)^2)
 val<-c(auc=weighted_auc(yy,pp,ww),brier=bw/n,brier_hajek=bw/sum(ww),observed=sum(ww*yy)/sum(ww),predicted=mean(p),weight_sum_ratio=sum(ww)/n,cal_intercept=NA_real_,cal_slope=NA_real_)
 if(calibrate){
  cd<-data.frame(y=yy,lp=qlogis(pmin(pmax(pp,1e-8),1-1e-8)))
  f0<-quiet_glm(y~offset(lp),cd,quasibinomial(),ww);f1<-quiet_glm(y~lp,cd,quasibinomial(),ww)
  val['cal_intercept']<-coef(f0)[1];val['cal_slope']<-coef(f1)['lp']
 }
 val
}
evaluate<-function(d,pred,pr,validation,variant='common',weighting='conditional',replicate=0,calibrate=TRUE){
 ans<-list()
 for(h in 1:2)for(m in dimnames(pred)[[2]]){
  v<-metric(d,pred[,m,h],pr,h,weighting,calibrate)
  ans[[length(ans)+1]]<-data.frame(cohort=unique(d$cohort),validation=validation,variant=variant,weighting=weighting,replicate=replicate,horizon=h,model=m,as.list(v),row.names=NULL)
 }
 do.call(rbind,ans)
}
run_cv<-function(d,raw=FALSE){
 pred<-array(NA_real_,c(nrow(d),3,2,5),dimnames=list(NULL,c('C','M0','M1'),c('visit1','visit2'),paste0('r',1:5)))
 for(r in 1:5){
  folds<-make_folds(d,seed+100*r+match(unique(d$cohort),names(paths)),TRUE)
  for(k in 1:5){fit<-fit_pipeline(d[folds!=k,],raw);pred[folds==k,,,r]<-predict_pipeline(fit,d[folds==k,])}
  cat('CV',unique(d$cohort),if(raw)'raw' else 'common','repeat',r,'complete\n');flush.console()
 }
 pred
}
calibration_groups<-function(d,p,pr,h,model){
 y<-endpoint(d,h);w<-get_weights(d,pr,h);gr<-pmin(5,ceiling(rank(p,ties.method='first')/length(p)*5));ans<-list()
 for(g in 1:5){ix<-which(gr==g);ok<-ix[!is.na(y[ix])];ans[[g]]<-data.frame(cohort=unique(d$cohort),horizon=h,model=model,group=g,n=length(ix),observed_n=length(ok),predicted=mean(p[ix]),observed=sum(w[ok]*y[ok])/sum(w[ok]),ess=sum(w[ok])^2/sum(w[ok]^2))}
 do.call(rbind,ans)
}
# Paired fixed-prediction bootstrap, with the same participants sampled across repeats.
fixed_bootstrap<-function(d,pred,pr,B=1000,seed=1){
 set.seed(seed);is_cv<-length(dim(pred))==4;nr<-if(is_cv)dim(pred)[4] else 1L;ans<-list()
 for(b in seq_len(B)){
  ix<-sample.int(nrow(d),replace=TRUE);dd<-d[ix,];qq<-pr[ix,,drop=FALSE]
  for(h in 1:2){
   v<-sapply(1:nr,function(r){pp<-if(is_cv)pred[ix,,,r,drop=FALSE] else pred[ix,,,drop=FALSE];dim(pp)<-c(length(ix),3,2);dimnames(pp)<-list(NULL,c('C','M0','M1'),c('visit1','visit2'));a<-metric(dd,pp[,'M0',h],qq,h,calibrate=FALSE);z<-metric(dd,pp[,'M1',h],qq,h,calibrate=FALSE);c(delta_auc=z['auc']-a['auc'],delta_brier=z['brier']-a['brier'])})
   a<-rowMeans(matrix(v,nrow=2));ans[[length(ans)+1]]<-data.frame(cohort=unique(d$cohort),validation=if(is_cv)'local_CV' else 'frozen_external',bootstrap=b,horizon=h,delta_auc=a[1],delta_brier=a[2])
  }
 }
 do.call(rbind,ans)
}
write_table<-function(d,name)write.csv(d,file.path(root,'results',name),row.names=FALSE,na='')
if(mode=='full'){
 data<-lapply(names(paths),load_data);names(data)<-names(paths)
 desc<-do.call(rbind,lapply(data,function(d){data.frame(cohort=unique(d$cohort),eligible=nrow(d),first_observed=sum(!is.na(d$y1)),second_at_risk=sum(!is.na(d$y1)&d$y1==0),second_observed=sum(!is.na(d$y1)&d$y1==0&!is.na(d$y2)),events1=sum(d$y1==1,na.rm=TRUE),events2=sum(endpoint(d,2)==1,na.rm=TRUE),known_endpoint2=sum(!is.na(endpoint(d,2))),missing_covariate=sum(!complete.cases(d[,features])),median_t1=median(d$t1,na.rm=TRUE),median_t2=median(d$t2,na.rm=TRUE))}))
 write_table(desc,'sample_description.csv')
 probs<-lapply(seq_along(data),function(i)observation_probs(data[[i]],seed+i));names(probs)<-names(data)
 fullfit<-fit_pipeline(data$HRS);pred<-lapply(data,function(d)predict_pipeline(fullfit,d))
 res<-list();cal<-list();wdiag<-list();cvs<-list();rawcvs<-list();boots<-list();ccobjects<-list()
 for(k in names(data)){
  d<-data[[k]];pr<-probs[[k]];val<-if(k=='HRS')'apparent_development' else 'frozen_external'
  for(wt in c('conditional','untruncated','marginal','unweighted'))res[[length(res)+1]]<-evaluate(d,pred[[k]],pr,val,weighting=wt)
  for(h in 1:2){ww<-get_weights(d,pr,h);un<-get_weights(d,pr,h,'untruncated');wdiag[[length(wdiag)+1]]<-data.frame(cohort=k,horizon=h,known=sum(ww>0),ess=sum(ww)^2/sum(ww^2),min_probability=min(pr),max_weight=max(ww),max_untruncated=max(un),weight_sum_ratio=sum(ww)/nrow(d),probability_clips=sum(pr<.02|pr>.98),weights_above20=sum(un>20))}
  for(m in c('M0','M1'))for(h in 1:2)cal[[length(cal)+1]]<-calibration_groups(d,pred[[k]][,m,h],pr,h,m)
  cvs[[k]]<-run_cv(d,FALSE);rawcvs[[k]]<-run_cv(d,TRUE)
  for(r in 1:5){res[[length(res)+1]]<-evaluate(d,cvs[[k]][,,,r],pr,'local_CV',replicate=r);res[[length(res)+1]]<-evaluate(d,rawcvs[[k]][,,,r],pr,'local_CV',variant='raw_depression_cohort_terms',replicate=r)}
  boots[[length(boots)+1]]<-fixed_bootstrap(d,cvs[[k]],pr,1000,seed+match(k,names(data)))
  if(k!='HRS')boots[[length(boots)+1]]<-fixed_bootstrap(d,pred[[k]],pr,1000,seed+10+match(k,names(data)))
  cat('Completed evaluations',k,'\n');flush.console()
 }
 # Complete-case sensitivity uses a new frozen HRS fit.
 cc<-lapply(data,function(d)d[complete.cases(d[,features]),,drop=FALSE]);ccfit<-fit_pipeline(cc$HRS)
 for(k in names(cc)){
  cp<-predict_pipeline(ccfit,cc[[k]]);cq<-observation_probs(cc[[k]],seed+500+match(k,names(cc)))
  res[[length(res)+1]]<-evaluate(cc[[k]],cp,cq,if(k=='HRS')'apparent_development' else 'frozen_external',variant='complete_predictor_cases')
  ccobjects[[k]]<-list(pred=cp,probs=cq)
 }
 write_table(do.call(rbind,res),'performance_all.csv');write_table(do.call(rbind,cal),'calibration_groups.csv');write_table(do.call(rbind,wdiag),'observation_weights.csv');write_table(do.call(rbind,boots),'conditional_bootstrap.csv')
 cf<-do.call(rbind,lapply(names(fullfit$fits),function(m)data.frame(model=m,term=names(coef(fullfit$fits[[m]])),estimate=as.numeric(coef(fullfit$fits[[m]])))))
 write_table(cf,'frozen_HRS_coefficients.csv')
 write_json(fullfit$imputer,file.path(root,'results/frozen_HRS_imputation.json'),pretty=TRUE,auto_unbox=TRUE)
 saveRDS(list(data=data,probs=probs,fit=fullfit,pred=pred,cvs=cvs,rawcvs=rawcvs,cc=ccobjects),file.path(root,'data_private/analysis_objects.rds'))
 writeLines(unique(allwarnings),file.path(root,'audit/model_warnings_full.txt'));writeLines(capture.output(sessionInfo()),file.path(root,'audit/R_session.txt'))
 writeLines(paste('Completed main module',format(Sys.time(),'%Y-%m-%d %H:%M:%S %z'),'GLM fits',fit_count),file.path(root,'audit/main_complete.txt'))
 cat('MAIN MODULE COMPLETE; GLM fits',fit_count,'\n')
}
if(mode=='bootstrap'){
 ob<-readRDS(file.path(root,'data_private/analysis_objects.rds'));data<-ob$data;ans<-list();fail<-list();B<-500L
 for(b in seq_len(B)){
  set.seed(seed+100000+b);tr<-data$HRS[sample.int(nrow(data$HRS),replace=TRUE),]
  fit<-tryCatch(fit_pipeline(tr),error=function(e)e)
  if(inherits(fit,'error')){fail[[length(fail)+1]]<-data.frame(bootstrap=b,cohort='HRS',reason=conditionMessage(fit));next}
  for(k in c('ELSA','SHARE')){
   dd<-data[[k]][sample.int(nrow(data[[k]]),replace=TRUE),]
   rr<-tryCatch({pp<-predict_pipeline(fit,dd);qq<-observation_probs(dd,seed+200000+b+1000*match(k,names(data)));ev<-evaluate(dd,pp,qq,'two_sample_bootstrap',replicate=b);ev$bootstrap<-b;ev},error=function(e)e)
   if(inherits(rr,'error'))fail[[length(fail)+1]]<-data.frame(bootstrap=b,cohort=k,reason=conditionMessage(rr)) else ans[[length(ans)+1]]<-rr
  }
  if(b%%25==0){cat('Two-sample bootstrap',b,'of',B,'\n');flush.console();write_table(do.call(rbind,ans),'two_sample_bootstrap.csv')}
 }
 write_table(do.call(rbind,ans),'two_sample_bootstrap.csv')
 if(length(fail))write_table(do.call(rbind,fail),'bootstrap_failures.csv')
 writeLines(unique(allwarnings),file.path(root,'audit/model_warnings_bootstrap.txt'));if(length(warning_records))write.csv(do.call(rbind,warning_records),file.path(root,'audit/bootstrap_warning_details.csv'),row.names=FALSE)
 writeLines(paste('Completed two-sample bootstrap',format(Sys.time(),'%Y-%m-%d %H:%M:%S %z'),'GLM fits',fit_count,'failed units',length(fail)),file.path(root,'audit/bootstrap_complete.txt'))
 cat('BOOTSTRAP COMPLETE; failed units',length(fail),'\n')
}
