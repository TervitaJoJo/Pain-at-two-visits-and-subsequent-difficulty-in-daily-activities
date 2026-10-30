Sys.setenv(LANGUAGE='en')
args<-commandArgs(trailingOnly=TRUE)
cohort<-args[1]
original<-normalizePath('analysis/pain_transition_diabetes_adl_20261029',winslash='/',mustWork=TRUE)
code<-readLines(file.path(original,'code/02_analyze.R'),warn=FALSE)
end<-grep("^cat\\(cohort,'building imputations",code)[1]-1L
eval(parse(text=code[1:end]))
rev<-file.path(original,'revision_20261029'); out<-file.path(rev,'results');audit<-file.path(rev,'audit')
seed<-if(cohort=='HRS')202610295L else 202610296L
readrev<-function(){
 d<-read.csv(file.path(rev,'data_private',paste0(cohort,'_longitudinal.csv')),colClasses=c(exposure='character'))
 d$exposure<-factor(d$exposure,levels=c('00','01','10','11'));d$education<-factor(d$education,levels=1:3);d$entry<-factor(d$entry,levels=c(2010,2012))
 for(v in c('female','smoking','arthritis','cvd')) d[[v]]<-factor(d[[v]],levels=0:1)
 d$diabetes_treatment<-factor(d$diabetes_treatment,levels=0:2);d
}
d<-readrev()
origimp<-readRDS(file.path(original,'data_private',paste0(cohort,'_longitudinal_mice50.rds')))
origcomp<-lapply(1:50,function(i){dd<-d;cc<-complete(origimp,i);for(v in names(cc)) dd[[v]]<-cc[[v]];dd})
form<-as.formula(paste('y ~',rhs,'+ offset(log(years))'))
origfits<-readRDS(file.path(original,'data_private',paste0(cohort,'_models.rds')))$primary
pool_contrasts(origfits,'primary_history_added',types=c('10_vs_00','11_vs_01'))
L<-rbind(contrast_vector(names(origfits[[1]]$b),'10_vs_00'),contrast_vector(names(origfits[[1]]$b),'11_vs_01'))
Q<-sapply(origfits,function(f)as.vector(L%*%f$b));U<-array(NA_real_,c(2,2,50))
for(i in 1:50) U[,,i]<-L%*%origfits[[i]]$V%*%t(L)
joint<-getFromNamespace('.D1','mitml')(Q,U,df.com=origfits[[1]]$n-origfits[[1]]$k)
write.csv(data.frame(cohort=cohort,hypothesis='beta10=0; beta11-beta01=0',F=as.numeric(joint$F),df1=joint$k,df2=joint$v,p=pf(joint$F,joint$k,joint$v,lower.tail=FALSE),riv=joint$r,method='mitml D1 with participant-clustered robust within-imputation covariance'),file.path(out,paste0(cohort,'_history_joint.csv')),row.names=FALSE)
cat(cohort,'history comparison done; fitting proxy and nominal-time analyses\n');flush.console()
fits<-lapply(1:50,function(i){dd<-origcomp[[i]];bad1<-!is.na(dd$proxy1)&dd$proxy1==1;bad2<-bad1|(!is.na(dd$proxy2)&dd$proxy2==1);dd$y1[bad1]<-NA;dd$y2[bad2]<-NA;fit_glm(intervals(dd),form,binomial('cloglog'),'proxy_censored_original_MI',i)})
pool_contrasts(fits,'proxy_censored_original_MI',types=c('11_vs_00','01_vs_00','10_vs_00','10_vs_11','11_vs_01'))
if(cohort=='HRS'){
 fits<-lapply(1:50,function(i){dd<-intervals(origcomp[[i]]);dd$years<-2;fit_glm(dd,form,binomial('cloglog'),'ADL_only_nominal_original_MI',i)})
 pool_contrasts(fits,'ADL_only_nominal_original_MI',types=c('11_vs_00','01_vs_00','10_vs_00','10_vs_11','11_vs_01'))
}
obs<-intervals(d);d$observed_event<-as.integer(d$pid%in%obs$pid[obs$y==1]);cnt<-table(obs$pid);yrs<-tapply(obs$years,obs$pid,sum)
d$n_intervals<-as.numeric(cnt[as.character(d$pid)]);d$n_intervals[is.na(d$n_intervals)]<-0
d$followup_years<-as.numeric(yrs[as.character(d$pid)]);d$followup_years[is.na(d$followup_years)]<-0
extra<-c('diabetes_treatment','landmark_function3',if(cohort=='HRS')c('bmi','diabetes_duration'))
cols<-c('exposure',covars,'entry','function3',extra,'observed_event','n_intervals','followup_years')
dat<-d[,cols];met<-setNames(rep('',length(cols)),cols)
for(v in cols)if(anyNA(dat[[v]]))met[v]<-if(v%in%c('education','diabetes_treatment'))'polyreg' else if(v%in%c('female','smoking','arthritis','cvd'))'logreg' else 'pmm'
pred<-matrix(1,length(cols),length(cols),dimnames=list(cols,cols));diag(pred)<-0
# Restrict duration donors to values compatible with recipient age (calendar rounding +1 year).
mice.impute.duration_pmm <- function(y,ry,x,wy=NULL,...) {
 if(is.null(wy))wy<-!ry
 values<-mice::mice.impute.pmm(y,ry,x,wy=wy,...)
 age_limit<-x[wy,'age']+1
 for(attempt in 1:100) {
  invalid<-values>age_limit|values<0
  if(!any(invalid))return(values)
  proposal<-mice::mice.impute.pmm(y,ry,x,wy=wy,...)
  values[invalid]<-proposal[invalid]
 }
 stop('Cannot impute a biologically admissible duration using predictive mean matching')
}
if(cohort=='HRS')met['diabetes_duration']<-'duration_pmm'
cat(cohort,'starting expanded imputations\n');flush.console()
imp<-mice(dat,m=50,maxit=20,method=met,predictorMatrix=pred,seed=seed,printFlag=FALSE,donors=5)
saveRDS(imp,file.path(rev,'data_private',paste0(cohort,'_revision_mice50.rds')))
if(!is.null(imp$loggedEvents))write.csv(imp$loggedEvents,file.path(audit,paste0(cohort,'_revision_mice_events.csv')),row.names=FALSE)
write.csv(mice::convergence(imp),file.path(audit,paste0(cohort,'_revision_mice_convergence.csv')),row.names=FALSE)
pdf(file.path(rev,'figures',paste0(cohort,'_revision_imputation_chains.pdf')),width=10,height=7);print(plot(imp));dev.off()
completed<-lapply(1:50,function(i){dd<-d;cc<-complete(imp,i);for(v in cols){known<-!is.na(dat[[v]]);stopifnot(all(as.character(cc[[v]][known])==as.character(dat[[v]][known])));dd[[v]]<-cc[[v]]};stopifnot(!anyNA(dd[,c(covars,extra)]));if(cohort=='HRS')stopifnot(all(dd$diabetes_duration>=0 & dd$diabetes_duration<=dd$age+1));dd})
for(v in names(met)[met!='']){
 vals<-unlist(lapply(completed,function(dd)as.numeric(as.character(dd[[v]][is.na(dat[[v]])]))))
 impdiag[[length(impdiag)+1]]<-data.frame(cohort=cohort,variable=v,method=met[v],missing=sum(is.na(dat[[v]])),imputed_min=min(vals),imputed_max=max(vals),imputed_mean=mean(vals),imputed_sd=sd(vals))
}
forms<-list(augmented_imputation_base=form,treatment_adjusted_MI=update(form,.~.+diabetes_treatment),landmark_function_adjusted_MI=update(form,.~.+landmark_function3))
if(cohort=='HRS')forms$diabetes_context_MI<-update(form,.~.+diabetes_treatment+bmi+diabetes_duration)
allfits<-list()
for(variant in names(forms)){
 cat(cohort,variant,'\n');flush.console()
 fits<-lapply(1:50,function(i)fit_glm(intervals(completed[[i]]),forms[[variant]],binomial('cloglog'),variant,i))
 pool_contrasts(fits,variant,types=c('11_vs_00','01_vs_00','10_vs_00','10_vs_11','11_vs_01'));allfits[[variant]]<-fits
}
write.csv(do.call(rbind,results),file.path(out,paste0(cohort,'_revision_contrasts.csv')),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(audit,paste0(cohort,'_revision_model_diagnostics.csv')),row.names=FALSE)
write.csv(do.call(rbind,impdiag),file.path(audit,paste0(cohort,'_revision_imputation_summary.csv')),row.names=FALSE)
if(length(modelwarnings))write.csv(do.call(rbind,modelwarnings),file.path(audit,paste0(cohort,'_revision_model_warnings.csv')),row.names=FALSE)
saveRDS(list(fits=allfits,joint=joint,results=results),file.path(rev,'data_private',paste0(cohort,'_revision_models.rds')))
writeLines(capture.output(sessionInfo()),file.path(audit,paste0(cohort,'_revision_R_session.txt')))
cat(cohort,'REVISION ANALYSES COMPLETE\n')
