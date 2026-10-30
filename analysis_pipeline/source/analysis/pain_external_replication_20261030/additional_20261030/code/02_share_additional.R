source('analysis/pain_external_replication_20261030/additional_20261030/code/00_common.R')
cohort<-'SHARE';completed<-load_completed(cohort);form<-get_formula(cohort)
primary<-readRDS(file.path(EXT,'data_private/SHARE_primary_fits.rds'));pool_effects(primary,'primary_history_added')
L<-rbind(contrast(names(primary[[1]]$b),'10_vs_00'),contrast(names(primary[[1]]$b),'11_vs_01'))
Q<-sapply(primary,function(f)as.vector(L%*%f$b));U<-array(NA_real_,c(2,2,50));for(i in 1:50)U[,,i]<-L%*%primary[[i]]$V%*%t(L)
joint<-getFromNamespace('.D1','mitml')(Q,U,df.com=primary[[1]]$n-primary[[1]]$k)
write.csv(data.frame(cohort=cohort,F=joint$F,df1=joint$k,df2=joint$v,p=pf(joint$F,joint$k,joint$v,lower.tail=FALSE),riv=joint$r),file.path(P,'results/SHARE_history_joint.csv'),row.names=FALSE)
for(variant in c('proxy_censored_MI','proxy_censored_strict_MI','ADL_only_nominal_MI','ADL_or_death_nominal_MI')){
 cat(variant,'\n');flush.console()
 fs<-lapply(1:50,function(i){d<-completed[[i]]
 if(grepl('proxy',variant)){
  bad1<-!is.na(d$proxy1)&d$proxy1==1;bad2<-!is.na(d$proxy2)&d$proxy2==1
  if(grepl('strict',variant)){bad1<-bad1|is.na(d$proxy1);bad2<-bad2|is.na(d$proxy2)}
  d$y1[bad1]<-NA;d$y2[bad1|bad2]<-NA
 }
 fit_one(intervals(d,composite=variant=='ADL_or_death_nominal_MI',nominal=grepl('nominal',variant)),form,variant,i)
 });pool_effects(fs,variant)
}
weight_diags<-list();weighted_data<-list()
for(i in 1:50){
 r<-intervals(completed[[i]],risk=TRUE)
 denform<-as.formula(paste('observed ~',paste(c('exposure',covars,'country','interval'),collapse='+')))
 den<-glm(denform,data=r,family=binomial())
 num<-glm(observed~exposure+interval+country,data=r,family=binomial())
 stopifnot(den$converged,num$converged)
 pd<-predict(den,type='response');pn<-predict(num,type='response');stopifnot(nobs(den)==nrow(r),nobs(num)==nrow(r),length(pd)==nrow(r),all(pd>0&pd<1))
 w<-ave(pn/pd,r$pid,FUN=cumprod);keep<-r$observed==1;q<-r[keep,];w<-w[keep];lim<-quantile(w,c(.01,.99));wt<-pmin(pmax(w,lim[1]),lim[2]);weighted_data[[i]]<-list(q=q,w=w,wt=wt)
 for(j in 0:2){ix<-if(j==0)rep(TRUE,length(w))else q$interval==j;ir<-if(j==0)rep(TRUE,nrow(r))else r$interval==j
 weight_diags[[length(weight_diags)+1]]<-data.frame(imputation=i,interval=j,risk_intervals=sum(ir),observed_intervals=sum(ix),min_den=min(pd[ir]),max_den=max(pd[ir]),raw_min=min(w[ix]),raw_max=max(w[ix]),trunc_min=min(wt[ix]),trunc_max=max(wt[ix]),weight_mean=mean(wt[ix]),weight_sd=sd(wt[ix]),ESS=sum(wt[ix])^2/sum(wt[ix]^2),truncated=sum(w[ix]!=wt[ix]))
 }
}
for(v in c('observation_weighted_MI','observation_weighted_untruncated_MI')){
 fs<-lapply(1:50,function(i){z<-weighted_data[[i]];fit_one(z$q,form,v,i,if(v=='observation_weighted_MI')z$wt else z$w)});pool_effects(fs,v)
}
write.csv(do.call(rbind,weight_diags),file.path(P,'audit/SHARE_weight_diagnostics.csv'),row.names=FALSE)
# Additional imputation model, including clinically measured baseline variables.
d<-format_data(read.csv(file.path(P,'data_private/SHARE_extended.csv'),colClasses=c(exposure='character')),cohort)
obs<-intervals(d);d$observed_event<-as.integer(d$pid%in%obs$pid[obs$y==1]);ct<-table(obs$pid);yrs<-tapply(obs$years,obs$pid,sum)
d$n_intervals<-as.numeric(ct[as.character(d$pid)]);d$n_intervals[is.na(d$n_intervals)]<-0;d$followup_years<-as.numeric(yrs[as.character(d$pid)]);d$followup_years[is.na(d$followup_years)]<-0
d$diabetes_medication<-factor(d$diabetes_medication,levels=0:1)
cols<-c('exposure',covars,'country','function3','landmark_function3','diabetes_medication','bmi','diabetes_duration','observed_event','n_intervals','followup_years')
dat<-d[,cols];met<-setNames(rep('',length(cols)),cols)
for(v in cols)if(anyNA(dat[[v]]))met[v]<-if(v=='education')'polyreg' else if(is.factor(dat[[v]]))'logreg' else 'pmm'
mice.impute.duration_pmm<-function(y,ry,x,wy=NULL,...){
 if(is.null(wy))wy<-!ry;values<-mice::mice.impute.pmm(y,ry,x,wy=wy,...)
 for(attempt in 1:100){bad<-values<0|values>x[wy,'age']+1;if(!any(bad))return(values);proposal<-mice::mice.impute.pmm(y,ry,x,wy=wy,...);values[bad]<-proposal[bad]};stop('Invalid duration imputation')
}
met['diabetes_duration']<-'duration_pmm';pred<-matrix(1,length(cols),length(cols),dimnames=list(cols,cols));diag(pred)<-0
cat('Starting expanded MI 50 x 20\n');flush.console()
cache<-file.path(P,'data_private/SHARE_extended_mice50.rds')
imp<-if(file.exists(cache)&&identical(Sys.getenv('REUSE_MI_CACHE'),'1'))readRDS(cache)else mice(dat,m=50,maxit=20,method=met,predictorMatrix=pred,donors=5,seed=202610303,printFlag=FALSE)
saveRDS(imp,file.path(P,'data_private/SHARE_extended_mice50.rds'))
write.csv(mice::convergence(imp),file.path(P,'audit/SHARE_extended_mice_convergence.csv'),row.names=FALSE)
if(!is.null(imp$loggedEvents))write.csv(imp$loggedEvents,file.path(P,'audit/SHARE_extended_mice_events.csv'),row.names=FALSE)
cc<-lapply(1:50,function(i){z<-d;v<-complete(imp,i);for(n in names(v)){known<-!is.na(dat[[n]]);stopifnot(all(as.character(dat[[n]][known])==as.character(v[[n]][known])));z[[n]]<-v[[n]]};stopifnot(all(z$diabetes_duration>=0&z$diabetes_duration<=z$age+1));z})
forms<-list(augmented_imputation_base=form,treatment_adjusted_MI=update(form,.~.+diabetes_medication),treatment_BMI_MI=update(form,.~.+diabetes_medication+bmi),diabetes_context_MI=update(form,.~.+diabetes_medication+bmi+diabetes_duration),landmark_function_adjusted_MI=update(form,.~.+landmark_function3))
for(v in names(forms)){cat(v,'\n');flush.console();pool_effects(lapply(1:50,function(i)fit_one(intervals(cc[[i]]),forms[[v]],v,i)),v)}
extra<-c(covars,'diabetes_medication','bmi','diabetes_duration');z<-d[complete.cases(d[,extra]),]
pool_effects(list(fit_one(intervals(z),forms$diabetes_context_MI,'diabetes_context_CC',1)),'diabetes_context_CC')
for(cty in levels(d$country)){
 v<-paste0('omit_country_',cty);fs<-lapply(1:50,function(i)fit_one(intervals(completed[[i]][completed[[i]]$country!=cty,]),form,v,i));pool_effects(fs,v)
}
write.csv(do.call(rbind,results),file.path(P,'results/SHARE_additional_contrasts.csv'),row.names=FALSE)
write.csv(do.call(rbind,diags),file.path(P,'audit/SHARE_model_diagnostics.csv'),row.names=FALSE)
if(length(warnings_log))write.csv(do.call(rbind,warnings_log),file.path(P,'audit/SHARE_model_warnings.csv'),row.names=FALSE)
saveRDS(allfits,file.path(P,'data_private/SHARE_additional_fits.rds'))
writeLines(capture.output(sessionInfo()),file.path(P,'audit/R_sessionInfo.txt'))
cat('SHARE additional analyses completed\n')
