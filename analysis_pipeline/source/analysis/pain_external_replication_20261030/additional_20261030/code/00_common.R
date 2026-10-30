Sys.setenv(LANGUAGE='en')
suppressPackageStartupMessages({library(mice);library(sandwich);library(metafor)})
P<-normalizePath('analysis/pain_external_replication_20261030/additional_20261030',winslash='/',mustWork=TRUE)
EXT<-dirname(P);OLD<-file.path(dirname(EXT),'pain_transition_diabetes_adl_20261029')
covars<-c('age','female','education','smoking','arthritis','cvd','cesd')
format_data<-function(d,cohort){
 d$exposure<-factor(d$exposure,levels=c('00','01','10','11'));d$education<-factor(d$education,levels=1:3)
 for(v in c('female','smoking','arthritis','cvd'))d[[v]]<-factor(d[[v]],levels=0:1)
 if(cohort=='SHARE')d$country<-factor(d$country) else d$entry<-factor(d$entry,levels=c(2010,2012));d
}
get_formula<-function(cohort,extra=character())as.formula(paste('y ~',paste(c('exposure',covars,if(cohort=='SHARE')'country' else 'entry','interval',extra),collapse='+'),'+offset(log(years))'))
intervals<-function(d,risk=FALSE,composite=FALSE,nominal=FALSE){
 z<-d
 if(composite)for(j in 1:2){v<-paste0('y',j);dead<-z[[paste0('dead',j)]];z[[v]][!is.na(dead)&dead==1]<-1}
 a<-z;a$interval<-factor(1,levels=1:2);a$y<-a$y1;a$years<-if(nominal)2 else a$t1
 b<-z[!is.na(z$y1)&z$y1==0,,drop=FALSE];b$interval<-factor(2,levels=1:2);b$y<-b$y2;b$years<-if(nominal)2 else b$t2-b$t1
 q<-rbind(a,b);q$observed<-as.integer(!is.na(q$y));q<-q[order(q$pid,q$interval),];rownames(q)<-NULL
 if(!risk)q<-q[q$observed==1,];q
}
diags<-list();warnings_log<-list();results<-list();allfits<-list()
fit_one<-function(q,ff,variant,i,w=NULL){
 q$.w<-if(is.null(w))1 else w;msgs<-character()
 f<-withCallingHandlers(glm(ff,data=q,family=if(is.null(w))binomial('cloglog') else quasibinomial('cloglog'),weights=.w,control=glm.control(maxit=100)),warning=function(x){msgs<<-c(msgs,conditionMessage(x));invokeRestart('muffleWarning')})
 V<-vcovCL(f,cluster=q$pid,type='HC1',cadjust=TRUE)
 stopifnot(f$converged,f$rank==length(coef(f)),all(is.finite(coef(f))),all(is.finite(V)),all(diag(V)>0))
 diags[[length(diags)+1]]<<-data.frame(cohort=cohort,variant=variant,imputation=i,persons=length(unique(q$pid)),events=sum(q$y),intervals=nrow(q),parameters=length(coef(f)),converged=f$converged,max_leverage=max(hatvalues(f)),min_cov_eigen=min(eigen(V,symmetric=TRUE,only.values=TRUE)$values),warnings=length(msgs))
 if(length(msgs))warnings_log[[length(warnings_log)+1]]<<-data.frame(variant=variant,imputation=i,message=unique(msgs))
 list(b=coef(f),V=V,n=length(unique(q$pid)),k=length(coef(f)),events=sum(q$y),formula=ff)
}
contrast<-function(nms,type){v<-setNames(rep(0,length(nms)),nms);g<-strsplit(type,'_vs_')[[1]];for(j in 1:2)if(g[j]!='00')v[paste0('exposure',g[j])]<-if(j==1)1 else -1;v}
pool_effects<-function(fits,variant){
 rr<-lapply(c('11_vs_00','01_vs_00','10_vs_00','10_vs_11','11_vs_01'),function(type){
 v<-contrast(names(fits[[1]]$b),type);Q<-sapply(fits,function(f)sum(v*f$b));U<-sapply(fits,function(f)as.numeric(t(v)%*%f$V%*%v));n<-min(sapply(fits,`[[`,'n'));k<-max(sapply(fits,`[[`,'k'))
 if(length(fits)>1){ps<-pool.scalar(Q,U,n=n,k=k);est<-ps$qbar;se<-sqrt(ps$t);df<-ps$df;fmi<-ps$fmi}else{est<-Q;se<-sqrt(U);df<-n-k;fmi<-0}
 data.frame(cohort=cohort,variant=variant,contrast=type,persons=n,log_estimate=est,se=se,estimate=exp(est),lower=exp(est-qt(.975,df)*se),upper=exp(est+qt(.975,df)*se),p=2*pt(abs(est/se),df,lower.tail=FALSE),df=df,fmi=fmi)
 });rr<-do.call(rbind,rr);rr$p_holm_secondary<-NA_real_;rr$p_holm_secondary[-1]<-p.adjust(rr$p[-1],'holm');results[[length(results)+1]]<<-rr;allfits[[variant]]<<-fits;invisible(rr)
}
load_completed<-function(cohort){
 if(cohort=='SHARE'){d<-read.csv(file.path(P,'data_private/SHARE_extended.csv'),colClasses=c(exposure='character'));imp<-readRDS(file.path(EXT,'data_private/SHARE_mice50.rds'))}else{d<-read.csv(file.path(OLD,'data_private',paste0(cohort,'_longitudinal.csv')),colClasses=c(exposure='character'));imp<-readRDS(file.path(OLD,'data_private',paste0(cohort,'_longitudinal_mice50.rds')))}
 d<-format_data(d,cohort);lapply(1:imp$m,function(i){dd<-d;cc<-complete(imp,i);for(v in names(cc))dd[[v]]<-cc[[v]];dd})
}
