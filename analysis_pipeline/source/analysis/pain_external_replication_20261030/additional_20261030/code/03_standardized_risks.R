source('analysis/pain_external_replication_20261030/additional_20261030/code/00_common.R')
rows<-list();scalar_rows<-list();gradient_checks<-list();support<-list()
for(cohort in c('HRS','ELSA','SHARE')){
 cat(cohort,'standardization\n');flush.console();completed<-load_completed(cohort);form<-get_formula(cohort)
 primary<-if(cohort=='SHARE')readRDS(file.path(EXT,'data_private/SHARE_primary_fits.rds')) else readRDS(file.path(OLD,'data_private',paste0(cohort,'_models.rds')))$primary
 func<-lapply(1:50,function(i)fit_one(intervals(completed[[i]]),update(form,.~.+function3),'function_adjusted_MI',i))
 pool_effects(func,'function_adjusted_MI');pool_effects(primary,'primary_MI')
 d<-completed[[1]]
 for(v in c('t1','t2')){z<-d[[v]];z<-z[is.finite(z)];support[[length(support)+1]]<-data.frame(cohort=cohort,time=v,n=length(z),min=min(z),q10=quantile(z,.1),median=median(z),q90=quantile(z,.9),max=max(z))}
 for(variant in c('primary_MI','function_adjusted_MI')){
  fits<-if(variant=='primary_MI')primary else func;ff<-if(variant=='primary_MI')form else update(form,.~.+function3)
  for(i in 1:50){
   fit<-fits[[i]];dd<-completed[[i]];allest<-list()
   for(g in c('00','01','10','11')){
    dd$exposure<-factor(g,levels=c('00','01','10','11'));dd$years<-2;dd$interval<-factor(1,levels=1:2)
    X1<-model.matrix(delete.response(terms(ff)),dd)[,names(fit$b),drop=FALSE]
    dd$interval<-factor(2,levels=1:2);X2<-model.matrix(delete.response(terms(ff)),dd)[,names(fit$b),drop=FALSE]
    H1<-2*exp(X1%*%fit$b);H2<-2*exp(X2%*%fit$b);S<-exp(-H1-H2);risk<-mean(1-S)
    grad<-colMeans(as.numeric(S*H1)*X1+as.numeric(S*H2)*X2)
    allest[[g]]<-list(r=risk,g=grad)
    vv<-as.numeric(t(grad)%*%fit$V%*%grad)/(risk*(1-risk))^2
    scalar_rows[[length(scalar_rows)+1]]<-data.frame(cohort=cohort,variant=variant,imputation=i,estimand=paste0('risk_',g),scale='logit',Q=qlogis(risk),U=vv,n=fit$n,k=fit$k,eligible=nrow(dd))
    if(i==1){ng<-numDeriv::grad(function(b)mean(1-exp(-2*exp(X1%*%b)-2*exp(X2%*%b))),fit$b);err<-max(abs(ng-grad));stopifnot(err<1e-6);gradient_checks[[length(gradient_checks)+1]]<-data.frame(cohort=cohort,variant=variant,group=g,max_gradient_error=err)}
   }
   for(type in c('11_vs_00','01_vs_00','10_vs_00','11_vs_01')){
    gg<-strsplit(type,'_vs_')[[1]];a<-allest[[gg[1]]];b<-allest[[gg[2]]];gd<-a$g-b$g
    scalar_rows[[length(scalar_rows)+1]]<-data.frame(cohort=cohort,variant=variant,imputation=i,estimand=paste0('RD_',type),scale='identity',Q=a$r-b$r,U=as.numeric(t(gd)%*%fit$V%*%gd),n=fit$n,k=fit$k,eligible=nrow(dd))
   }
  }
 }
}
s<-do.call(rbind,scalar_rows)
for(key in unique(paste(s$cohort,s$variant,s$estimand))){
 z<-s[paste(s$cohort,s$variant,s$estimand)==key,];stopifnot(nrow(z)==50)
 ps<-pool.scalar(z$Q,z$U,n=min(z$n),k=max(z$k));lo<-ps$qbar-qt(.975,ps$df)*sqrt(ps$t);hi<-ps$qbar+qt(.975,ps$df)*sqrt(ps$t)
 tr<-if(z$scale[1]=='logit')plogis else identity
 rows[[length(rows)+1]]<-data.frame(cohort=z$cohort[1],variant=z$variant[1],estimand=z$estimand[1],horizon_years=4,eligible=z$eligible[1],estimate=tr(ps$qbar),lower=tr(lo),upper=tr(hi),se_pooled_scale=sqrt(ps$t),df=ps$df,fmi=ps$fmi,p=if(z$scale[1]=='identity')2*pt(abs(ps$qbar/sqrt(ps$t)),ps$df,lower.tail=FALSE)else NA_real_)
}
write.csv(do.call(rbind,rows),file.path(P,'results/standardized_net_risks.csv'),row.names=FALSE)
write.csv(s,file.path(P,'results/standardization_scalar_inputs.csv'),row.names=FALSE)
write.csv(do.call(rbind,gradient_checks),file.path(P,'audit/standardization_gradient_checks.csv'),row.names=FALSE)
write.csv(do.call(rbind,support),file.path(P,'audit/followup_support.csv'),row.names=FALSE)
write.csv(do.call(rbind,results),file.path(P,'results/primary_and_function_contrasts.csv'),row.names=FALSE)
write.csv(do.call(rbind,diags),file.path(P,'audit/function_model_diagnostics.csv'),row.names=FALSE)
cat('Standardization completed\n')
