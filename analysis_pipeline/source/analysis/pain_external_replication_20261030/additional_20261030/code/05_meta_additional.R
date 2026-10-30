source('analysis/pain_external_replication_20261030/additional_20261030/code/00_common.R')
d<-read.csv(file.path(P,'results/primary_and_function_contrasts.csv'))
res<-list();inputs<-list()
for(v in unique(d$variant))for(ct in unique(d$contrast)){
 z<-d[d$variant==v&d$contrast==ct,];stopifnot(nrow(z)==3,setequal(z$cohort,c('HRS','ELSA','SHARE')))
 f<-rma.uni(yi=z$log_estimate,sei=z$se,method='REML',test='adhoc',control=list(threshold=1e-12,maxiter=1000))
 res[[length(res)+1]]<-data.frame(variant=v,contrast=ct,k=3,logHR=as.numeric(coef(f)),se=as.numeric(f$se),HR=exp(as.numeric(coef(f))),lower=exp(f$ci.lb),upper=exp(f$ci.ub),p=f$pval,tau2=f$tau2,I2=f$I2,Q_p=f$QEp,HK_scale=f$s2w)
 z$weight_percent<-as.numeric(weights(f));inputs[[length(inputs)+1]]<-z
}
r<-do.call(rbind,res);r$p_holm_secondary<-NA_real_
for(v in unique(r$variant)){ix<-r$variant==v&r$contrast!='11_vs_00';r$p_holm_secondary[ix]<-p.adjust(r$p[ix],'holm')}
write.csv(r,file.path(P,'results/three_cohort_meta.csv'),row.names=FALSE)
write.csv(do.call(rbind,inputs),file.path(P,'results/three_cohort_meta_inputs.csv'),row.names=FALSE)
print(r[,c('variant','contrast','HR','lower','upper','p','p_holm_secondary','I2')],row.names=FALSE)
