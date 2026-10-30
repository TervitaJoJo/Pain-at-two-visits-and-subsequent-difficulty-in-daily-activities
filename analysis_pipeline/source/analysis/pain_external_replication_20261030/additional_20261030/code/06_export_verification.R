source('analysis/pain_external_replication_20261030/additional_20261030/code/00_common.R')
cohort<-'SHARE';ff<-readRDS(file.path(P,'data_private/SHARE_additional_fits.rds'));rows<-list()
for(v in names(ff))for(i in seq_along(ff[[v]]))for(ct in c('11_vs_00','01_vs_00','10_vs_00','10_vs_11','11_vs_01')){
 f<-ff[[v]][[i]];L<-contrast(names(f$b),ct)
 rows[[length(rows)+1]]<-data.frame(variant=v,imputation=i,contrast=ct,Q=sum(L*f$b),U=as.numeric(t(L)%*%f$V%*%L),n=f$n,k=f$k)
}
write.csv(do.call(rbind,rows),file.path(P,'results/SHARE_scalar_inputs.csv'),row.names=FALSE)
d<-format_data(read.csv(file.path(P,'data_private/SHARE_extended.csv'),colClasses=c(exposure='character')),cohort)
imp<-readRDS(file.path(P,'data_private/SHARE_extended_mice50.rds'));z<-complete(imp,1);for(v in names(z))d[[v]]<-z[[v]]
q<-intervals(d);form<-get_formula(cohort,c('diabetes_medication','bmi','diabetes_duration'));f<-glm(form,data=q,family=binomial('cloglog'),control=glm.control(maxit=100))
X<-model.matrix(f);V<-vcovCL(f,cluster=q$pid,type='HC1',cadjust=TRUE)
write.csv(cbind(data.frame(pid=q$pid,y=q$y,offset=log(q$years)),as.data.frame(X,check.names=FALSE)),file.path(P,'data_private/context_verification_design.csv'),row.names=FALSE)
write.csv(data.frame(term=names(coef(f)),estimate=coef(f),se=sqrt(diag(V))),file.path(P,'results/context_verification_coefficients.csv'),row.names=FALSE)
stopifnot(max(abs(coef(f)-ff$diabetes_context_MI[[1]]$b))<1e-12)
write.csv(data.frame(variable=names(imp$method),method=imp$method,missing=colSums(is.na(imp$data))),file.path(P,'audit/SHARE_imputation_methods.csv'),row.names=FALSE)
cat('Exported summary-only pooling inputs and PRIVATE independent-design check\n')
