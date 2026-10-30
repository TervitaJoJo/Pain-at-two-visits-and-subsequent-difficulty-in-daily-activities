library(mice)
a <- normalizePath('analysis/pain_transition_diabetes_adl_20261029/revision_20261029',winslash='/')
checks <- list()
for (cohort in c('HRS','ELSA')) {
 obj <- readRDS(file.path(a,'data_private',paste0(cohort,'_revision_models.rds')))
 for (variant in names(obj$fits)) {
  fs <- obj$fits[[variant]]
  z <- read.csv(file.path(a,'results',paste0(cohort,'_revision_contrasts.csv')))
  for(ct in c('11_vs_00','01_vs_00','10_vs_00','10_vs_11','11_vs_01')) {
   lev <- strsplit(ct,'_vs_')[[1]]
   vv <- setNames(rep(0,length(fs[[1]]$b)),names(fs[[1]]$b))
   if(lev[1]!='00') vv[paste0('exposure',lev[1])]<-1
   if(lev[2]!='00') vv[paste0('exposure',lev[2])]<- -1
   q <- vapply(fs,function(f)sum(vv*f$b),numeric(1))
   u <- vapply(fs,function(f)as.numeric(t(vv)%*%f$V%*%vv),numeric(1))
   meanq <- mean(q); total <- mean(u)+(1+1/length(q))*var(q)
   zz<-z[z$variant==variant & z$contrast==ct,]
   stopifnot(abs(meanq-zz$log_estimate)<1e-10,abs(sqrt(total)-zz$se)<1e-10)
   checks[[length(checks)+1]]<-data.frame(cohort,variant,contrast=ct,log_diff=meanq-zz$log_estimate,se_diff=sqrt(total)-zz$se)
  }
 }
 imp<-readRDS(file.path(a,'data_private',paste0(cohort,'_revision_mice50.rds')))
 if(cohort=='HRS') {
  bad<-sapply(1:50,function(i){d<-complete(imp,i);sum(d$diabetes_duration>d$age+1)})
  cat(cohort,'invalid imputed duration total',sum(bad),'maximum per imputation',max(bad),'\n')
 }
}
write.csv(do.call(rbind,checks),file.path(a,'audit/independent_rubin_checks.csv'),row.names=FALSE)
cat('MANUAL RUBIN CHECKS PASSED',length(checks),'\n')