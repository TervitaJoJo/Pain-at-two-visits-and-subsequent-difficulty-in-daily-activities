suppressPackageStartupMessages(library(mice))
p <- normalizePath('analysis/pain_transition_diabetes_adl_20261029',winslash='/')
for(cohort in c('HRS','ELSA')) {
 models <- readRDS(file.path(p,'data_private',paste0(cohort,'_models.rds')))$primary
 q <- list()
 for(i in seq_along(models)) {
  f <- models[[i]]
  for(pair in list(c('11','00'),c('01','00'),c('10','00'),c('10','11'))) {
   v <- setNames(rep(0,length(f$b)),names(f$b))
   if(pair[1]!='00') v[paste0('exposure',pair[1])]<-1
   if(pair[2]!='00') v[paste0('exposure',pair[2])]<--1
   q[[length(q)+1]] <- data.frame(cohort=cohort,imputation=i,contrast=paste(pair,collapse='_vs_'),Q=sum(v*f$b),U=as.numeric(t(v)%*%f$V%*%v),n=f$n,k=f$k)
  }
 }
 write.csv(do.call(rbind,q),file.path(p,'audit',paste0(cohort,'_primary_unpooled_contrasts.csv')),row.names=FALSE)
 # Audit display uses actual category labels, not factor storage indices.
 imps <- list()
 for(sample in c('longitudinal','cross_sectional')) {
  imp <- readRDS(file.path(p,'data_private',paste0(cohort,'_',sample,'_mice50.rds')))
  for(v in names(imp$method)[imp$method!='']) {
   miss <- is.na(imp$data[[v]])
   values <- unlist(lapply(1:imp$m,function(i) as.numeric(as.character(complete(imp,i)[[v]][miss]))))
   imps[[length(imps)+1]] <- data.frame(cohort=cohort,sample=sample,variable=v,method=imp$method[v],missing=sum(miss),imputed_mean=mean(values),imputed_sd=sd(values),imputed_min=min(values),imputed_max=max(values),logged_events=if(is.null(imp$loggedEvents)) 0 else nrow(imp$loggedEvents))
  }
 }
 write.csv(do.call(rbind,imps),file.path(p,'audit',paste0(cohort,'_imputation_summary.csv')),row.names=FALSE)
}
cat('Exported aggregate Rubin-pooling verification inputs and corrected diagnostic category labels.\n')
