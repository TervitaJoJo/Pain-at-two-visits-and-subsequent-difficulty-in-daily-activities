# Numerical and leakage checks; invoke with INPUT_ROOT OUTPUT_ROOT qa.
args<-commandArgs(trailingOnly=TRUE)
source(file.path(args[2],'code/01_prediction.R'))
ob<-readRDS(file.path(root,'data_private/analysis_objects.rds'))
checks<-list();add<-function(name,pass,detail=''){checks[[length(checks)+1]]<<-data.frame(check=name,pass=isTRUE(pass),detail=as.character(detail));if(!isTRUE(pass))stop(paste('QA failed:',name,detail))}
set.seed(867)
y<-sample(0:1,50,TRUE);p<-sample(seq(.1,.9,.1),50,TRUE);w<-runif(50,.5,3)
brute<-sum(outer(w[y==1],w[y==0])*outer(p[y==1],p[y==0],function(a,b)(a>b)+.5*(a==b)))/(sum(w[y==1])*sum(w[y==0]))
add('weighted AUC equals pairwise calculation including ties',abs(brute-weighted_auc(y,p,w))<1e-12)
add('perfect weighted AUC',abs(weighted_auc(y,y,w)-1)<1e-12)
add('constant-score AUC equals one half',abs(weighted_auc(y,rep(.5,50),w)-.5)<1e-12)
for(k in names(ob$data)){
 d<-ob$data[[k]];pred<-ob$pred[[k]];pp<-predict_pipeline(ob$fit,d)
 add(paste(k,'stored frozen predictions reproduce'),max(abs(pred-pp))<1e-12)
 d2<-d;d2$y1<-1;d2$y2<-0;d2$t1<-999;d2$t2<-999
 add(paste(k,'future outcome perturbation cannot change predictions'),max(abs(predict_pipeline(ob$fit,d2)-pp))<1e-12)
 revd<-d[nrow(d):1,];revpred<-predict_pipeline(ob$fit,revd)[nrow(d):1,,,drop=FALSE]
 add(paste(k,'predictions do not depend on validation row order'),max(abs(revpred-pp))<1e-12)
 add(paste(k,'every CV repeat covers all persons'),all(is.finite(ob$cvs[[k]]))&&all(is.finite(ob$rawcvs[[k]])))
 folds<-make_folds(d,seed+100+match(k,names(paths)),TRUE);tr<-d[folds!=1,];te<-d[folds==1,];f<-fit_pipeline(tr);pv<-predict_pipeline(f,te)
 add(paste(k,'recomputed held-out fold matches OOF predictions'),max(abs(pv-ob$cvs[[k]][folds==1,,,1]))<1e-12)
 add(paste(k,'first missing visit terminates endpoint observation'),all(is.na(endpoint(d,2)[is.na(d$y1)])))
 add(paste(k,'early events remain cases if later outcome missing'),all(endpoint(d,2)[!is.na(d$y1)&d$y1==1]==1))
 add(paste(k,'cumulative predictions nondecreasing'),all(pp[,,2]>=pp[,,1]))
 # Deidentified row ordinal is internal only, never placed in aggregate outputs.
 ex<-data.frame(row_ordinal=seq_len(nrow(d)),y=endpoint(d,2),weight=get_weights(d,ob$probs[[k]],2),M0=pp[,'M0',2],M1=pp[,'M1',2])
 write.csv(ex,file.path(root,'data_private',paste0(k,'_qa_predictions_PRIVATE.csv')),row.names=FALSE,na='')
}
d<-apply_imputer(ob$data$HRS,ob$fit$imputer);d$pattern<-factor(d$exposure,levels=c('00','01','10','11'));q<-interval_data(d)
ff<-as.formula(paste('y ~ factor(interval) +',common_rhs,'+ pattern'));f<-quiet_glm(ff,q,binomial('cloglog'))
add('history plus interaction equals unrestricted four patterns',max(abs(predict(f,type='response')-predict(ob$fit$fits$M1,type='response')))<1e-10)
# Validate no-pain contextual model is nested within M0 and M0 within M1.
add('M1 adds exactly two coefficients',length(coef(ob$fit$fits$M1))-length(coef(ob$fit$fits$M0))==2)
add('M0 adds exactly one coefficient',length(coef(ob$fit$fits$M0))-length(coef(ob$fit$fits$C))==1)
add('preprocessing field names exclude outcomes and future time',!any(names(ob$fit$imputer)%in%c('y1','y2','t1','t2','death_month','dead1','dead2')))
write.csv(do.call(rbind,checks),file.path(root,'audit/numerical_and_leakage_checks.csv'),row.names=FALSE)
cat('Passed',length(checks),'numerical and leakage checks\n')
