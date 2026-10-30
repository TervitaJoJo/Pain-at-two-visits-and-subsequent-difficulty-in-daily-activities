Sys.setenv(LANGUAGE='en')
args <- commandArgs(trailingOnly=TRUE)
cohort <- if(length(args)) args[1] else 'HRS'
stopifnot(cohort %in% c('HRS','ELSA'))
original <- normalizePath('analysis/pain_transition_diabetes_adl_20261029', winslash='/', mustWork=TRUE)
rev <- file.path(original,'revision_20261029')
out <- file.path(rev,'results'); audit <- file.path(rev,'audit')
code <- readLines(file.path(original,'code/02_analyze.R'), warn=FALSE)
end <- grep("^cat\\(cohort,'building imputations", code)[1]-1L
eval(parse(text=code[1:end]))

d <- read.csv(file.path(rev,'data_private',paste0(cohort,'_longitudinal.csv')), colClasses=c(exposure='character'))
d$exposure <- factor(d$exposure, levels=c('00','01','10','11'))
d$education <- factor(d$education, levels=1:3)
d$entry <- factor(d$entry, levels=c(2010,2012))
for(v in c('female','smoking','arthritis','cvd')) d[[v]] <- factor(d[[v]], levels=0:1)
d$diabetes_treatment <- factor(d$diabetes_treatment, levels=0:2)
out <- file.path(rev,'results'); audit <- file.path(rev,'audit')
form <- as.formula(paste('y ~',rhs,'+ offset(log(years))'))
all_types <- c('11_vs_00','01_vs_00','10_vs_00','10_vs_11','11_vs_01')
cc_rows <- list(); cc_joint_rows <- list()
add_cc <- function(fit, variant, types=all_types, persons_override=NA_real_) {
  z <- pool_contrasts(list(fit), paste0(variant,'_CC'), types=types)
  z$analysis_source <- 'complete_case'
  cc_rows[[length(cc_rows)+1]] <<- z
}

# Original-variable complete-case set for analyses that do not add new covariates.
cc_orig <- d[complete.cases(d[,covars]),,drop=FALSE]
fit <- fit_glm(intervals(cc_orig), form, binomial('cloglog'), 'primary_history_added_CC', 1)
add_cc(fit, 'primary_history_added', types=c('10_vs_00','11_vs_01'))

# Complete-case checks for each expanded-imputation model.
extra_treatment <- c(covars,'diabetes_treatment')
extra_landmark <- c(covars,'landmark_function3')
cc_treat <- d[complete.cases(d[,extra_treatment]),,drop=FALSE]
cc_landmark <- d[complete.cases(d[,extra_landmark]),,drop=FALSE]
fit <- fit_glm(intervals(cc_treat), update(form,.~.+diabetes_treatment), binomial('cloglog'), 'treatment_adjusted_CC', 1)
add_cc(fit, 'treatment_adjusted')
fit <- fit_glm(intervals(cc_landmark), update(form,.~.+landmark_function3), binomial('cloglog'), 'landmark_function_adjusted_CC', 1)
add_cc(fit, 'landmark_function_adjusted')
if(cohort=='HRS') {
  extra_context <- c(covars,'diabetes_treatment','bmi','diabetes_duration')
  cc_context <- d[complete.cases(d[,extra_context]),,drop=FALSE]
  fit <- fit_glm(intervals(cc_context), update(form,.~.+diabetes_treatment+bmi+diabetes_duration), binomial('cloglog'), 'diabetes_context_CC', 1)
  add_cc(fit, 'diabetes_context')
}

# Sensitivities based on original-variable complete cases.
dd <- cc_orig
bad1 <- !is.na(dd$proxy1) & dd$proxy1==1
bad2 <- bad1 | (!is.na(dd$proxy2) & dd$proxy2==1)
dd$y1[bad1] <- NA; dd$y2[bad2] <- NA
fit <- fit_glm(intervals(dd), form, binomial('cloglog'), 'proxy_censored_CC', 1)
add_cc(fit, 'proxy_censored_original')
if(cohort=='HRS') {
  nd <- intervals(cc_orig); nd$years <- 2
  fit <- fit_glm(nd, form, binomial('cloglog'), 'ADL_only_nominal_CC', 1)
  add_cc(fit, 'ADL_only_nominal_original')
}

# A direct comparison for the expanded-imputation base model: same original CC set.
fit <- fit_glm(intervals(cc_orig), form, binomial('cloglog'), 'augmented_imputation_base_CC', 1)
add_cc(fit, 'augmented_imputation_base')

# Complete-case version of the joint history test: H0 beta10=0 and beta11-beta01=0.
fit <- fit_glm(intervals(cc_orig), form, binomial('cloglog'), 'history_joint_CC', 1)
L <- rbind(contrast_vector(names(fit$b),'10_vs_00'), contrast_vector(names(fit$b),'11_vs_01'))
Q <- as.vector(L %*% fit$b); U <- L %*% fit$V %*% t(L)
Fstat <- as.numeric(t(Q) %*% solve(U) %*% Q / 2)
df2 <- fit$n - fit$k
cc_joint_rows[[1]] <- data.frame(cohort=cohort,hypothesis='beta10=0; beta11-beta01=0',F=Fstat,df1=2,df2=df2,p=pf(Fstat,2,df2,lower.tail=FALSE),persons=fit$n,parameters=fit$k,analysis_source='complete_case')

cc <- do.call(rbind,cc_rows)
write.csv(cc, file.path(out,paste0(cohort,'_complete_case_audit.csv')), row.names=FALSE)
write.csv(do.call(rbind,cc_joint_rows), file.path(out,paste0(cohort,'_history_joint_complete_case.csv')), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(audit,paste0(cohort,'_complete_case_audit_R_session.txt')))
cat(cohort,'complete-case audit complete\n')
