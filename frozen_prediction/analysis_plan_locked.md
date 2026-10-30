# Exploratory prediction increment: locked analysis plan

Status: exploratory extension to an already analysed association study. This plan was fixed before running any new prediction model; it is not original-study preregistration. Existing pain-history association results were already known. No manuscript changes are authorized for this module.

## Question and population
Does first-visit self-reported pain history improve held-out prediction beyond the second-visit report and available clinical predictors among adults with diabetes and no five-item basic ADL difficulty at either exposure visit? Use the existing HRS (development), ELSA and SHARE (external testing) eligible samples without changing the inclusion criteria. CHARLS remains exploratory in the existing association paper and is not used for model selection here.

## Target and time axis
Predict first observed ADL difficulty by the second subsequent scheduled survey visit. Use two discrete visit intervals with a separate interval baseline and complementary log-log link. Do not use an exact four-year label or assign survey observations as exact onset times. The first subsequent visit is a fixed secondary horizon. Stop at the first missing outcome. Death and missing visits censor ascertainment; the target is an attrition-adjusted, hypothetical net probability under conditional independent observation, not a competing-risk cumulative incidence or the unconditional real-world probability of living with disability. The suitability of this target for clinical use is limited. Record actual time distributions by cohort; the later SHARE second visit limits horizon transportability.

## Locked predictors and model comparison
Use identical persons, validation folds and outcome definitions for every pair.
- C (context only): age represented by natural splines with fixed internal knots 60 and 75 and boundary knots 50 and 100; sex; three-level education; current smoking; arthritis; cardiovascular disease; elevated depressive symptoms; first-visit function3 count.
- M0 (comparator): C + second-visit pain.
- M1 (extended): M0 + first-visit pain + first-by-second pain interaction (two added degrees of freedom; equivalent to four pain patterns).
All baseline clinical covariates originate at the first exposure visit and are available at the landmark. No post-landmark measurements are predictors. No variable selection, machine-learning search, data-driven threshold tuning or coefficient refitting in external validation.
Elevated symptoms: CES-D8 >=3 in HRS/ELSA, EURO-D >=4 in SHARE. This pragmatic screen is not proof of psychometric equivalence. Raw-score sensitivity uses each cohort's own depression scale and its entry-wave (HRS/ELSA) or country (SHARE) terms, with local cross-validation only. Do not apply raw CES-D coefficients to EURO-D.

## Missing predictors and validation
Use a deployable train-only imputer: mode for categorical/binary predictors and median for count/continuous predictors. External applications freeze HRS fill values as well as coefficients. These sparse missing predictors do not justify reusing the original outcome-informed MICE objects. No outcomes, validation statistics or future visits enter risk-model preprocessing. Complete-predictor cases form a prespecified sensitivity analysis; do not impute pain or outcomes.
Fit unweighted observed-interval discrete hazards. Validate HRS internally by five repetitions of fivefold participant-level cross-validation; use fixed seeds and pair models within folds. Repeat the same local validation in ELSA and SHARE to assess reproducibility, clearly distinguished from frozen external validation. The primary external analysis trains on all HRS and freezes imputation and all coefficients for ELSA and SHARE. Never recalibrate before reporting frozen performance.

## Observation adjustment and metrics
For each cohort separately, estimate probabilities of observing visit 1 and of observing visit 2 among observed event-free visit-1 participants using fivefold cross-fitted logistic models of the common predictors and the full pain pattern. Nuisance imputation is fitted inside each observation-training fold. The evaluation nuisance process may be estimated in the validation cohort; it never changes any risk prediction. Use the first-visit observation probability for early cases and its product with second-visit observation probability for late cases and controls. Clip each observation probability to [0.02,0.98], cap combined inverse weights at 20, and report clipping, effective sample size, and untruncated-weight sensitivity. This requires observation conditional on the measured predictors, an unverifiable assumption especially for death and SHARE attrition.
At both horizons compute IPCW cumulative/dynamic AUC; paired delta AUC (M1-M0); IPCW Brier score (Horvitz-Thompson denominator = all eligible persons), paired delta Brier (M1-M0, negative better), and normalized/Hajek Brier sensitivity; calibration intercept (ideal 0), slope (ideal 1), weighted observed proportion, mean predicted probability. Use marginal observation weights and unweighted observable-case scores as sensitivity analyses, not replacements selected by results. Calibration fits diagnose predictions and never update them. No NRI/IDI or decision curves without a defined clinical decision threshold.

## Uncertainty and repetitions
External frozen-model validation: 1000 paired participant bootstrap samples conditional on frozen predictions/estimated observation weights; label this conditioning explicitly. Additionally 500 two-sample bootstraps resample HRS development participants and each validation cohort, refit train-only imputation, risk models and cross-fitted validation observation models; report 95% percentile CIs as the preferred uncertainty for external delta metrics. Risk model and evaluation nuisance fits are repeated on resamples; original participant fold assignment groups duplicate rows within observation cross-fitting.
Internal validation: report the mean metric over five OOF prediction sets and its range across repeats. Conditional participant-bootstrap intervals (1000 samples on the same individuals across all repeat-specific prediction sets) do not include full model-fitting uncertainty and must be labelled accordingly, never presented as full nested-resampling confidence intervals. No repeat selection.

## Prespecified sensitivity analyses
1. First-follow-up visit prediction using the same two-interval fitted model.
2. Marginal IPCW, untruncated conditional IPCW and unweighted observable-case scoring.
3. Complete predictor cases, with refitting and external validation.
4. Original continuous depression score and cohort-specific nuisance terms, five-by-five local CV only.
No model without first-visit function is the primary comparator. No country or entry terms are silently refitted for frozen transfer. No additional predictors or alternative cutoffs will be searched to obtain positive results.

## Interpretation and reporting
Evaluate incremental discrimination together with Brier error, calibration, uncertainty and consistency across ELSA and SHARE. A positive training likelihood test is insufficient. There is no validated clinically important delta-AUC cutoff for this application. Tiny or inconsistent improvements, calibration failure, or CIs spanning no increment favour a short supplemental report rather than a prediction-led manuscript. All results, including adverse or null effects, are retained. Current source manuscript/package and frozen association code must remain byte-identical.

## Methodological basis verified before modelling
- Van Calster et al. Calibration: the Achilles heel of predictive analytics. BMC Medicine 2019. DOI 10.1186/s12916-019-1466-7; PMC6912996.
- Suresh et al. Survival prediction models: an introduction to discrete-time modeling. BMC Medical Research Methodology 2022. DOI 10.1186/s12874-022-01679-6; PMC9316420.
- Longitudinal Associations of Perceived Neighborhood Disorder and Lack of Social Cohesion With Depression Among Adults Aged 50 Years or Older. American Journal of Epidemiology 2020;189:343. DOI 10.1093/aje/kwz209. Pragmatic CES-D8 >=3 and EURO-D >=4 thresholds.
- On the estimation of inverse-probability-of-censoring weights for the evaluation of survival prediction error. PMC11785332.
