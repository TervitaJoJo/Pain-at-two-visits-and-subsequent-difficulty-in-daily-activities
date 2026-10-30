# Exploratory prediction extension
This is a separate extension to the frozen association pipeline. It was specified after the association results were known, before fitting the new prediction models. The historical locked plan is retained verbatim, including its original restriction against manuscript editing at that stage. Later authorized manuscript integration does not change that chronology.

## Included and excluded
Included: anonymized source code, locked plan, frozen aggregate outputs, numerical QA, runtime versions and bootstrap-warning details. No participant data, model objects containing data, individual predictions, user paths or source-paper full texts are distributed. The earlier five code projects are unchanged. Main manuscript claims remain about association; this module did not establish added prediction value.

## Aggregate-only reproduction
Requires Python 3.14, NumPy and pandas. From this project directory run:
```
python code/04_summary_only.py .
```
This recomputes paired increments, bootstrap intervals and cross-validation summaries from frozen aggregate CSVs, and checks the bootstrap hash. It does not rerun participant-level validation.

## Full reproduction with authorized cohort data
Use the existing extraction pipeline to recreate the participant tables. Put the three prepared inputs under an ANALYSIS_ROOT with these paths:
- pain_transition_diabetes_adl_20261029/data_private/HRS_longitudinal.csv
- pain_transition_diabetes_adl_20261029/data_private/ELSA_longitudinal.csv
- pain_external_replication_20261030/data_private/SHARE_primary_longitudinal.csv
The original locked input hashes are in audit/plan_lock.json. Required columns are pid (within-cohort unique), exposure (00/01/10/11), age, female, education (1/2/3), smoking, arthritis, cvd, cesd, function3, y1, y2, t1, t2 and entry for HRS/ELSA or country for SHARE. Missing covariates and outcomes retain NA; exposure is never imputed. CES-D8 is used in HRS/ELSA and EURO-D in SHARE. Outcome coding and selection rules are inherited from the extraction pipeline.

Create a private OUTPUT_ROOT with code, results, audit and data_private folders. Copy the supplied code there; copy audit/bootstrap_first_run_sha256.txt there to compare the original frozen bootstrap. R 4.6.0 with jsonlite 2.0.0 and base splines was used. Run sequentially:
```
Rscript OUTPUT_ROOT/code/01_prediction.R ANALYSIS_ROOT OUTPUT_ROOT full
Rscript OUTPUT_ROOT/code/01_prediction.R ANALYSIS_ROOT OUTPUT_ROOT bootstrap
Rscript OUTPUT_ROOT/code/02_validate.R ANALYSIS_ROOT OUTPUT_ROOT qa
python OUTPUT_ROOT/code/03_summarize.py OUTPUT_ROOT
```
The final command checks 18 independent Python metrics plus the original bootstrap hash. Numerical results or byte hashes can differ across runtime versions; inspect any discrepancy rather than silently accepting it. Full reproduction writes protected individual data and predictions under data_private; never submit that folder. It also creates runtime logs in the private working output. None are part of this distributed package.

## Interpretation and checks
M0 contains current (second-visit) pain and common covariates including earlier function. M1 adds first-visit pain and its interaction with current pain. HRS models and imputations are frozen for ELSA/SHARE. Local repeated cross-validation is a distinct analysis. The primary horizon is the second scheduled follow-up, not a shared exact four-year horizon. Death/loss are censored; observation weights target hypothetical net performance under conditional independent observation. Positive delta AUC and negative delta Brier favor M1.
Two bootstrap logistic observation models generated extreme-probability warnings; both converged. Predefined clipping and caps were retained, with no risk-model warning or failed cohort unit identified. Archived QA includes 31 R checks and 19 independent Python/hash checks. Reported external confidence intervals include development and validation resampling (500 two-sample bootstraps). No threshold tuning or NRI/IDI claims are made.
