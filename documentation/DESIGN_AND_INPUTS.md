# Design, inputs and execution boundaries

## Cohort structure and timing

HRS, ELSA and SHARE form the main synthesis. CHARLS is retained as measurement-adapted exploratory replication because the two pain measurements have different recall wording. This distinction was documented from measurement availability; nonsignificant estimates remain in the package.

Adults with diabetes and no five-item ADL difficulty at the landmark contribute to incident-ADL models. Two preceding pain reports define 00 (neither visit), 01 (newly reported), 10 (no longer reported), and 11 (both visits). These are observed report patterns, not verified continuous pain trajectories or treatment-induced remission. Separate cross-sectional HRS/ELSA analyses are preserved as supplementary context.

The HRS exposure windows use waves 10/11 or 11/12 with the next two visits as outcomes. ELSA uses waves 5/6 or 6/7. The earliest eligible window is selected once per person. SHARE uses waves 5/6 with waves 7/8 outcomes. The extracted SHARE wave-6/7 candidate window is retained only as a measurement audit; it was rejected before effect modeling because current-smoking routing made this covariate almost entirely unavailable. CHARLS uses waves 1/2 with waves 3/4 outcomes.

## Required authorized input layout

Below the configurable raw data root:

```text
HRS/RAND HRS Data/Longitudinal and Cross-Wave Data Products/randhrs1992_2022v1.dta
HRS/public_survey_data/<YEAR>/h<YY>sta/H<YY>C_R.dct
HRS/public_survey_data/<YEAR>/h<YY>da.zip
ELSA/Raw_data/UKDA-5050-stata/stata/stata13_se/h_elsa_g3.dta
CHARLS/Harmonized_CHARLS/H_CHARLS_D_Data.dta
CHARLS/2011/<health-functioning Stata file>
CHARLS/2013/<health-functioning Stata file>
SHARE/<release folders>/sharew<WAVE>_rel9-0-0_<MODULE>.dta
```

Obtain data from each cohort's official access system under its current permissions. This package distributes no cohort data. Preserve official filenames; the code discovers SHARE files recursively and expects a unique relevant release file. HRS pain extraction uses public core data dictionaries and fixed-width archives for the exposure years. The extraction scripts specify every requested field and reject missing fields or invalid identifier alignment. Metadata dictionaries are reconstructed from the authorized harmonized files at runtime, with no dependency on unrelated research projects.

Main SHARE modules are cv_r, dn, ph, br, mh and gv_health as selected in the original builder. Additional extraction uses historical ph diagnosis records from waves 5, 4, 2 and 1, and wave-7/8 cv_r, ph and xt follow-up/death modules. All historical diagnosis information used is no later than baseline.

## Added SHARE variables

- Diabetes medication: ph011d6, taking diabetes drugs at least weekly, binary. This is not the multi-category treatment variable available in HRS/ELSA.
- BMI: gv_health.bmi, valid range 10 to 70.
- Diabetes duration: age minus diagnosis age ph009_5. Use the wave-5 report first, then earlier waves 4, 2, 1 if missing. Reject impossible duration; imputation draws are constrained to a plausible age-based upper bound.
- Landmark function: ph049d2, ph048d3 and ph048d5. Adjustment is a conditional sensitivity because function may lie after part of the pain exposure history.
- Follow-up proxy: ph054_, 1 means respondent alone, 2 or 3 means any proxy participation. Censor at the first proxy visit and all subsequent follow-up; strict sensitivity also excludes unknown status.
- Death: contemporaneous cv_r.deceased or exit-interview xt009_ year, supplemented by xt008_ month when valid. Do not use subsequently updated vital status as if known at an earlier wave. Known death is absorbing.

The scripts are the authoritative field-level dictionary. SHARE eligible n=4215; medication missing 2, BMI missing 82, duration missing 1175, landmark function missing 0. Earlier diagnosis records reduce duration missingness from approximately 70.7% using wave 5 alone to 27.9%.

## Models and missingness

Primary discrete-time complementary-log-log models use the actual interval length as an offset, participant-cluster covariance, and age, sex, education, smoking, arthritis, cardiovascular disease and depression adjustment; HRS/ELSA include entry window and SHARE includes country. Fifty chained-equation imputations address eligible covariate missingness. Exposure, outcomes and death are not imputed. Do not reconnect follow-up after a missed first outcome visit or known death.

Added SHARE observation weights cumulate across both follow-up intervals. Numerator: exposure, interval and country. Denominator: original covariates, exposure, interval and country. The observation-probability models have no event-time offset. Weight predictions are verified for every at-risk interval, including records with missing outcome times. Truncated (pooled observed-interval 1st/99th percentiles) and untruncated fits are retained. The fitted weights are treated as fixed in outcome-model variance estimates; this remains a sensitivity analysis.

The standardized four-year probabilities use two nominal two-year intervals and each cohort's empirical eligible covariate distribution. Risks are pooled on the logit scale and risk differences on the identity scale. Delta-method uncertainty conditions on the empirical standardization population and interval boundaries. These are model-based NET ADL probabilities under the modeled outcome/censoring assumptions, not observed proportions, death-adjusted cumulative incidences, causal intervention effects or risk prediction performance metrics. Averages are within each cohort; the three cohorts' absolute probabilities are not standardized to one common international population.

The SHARE ADL-or-death composite uses an observed ADL event or known death at the corresponding visit, with nominal intervals. It addresses a different outcome and is not combined with primary ADL-only estimates. ADL onset is interval-observed and harmonized mortality information is incomplete across cohorts. This analysis does not fully resolve competing mortality or informative censoring.

Meta-analysis uses REML and modified Hartung-Knapp inference with scale at least one. Four secondary contrasts undergo Holm adjustment within each model variant. With three main cohorts, between-cohort heterogeneity and prediction intervals are imprecise; retain cohort-specific CIs. Neither secondary nominal significance nor a pooled primary association demonstrates that the first pain report adds predictive utility beyond the second.

## Chronology and implementation

The additional-analysis specification states explicitly that primary and four-cohort results were already known. It is a time-stamped additional specification, not original preregistration. Earlier protocols are included with their original chronology; machine paths alone are anonymized in the public copies. Existing source filenames containing `posthoc` are historical labels and do not establish preregistration status.

An initial observation-weight implementation was superseded after checking prediction lengths revealed an unintended time offset. The frozen code explicitly constructs the offset-free observation formula and requires all at-risk rows to be modeled. The reported weights are from the corrected successful run. No superseded result or repair-only script is used in the runnable pipeline. Earlier bounded-duration corrections are documented in the retained protocol note.

No new mediation analysis is represented here. These analyses quantify and test associations and robustness; they do not identify a pain-to-disability mechanism.
