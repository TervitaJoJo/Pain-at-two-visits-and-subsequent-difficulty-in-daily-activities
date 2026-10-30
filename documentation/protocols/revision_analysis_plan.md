# Author-side revision analysis plan

Locked before viewing new associations. Original primary and secondary analyses remain unchanged; all additions are post-primary-results, author-side-review analyses, not original prespecified analyses.

1. Reconstruct the exact original samples and verify every previously exported variable against the original person-level files; no re-selection.
2. At FIRST pain visit, derive treatment categories: neither oral medication nor insulin (both explicit no); oral only (insulin no, oral yes); insulin with or without oral medication (insulin yes). Other values missing. Audit by cohort/window. Add treatment to the primary adjustment set in both cohorts.
3. HRS only: add first-visit BMI (positive reported value) and duration (interview calendar year minus diagnosis year, nonnegative and no greater than age + 1) alongside treatment. No post-first-visit ELSA BMI substitution. Report definition and missingness.
4. Separate sensitivity: add function3 measured at SECOND pain visit. This conditions on function at follow-up start and may adjust for downstream processes; it does not replace the original model.
5. Reuse original primary model coefficient vectors and complete robust covariance matrices to calculate 11 vs 01. Joint D1 pooled Wald test H0: beta10=0 and beta11-beta01=0, with finite complete-data df=persons-parameters, using mitml implementation. No prediction or causal claim.
6. Report proxy participants/intervals/events by exposure; censor before first proxy outcome and never reconnect later outcomes. Retain original main analysis with proxies.
7. Describe mutually exclusive missing-outcome reasons from available statuses. HRS death and other missingness separated without subtracting unlike denominators. Add ADL-only nominal-time sensitivity for comparison with original nominal-time composite; neither resolves competing risks.
8. New-variable analyses: 50 imputations, 20 iterations, fixed seeds, person-level imputation before interval expansion. Same predictor and outcome-summary strategy as original with additional measured variables. No exposure or outcome imputation. Refit an augmented-imputation base model to document effects of the expanded imputation model. Proxy sensitivity uses the original imputations and is explicitly conditional on the original imputation strategy.
9. Inspect diagnostics, missingness, category support, convergence and all planned results. Report all additions regardless of direction or P value; no outcome/threshold/subgroup search.
