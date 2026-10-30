# Four-cohort primary-contrast synthesis

Requested on 2026-10-30, after individual-cohort results were available. This is an additional, post hoc synthesis, not an original preregistered analysis.

- Contrast: pain at both exposure visits (11) versus neither visit (00), incident five-item ADL difficulty.
- Retain HRS, ELSA, SHARE and CHARLS, including nonsignificant estimates. One independent estimate per cohort.
- Use original primary_MI pooled log hazard ratio and Rubin total standard error, not rounded or reverse-derived CI values.
- Main: inverse-variance random-effects meta-analysis, REML tau-squared, modified Hartung-Knapp inference (variance scale max[1,q], t with k-1 df). This safeguards against a smaller standard error than conventional random-effects inference with very few cohorts.
- Report conventional Hartung-Knapp and normal-Wald random-effects inference and common-effect inverse-variance pooling as method sensitivities; do not select among them based on significance.
- Report Q, its df and P, REML-based I-squared, tau-squared, cohort weights and a descriptive 95% prediction interval (t k-2 df, tau-squared plus modified-HK mean-effect variance). Prediction interval is approximate and unstable with k=4.
- Sensitivities: exclude CHARLS for pain-recall incompatibility; leave each cohort out once. Same modified-HK REML inference. No funnel asymmetry/Egger test with four cohorts.
- Do not pool secondary pain contrasts in this request. No individual-data models will be refitted. No manuscript changes.
- Country and depression/measurement differences remain; adjusted associations are conceptually rather than perfectly harmonized. Pooled significance cannot establish incremental prediction, mechanism, treatment benefit, or significance in every cohort.
