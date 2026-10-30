# Exploratory prediction module

Status: post hoc relative to the association study; module specification fixed before fitting these prediction models. Existing association results were already known.

Question: does first-visit pain history add to second-visit pain plus a common clinical model that already includes first-visit function?

Primary horizon is the second subsequent scheduled visit, not a common exact four-year horizon. HRS development predictions are frozen for ELSA/SHARE. AUC differences are M1 minus M0; negative Brier differences favor M1.

| Cohort | AUC M0 | AUC M1 | Change (95% CI) | Brier change (95% CI) |
|---|---:|---:|---|---|
| ELSA | 0.7484 | 0.7391 | -0.0093 (-0.0294, +0.0041) | +0.00035 (-0.00232, +0.00698) |
| SHARE | 0.6890 | 0.6826 | -0.0063 (-0.0178, +0.0003) | +0.00288 (+0.00036, +0.00762) |

External confidence intervals use 500 paired two-sample bootstrap replicates, refitting development preprocessing and models and validation observation models. No failed cohort units. Two observation-model fits produced extreme-probability warnings; both converged. Predefined probability clipping and weight caps were retained. The diagnostic rerun reproduced every aggregate bootstrap result byte for byte.

HRS repeated local cross-validation gave a mean AUC change of +0.00057. Some Brier sensitivities improved slightly in ELSA; the result is no stable added value, not universal worsening or proof of equivalence.

Death and loss are censored. Predictions and weighted scoring target hypothetical net probabilities under conditional independent observation; they are not competing-risk cumulative incidence. Calibration estimates are descriptive and did not recalibrate the frozen predictions.

QA: 31 R numerical/leakage checks and 19 independent Python/hash checks passed. The private prediction CSVs are not distributable and are excluded from the submission package.

Disposition: retain as supplementary boundary analysis with brief main-text disclosure. Do not promote to an independent clinical prediction claim.