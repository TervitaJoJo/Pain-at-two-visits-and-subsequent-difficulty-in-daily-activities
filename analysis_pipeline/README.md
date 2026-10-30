# Anonymous reproducibility package

Freeze date: 2026-10-30. This snapshot includes the original HRS/ELSA workflow, its treatment/function/proxy sensitivities, SHARE/CHARLS replication, and the additional analyses completed on the freeze date. Main synthesis uses HRS, ELSA and SHARE. CHARLS remains a measurement-adapted exploratory replication; its null-compatible estimates are retained.

This package is a computational supplement. The manuscript has not yet been revised to the 3+1 structure. Original HRS/ELSA figure scripts reproduce the existing figures; newly generated analyses do not replace the submission figures automatically.

## Run

Use Python and R versions recorded in `../documentation/environment.json`. Install the Python requirements in your own environment and install the stated R packages before execution. The runner does not install software or retrieve restricted cohort data.

From this folder:

```text
python run_pipeline.py --stage summary-check --work-dir "<EMPTY_SUMMARY_WORKSPACE>" --rscript "<RSCRIPT_EXECUTABLE>"
python run_pipeline.py --stage all --data-root "<AUTHORIZED_COHORT_DATA_ROOT>" --work-dir "<EMPTY_ANALYSIS_WORKSPACE>" --rscript "<RSCRIPT_EXECUTABLE>"
```

Replace angle-bracket placeholders with local paths. Choose a short work-directory path on Windows to avoid its long-path limit. All generated data, imputations, model objects, figures, logs and reports are written to the external work directory. Keep this frozen code directory unchanged. Do not publish the generated execution workspace: it contains participant-level records and local paths.

`summary-check` copies only frozen aggregate reference inputs into a fresh workspace, repeats all eight four-cohort synthesis variants and ten three-cohort summaries, performs 333 numerical checks, and renders the seven figure scripts. It does not require participant data and does not perform the two checks that require an individual design matrix.

`all` runs: prepare -> primary -> revision -> external -> additional -> synthesis -> verify -> figures. It refits the original and expanded multiple-imputation models with their original seeds and 50 imputations. Expect substantial runtime. Individual `--stage` options resume completed stages in the same runner-created workspace. The runner disables SHARE cache reuse by default. Never use a raw-data directory, manuscript directory or unrelated existing project as the work directory.

## Folder map

- `source/analysis/pain_transition_diabetes_adl_20261029`: HRS/ELSA eligibility, longitudinal/cross-sectional data, primary and sensitivity models, scalar pooling, verification, aggregate tables.
- Its `revision_20261029`: treatment, HRS treatment/BMI/duration, landmark function, proxy censoring, history test, 11 vs 01, nominal-time outcomes and complete-case checks.
- `source/analysis/pain_external_replication_20261030`: raw CHARLS/SHARE extraction, primary MI, function sensitivity, independent model checking and cohort compilation.
- Its `meta_analysis`: four-cohort and leave-one-cohort-out synthesis. The `omit_CHARLS` analysis is the three-main-cohort primary result.
- Its `additional_20261030`: expanded SHARE clinical covariates, cumulative follow-up observation weights, proxy and death-sensitive outcomes, leave-one-country-out, three-cohort secondary synthesis, and standardized net ADL probabilities.
- `source/presentation`: original HRS/ELSA visualization sources and rendering dependencies. Additional absolute-risk and meta plots live in their analysis folders.
- Sibling `frozen_primary`, `frozen_external`, `frozen_additional`: selected aggregate output snapshots used by summary-check.
- Sibling `documentation`: design, input map, validation scope, software versions and SHA-256 manifest.

No participant datasets, fitted model objects, original respondent identifiers, credentials, author details or machine-specific paths are included. Each folder immediately below `code` contains no more than 100 files recursively. There are no files directly under `code`.

Run `python verify_freeze.py` to check every frozen file and enforce the recursive file-count limit. Existing documents in the submission folder are outside this code snapshot and were not modified.
