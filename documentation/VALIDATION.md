# Validation scope

Validation date: 2026-10-30.

The new statistical analyses were completed in the private analysis workspace. There were 1251 SHARE model fits and 150 function-adjusted fits; all converged with no fitted-model warnings. The expanded SHARE MICE run had no logged events. Its final imputed-variable PSRF values ranged from approximately 0.986 to 1.015. This is a numerical check, not proof of the missing-data assumptions.

The private independent verifier passed 335 checks. These include explicit Rubin/Barnard-Rubin scalar pooling, confidence intervals, standardized probabilities and differences, independent REML and modified Hartung-Knapp calculations for ten synthesis rows, specification integrity, and a participant-cluster sandwich covariance/score check for the expanded clinical-context model.

For the anonymous package:

1. All Python sources were syntax-compiled without writing bytecode into the code tree. All 13 R sources were parsed successfully.
2. All four cohorts' extraction steps were executed in a separate workspace using authorized raw data. Twelve generated datasets, including the historical cross-sectional datasets, both HRS/ELSA revision datasets, SHARE candidate/primary windows, CHARLS and expanded SHARE, matched the frozen private originals exactly, including values, missingness and row order.
3. Summary-only execution required no participant data. It reproduced eight four-cohort synthesis variants and ten three-cohort summaries, passed 333 independent checks, and rendered all seven figure scripts. The two private-design checks are intentionally unavailable in this mode; absence of private inputs is not counted as a successful private-data check.
4. All 94 bundled CSV files were reviewed as aggregate estimates, counts or model/imputation diagnostics. None contains participant identifiers, participant-level rows, model matrices, imputed datasets or fitted-model objects. The generic submission audit flags CSV extensions conservatively; its 94 such warnings were reviewed, not silently removed. There were zero remaining automatic errors for personal information, credentials or forbidden identifying text and no absolute-path warnings.
5. Source files are rebased to user-supplied data and execution paths. No scientific model formula or seed was changed for packaging. Adaptation details and source hashes are recorded separately.
6. The portable package was not subjected to a second full refit of every historical 50-imputation model. Full extraction, source parsing, summary recomputation and plotting were repeated. Original fitted results remain linked to their source code and pre-packaging verification outputs. Raw-data execution of all stages is provided for authorized replication.

The original manuscript and submission files are outside the computational supplement. Their before/after SHA-256 comparison is stored in the private audit directory. The new code and results have their own manifest; the earlier submission manifest was not edited.

The snapshot is tamper-evident by SHA-256, not a claim of original preregistration or a cryptographically signed public timestamp. `analysis_pipeline/verify_freeze.py` checks the frozen files and the recursive 100-file limit for every first-level code folder. Running the analysis writes only to a separate execution workspace.
