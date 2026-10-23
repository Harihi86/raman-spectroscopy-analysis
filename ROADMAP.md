# Development roadmap

| Priority / status | Objective and expected behavior | Dependencies | Acceptance criteria / verification |
|---|---|---|---|
| P0 / Implemented | Explicit imports and input provenance without experiment profiles | Tabular input schema | Unit/duplicate checks, opt-in invalid rows, immutable raw inputs and parsed-byte hashes |
| P0 / Implemented | Preserve failed or unidentifiable fits | Bounded least squares and diagnostics | Nonconvergence, affine spectra, covariance/bound/residual concerns remain visible |
| P1 / Implemented | Batch manifests and independent pressure-series tracking | Per-spectrum API, JSON Schema, global assignment | Repeated/missing/decreasing pressures, births/gaps/crossings, partial batches and exhaustive assignment oracle |
| P1 / Implemented | Usable bilingual delivery | CLI, result artifacts and documented settings | Chinese help/README/guide, offline HTML index, CSV/JSON, complete synthetic demonstrations and saved-fit retracking |
| P0 / Exercised locally | Measured-spectrum applicability checks | Authorized source workbooks and declared import metadata | Record every input, original-byte hashes, full-axis/windowed outcomes and sensitivity; no physical-count accuracy claim without labels |
| P1 / Implemented | Automatic seeds, family and count comparison | Explicit baseline/window and likelihood assumptions | Include zero peaks, BIC/AICc, multi-start diagnostics, search caps, ambiguous counts and original-sample fits |
| P0 / Planned | Independent measured peak-count reference and full legacy migration | Reviewed labels, compatible legacy environment and data rights | Assess count accuracy/reproducibility; migrate deprecated integration with reference checks |
| P1 / Planned | Validated instrument loaders and data-quality diagnostics | Actual format specifications | Shared/pair layouts explicitly described; cosmic-ray/saturation flags without silent edits |
| P2 / Planned | Baseline/calibration uncertainty and model diagnostics | Calibrated errors and alternative baselines | Bootstrap coverage checked on seeded simulations and authorized references |
| P1 / Implemented; measured applicability reported | Lorentzian/exact Voigt comparison, adaptive arPLS and saved-count stability | Validated profile derivatives, family-specific penalties and unchanged diagnostic controls | Analytical limits, covariance oracles, all-family residual exports and full-cohort baseline sensitivity; independent physical-count accuracy remains unassessed |
| P2 / Planned | Reviewed annotations and safer difficult-trajectory handling | Saved candidate identities and scientific review | Explicit manual anchors, uncertainty assessment and robustness against non-linear shifts; no claimed calibrated probability without evidence |
| P2 / Planned | Resume interrupted batches and scale reporting | Stable journals and input/configuration hashes | Reuse only matching completed artifacts; changed inputs cannot inherit stale success |
| P3 / Experimental | Waterfalls, heatmaps and automated scientific reports | Validated batch metadata and measurement references | Observed/interpolated quantities distinguished; no automatic phase identification |

The 0.4 milestone adds alternative symmetric peak families, an adaptive baseline and explicit count-stability reporting. The next scientific milestone is independently reviewed reference labels, artifact diagnostics and propagation of baseline/model uncertainty, guided by the measured failures.
