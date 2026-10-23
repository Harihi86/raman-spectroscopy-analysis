# Development roadmap

| Priority / status | Objective and expected behavior | Dependencies | Acceptance criteria / verification |
|---|---|---|---|
| P0 / Implemented | Explicit imports and input provenance without experiment profiles | Tabular input schema | Unit/duplicate checks, opt-in invalid rows, immutable raw inputs and parsed-byte hashes |
| P0 / Implemented | Preserve failed or unidentifiable fits | Bounded least squares and diagnostics | Nonconvergence, affine spectra, covariance/bound/residual concerns remain visible |
| P1 / Implemented | Batch manifests and independent pressure-series tracking | Per-spectrum API, JSON Schema, global assignment | Repeated/missing/decreasing pressures, births/gaps/crossings, partial batches and exhaustive assignment oracle |
| P1 / Implemented | Usable bilingual delivery | CLI, result artifacts and documented settings | Chinese help/README/guide, offline HTML index, CSV/JSON, complete synthetic demonstrations and saved-fit retracking |
| P0 / Planned | Validate authorized measured spectra and migrate legacy integration safely | Releasable data and instrument metadata | Reproduce approved references; migrate deprecated NumPy integration with matched regression checks |
| P1 / Planned | Validated instrument loaders and data-quality diagnostics | Actual format specifications | Shared/pair layouts explicitly described; cosmic-ray/saturation flags without silent edits |
| P2 / Planned | Baseline/calibration uncertainty and model diagnostics | Calibrated errors and alternative baselines | Bootstrap coverage checked on seeded simulations and authorized references |
| P2 / Planned | Voigt/Lorentzian models and explicit model selection | Validated line profiles and scoring rules | Analytical limits, identifiability tests and held-out comparisons |
| P2 / Planned | Reviewed annotations and safer difficult-trajectory handling | Saved candidate identities and scientific review | Explicit manual anchors, uncertainty assessment and robustness against non-linear shifts; no claimed calibrated probability without evidence |
| P2 / Planned | Resume interrupted batches and scale reporting | Stable journals and input/configuration hashes | Reuse only matching completed artifacts; changed inputs cannot inherit stale success |
| P3 / Experimental | Waterfalls, heatmaps and automated scientific reports | Validated batch metadata and measurement references | Observed/interpolated quantities distinguished; no automatic phase identification |

The 0.2 milestone completes the first batch/tracking workflow and Chinese user
materials. The next scientific milestone is reference validation on authorized
measurements, with model/uncertainty work guided by observed limitations.
