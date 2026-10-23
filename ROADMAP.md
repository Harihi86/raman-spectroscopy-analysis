# Development roadmap

| Priority / status | Objective and expected behavior | Dependencies | Acceptance criteria / verification |
|---|---|---|---|
| P0 / Implemented | Remove implicit layout, pressure and missing-value assumptions from supported API | Explicit tabular schema | CSV/XLSX import, unit checks, duplicate rejection, opt-in missing rows and raw-preservation tests |
| P0 / Implemented | Fail clearly on optimizer uncertainty | Bounded least squares and diagnostics | Nonconvergence remains FAIL; rank/bound concerns remain INCONCLUSIVE; no maximum fallback |
| P0 / Planned | Validate real spectra and legacy migrations | Releasable data and instrument metadata | Curated measurement references reproduce approved fits; migrate deprecated integration calls with regression checks |
| P1 / Implemented | Configurable Gaussian count, baseline, numerical/figure exports | SciPy/NumPy and explicit configuration | Known 1/2/3-peak analytical fixtures and weighted local covariance tests |
| P1 / Planned | Batch manifests and pressure peak tracking | Per-spectrum provenance and explicit acquisition IDs | Missing/duplicate pressures do not reorder identity; synthetic crossings and heating labels tested |
| P1 / Planned | Shared-axis/pair-axis loaders and robust measurement diagnostics | Documented format schemas | Ambiguous files rejected; invalid/cosmic-ray/saturated data flagged without silent modification |
| P2 / Planned | Baseline/calibration uncertainty and model diagnostics | Calibrated errors and alternative baselines | Bootstrap coverage checked on seeded simulations; sensitivity to fit windows, baseline and peak count reported |
| P2 / Planned | Voigt/Lorentzian models and model selection | Validated line profiles and scoring rules | Analytical limiting cases, identifiability checks and held-out synthetic comparisons |
| P2 / Planned | Waterfalls, heatmaps and pressure trends | Batch metadata and tested peak tracking | No unlabelled interpolation or zero fill; observed and interpolated values visibly distinct |
| P3 / Experimental | Automated scientific reports and heating comparisons | Validated pressure series, provenance and domain review | Reports separate measurements, fits and interpretations; no automatic phase identification |

Next milestone: a small batch-manifest API with synthetic pressure-series tests,
followed by authorized measured-spectrum validation. A repository-ready core is
not a claim of universal spectroscopy capability.
