# Changelog

## 0.3.0 — 2026-10-23

- Added signal-derived Gaussian seeds, zero-to-N count comparison, BIC/AICc tables, deterministic multi-start search and explicit count ambiguity/limits.
- Added `demo-batch --automatic` and an automatic manifest with no generator-supplied peak parameters.
- Added explicit headerless input and opt-in baseline-only interpolation for irregular sampling.
- Added an analytic fitting Jacobian checked against independent central finite differences.
- Exercised the new API on locally authorized real pressure workbooks, with full-axis/windowed comparisons and baseline sensitivity; measured physical peak-count accuracy remains unassessed without independent labels.
- Updated Chinese/English automatic-selection guidance, outputs and measured-data evidence boundaries.

## 0.2.0 — 2026-10-23

- Added versioned, schema-validated batch manifests with shared defaults and per-acquisition overrides.
- Added pressure-aware global peak assignment, explicit independent series, gaps, births and ambiguous identities.
- Added saved-fit retracking, interruption journals, portable HTML/CSV/JSON outputs and pressure figures.
- Added Chinese README, step-by-step guide, troubleshooting, console messages and configuration references.
- Added synthetic independent-sweep and crossing demonstrations, an exhaustive assignment oracle and integration tests.
- Snapshot input bytes once so provenance hashes refer to the bytes actually parsed.
- Console summaries are now human-readable by default; use `--json` for machine output.

## 0.1.0 — 2026-10-23

- Published an independent single-spectrum API, preserved selected SHproject source and synthetic examples.
- Added explicit imports, baseline choices, bounded Gaussian fits, covariance/residual diagnostics and reproducible outputs.
- Corrected a cross-platform false acceptance of affine spectra and retained the regression tests.
