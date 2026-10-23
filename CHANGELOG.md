# Changelog

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
