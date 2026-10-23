# Observed validation / 验证记录

Date: 2026-10-23. Version: 0.4.0. Local host: macOS arm64, CPython 3.12.0.

- Full suite: **142 passed**, no failures/errors/skips.
- Ruff check/format for source, tests and independent checking tools: PASS.
- Gaussian/Lorentzian/exact Voigt analytic derivatives match central finite differences with and without measurement weighting, including near-pure limits.
- Profile heights, full-line areas, exact Voigt half-height widths and covariance cross terms are checked against quadrature/numerical-gradient oracles. An unsorted two-Voigt fit checks all width blocks and covariance ordering.
- Automatic comparisons count 3 parameters per Gaussian/Lorentzian peak and 4 per Voigt peak; known measurement errors enter likelihood scores. One shared zero-peak model avoids repeated background candidates.
- Wrong-family fits, count caps, incomplete searches and cross-family disagreement remain visible. Stability tests distinguish supported agreement, agreement with diagnostics, sensitivity and blocked records, and reject changed thresholds or repeated reports.
- Adaptive arPLS tests exercise noise-offset behavior, affine inputs, explicit irregular-axis interpolation and nonconvergence. A finite baseline curve cannot silently turn nonconvergence into PASS.
- A fresh environment installed the wheel from its built artifact (not an editable source path), ran the three-profile synthetic batch, and passed independent checks of all nine family-best curves and scores. `pip check`: PASS.
- All three synthetic inputs select their generating family and one component in the global score comparison; cross-family disagreement remains INCONCLUSIVE. This is a software fixture, not measured accuracy.
- Chinese offline family/count pages and generated comparison figures were visually inspected.

## Measured-spectrum results

- Same 184 declared-window pressure spectra, six predeclared baseline variants (AsLS/arPLS × lambda 10^5/10^6/10^7), three families: 1104 complete spectral analyses and 3312 family-count searches, including all unresolved candidates.
- Existing count cap 5, three starts, 1500 evaluations, score margin 6, proposal SNR 4 and residual-correlation limit 0.35 were retained. Supplied-error chi-square and covariance-condition rules were unchanged.
- Across baseline variants, 163/184 overall best counts changed. Across all baseline/family best counts: 182 SENSITIVE, 2 AGREEMENT_ONLY at the five-peak search cap, 0 SUPPORTED_STABLE. Numerical success does not establish physical mode counts.
- All 1104 saved outputs passed independent numerical checks, including reconstruction of all 3312 family-best curves and likelihood penalties. All 1104 direct raw-column comparisons passed; 12 source workbook hashes were unchanged; no valid pairs were dropped.
- The original Gaussian/AsLS central-lambda control reproduces the previous selected count and RSS for 184/184 spectra (RSS relative/absolute tolerance 10^-10).
- An independent banded-Cholesky arPLS reference checked all 552 curves; maximum error relative to raw intensity scale was 7.70e-9 against a predeclared numerical tolerance of 2e-7. All convergence flags and iteration counts agreed. This checks numerical implementation, not true backgrounds.
- Adaptive baselines can absorb broad spectral structure. Symmetric profiles, residual assumptions and bounded local search still limit interpretation. Physical count accuracy remains NOT_ASSESSED without independently reviewed labels.
- Optional new full-axis family runs were interrupted for runtime and excluded. The primary 184-spectrum windowed cohort and diagnostic controls were unchanged; historical full-axis Gaussian results remain separately versioned.

[Completed measured protocol and aggregate evidence](real-validation.zh-CN.md) · [Machine-readable evidence](measured-validation-summary.json).
The measured record distinguishes unchanged diagnostic controls, baseline/family sensitivity and missing physical reference labels. No public CI access to private measurements is claimed.

The existing [Actions matrix](https://github.com/Harihi86/raman-spectroscopy-analysis/actions/workflows/ci.yml) runs Linux/Windows × Python 3.11/3.12. Check the actual published commit's run; a workflow definition or older green result is narrower evidence than a completed current run.

The records below are retained historical evidence, with their original versions and scopes.

# Historical version 0.3 validation

Date: 2026-10-23. Version: 0.3.0. Local host: macOS arm64, CPython 3.12.0.

- Full suite: **104 passed**, no failures/errors/skips.
- Wheel/source archive built successfully; a fresh virtual environment installed the wheel, ran the automatic crossing demo and manifest validation, and passed saved-score/hash checks and `pip check`.
- Ruff check/format for source, tests and tools: PASS.
- Automatic detection fixtures cover zero to three peaks without truth-informed seeds, overlapping peaks, known heteroscedastic errors, BIC/AICc likelihood formulas, ambiguity, count caps and incomplete searches.
- The analytic Gaussian Jacobian is independently checked by central finite differences, with and without weighting.
- The saved-output checker independently recomputes selected Gaussian log-likelihood scores from exported observations/errors and detects changed artifacts.
- Headerless CSV/XLSX preserve the first numeric row; explicit baseline-only interpolation supports irregular axes without resampling fitted observations.
- A synthetic crossing batch fits and tracks both branches automatically without generator-supplied peak seeds.
- Local measured checks exercise 12 pressure workbooks / 184 distinct real spectra, all in full-axis and declared-window runs, plus 36 prespecified spectra at two alternative baseline lambdas.
- Direct raw-column checks: 368 / 368; 12 source workbook hashes unchanged; no dropped raw pairs.
- Full-axis outcomes: 0 PASS, 184 INCONCLUSIVE. Declared-window outcomes: 6 PASS, 178 INCONCLUSIVE. No BLOCKED inputs.
- Best-count sensitivity: 94/184 across window choices and 23/36 across baseline lambdas. Physical peak-count accuracy is **NOT_ASSESSED** without independent labels.
- Current source exactly replayed three selected measured records, including a configured PASS and two INCONCLUSIVE cases. The full search traces, peak tables and RSS matched the saved results.
- Chinese count-comparison pages, the local measured report and selected measured fit/residual plots were visually checked.

[Detailed measured protocol and limits](real-validation-v03.zh-CN.md) · [Aggregate JSON](measured-validation-v03-summary.json).
These are local real-data observations, not a claim that public CI accessed private measurements.
Build/install evidence and platform results should be checked for the actual published commit in [Actions](https://github.com/Harihi86/raman-spectroscopy-analysis/actions/workflows/ci.yml).

The version 0.2 record below is retained as historical evidence; its counts and limitations do not describe the new measured-data checks.

# Historical version 0.2 validation

Date: 2026-10-23. Version: 0.2.0. Local host: macOS arm64, CPython 3.12.0.
This page records executed checks, with synthetic validation kept separate from measured-spectrum evidence.

## Local checks

- `python -m pytest -q`: **72 passed**, 0 failures/errors/skips.
- `python -m ruff check src tests`: PASS.
- `python -m ruff format --check src tests`: PASS.
- `python -m build`: wheel and source archive built successfully.
- A fresh virtual environment installed the wheel, loaded the bundled schemas, validated the portable batch inputs, ran both batch demonstrations and retracked saved fits.
- `pip check` in the wheel environment: no broken requirements.
- `pip-audit --skip-editable` on the validation environment: no known vulnerabilities reported at the date above. This is a dated advisory lookup, not a security guarantee.
- Output artifact hashes and parsed input-byte provenance: verified.
- Chinese offline HTML and generated crossing plot: visually checked; local Markdown links checked.
- Source distribution contains Chinese documentation, synthetic batch fixtures, tests and preserved legacy processing source. Python source parses under the Python 3.11 syntax target.

## Synthetic end-to-end results

| Fixture | Spectra | Independent series | Adopted tracks | Outcome | Maximum absolute center error |
|---|---:|---:|---:|---|---:|
| Independent sweeps, including repeated/decreasing pressures | 11 | 2 | 4 | PASS | 0.065537 cm^-1 |
| Frequency-order crossing, sampled on both sides | 6 | 1 | 2 | PASS | 0.039662 cm^-1 |

The branch centers are compared with independently saved generator truth.
Fitting priors are explicitly informed by the synthetic construction, so these
results do not measure blind peak discovery or establish physical mode identities.
The crossing diagram's lines do not imply a measured spectrum at the unsampled intersection.

The tests additionally cover exact/near assignment ambiguity, an exhaustive small
assignment oracle, optional width costs, peak birth/disappearance, gaps and
expiry, missing pressure, independent interleaved series and exclusion of
untrusted fits. Batch tests cover relative/Chinese paths, a shared Excel workbook
with explicit columns, exact agreement with independently executed single-spectrum
fits, duplicate IDs, malformed manifests, absent/corrupt inputs, partial outcomes,
interruption journals and final export failures. A failed final export cannot
leave a completed report eligible for retracking.

## Remote CI evidence

The existing workflow runs the full suite, lint/format checks, a single-spectrum
demo and packaging on Linux and Windows with Python 3.11/3.12. Both batch
scenarios execute as integration tests in every matrix job. Consult the
[actual Actions runs](https://github.com/Harihi86/raman-spectroscopy-analysis/actions/workflows/ci.yml)
for the commit being assessed. A workflow definition or an earlier green commit
alone does not establish a pass for later changes.

## Scientific limits and retained history

At the version 0.2 milestone, original measurement datasets and full legacy experiment pipelines were **NOT_RUN**. No measured-spectrum validation or universal peak-identity accuracy
is claimed. The retained legacy AsLS routine is differentially checked on a
nonflat synthetic fixture; this is narrower than validating the full old workflow.

The v0.1 single-spectrum example remains a regression test. Its generated
centers were 350 and 720 cm^-1; the observed fitted centers were 350.073450 and
720.168367 cm^-1, RMSE 0.154795 a.u. An early Linux run exposed a false PASS for
flat/affine spectra. The acceptance rule was corrected and the cross-platform
regressions remain in the suite.
