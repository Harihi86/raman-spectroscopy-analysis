# Observed validation / 验证记录

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

Original measurement datasets and full legacy experiment pipelines remain
**NOT_RUN**. No measured-spectrum validation or universal peak-identity accuracy
is claimed. The retained legacy AsLS routine is differentially checked on a
nonflat synthetic fixture; this is narrower than validating the full old workflow.

The v0.1 single-spectrum example remains a regression test. Its generated
centers were 350 and 720 cm^-1; the observed fitted centers were 350.073450 and
720.168367 cm^-1, RMSE 0.154795 a.u. An early Linux run exposed a false PASS for
flat/affine spectra. The acceptance rule was corrected and the cross-platform
regressions remain in the suite.
