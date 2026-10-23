# Observed validation

Date: 2026-10-23. Host: macOS arm64, CPython 3.12.0.

- `python -m pytest -q`: **25 passed**, 0 failed, 0 errors, 0 skipped.
- `python -m ruff check src tests`: PASS.
- `python -m ruff format --check src tests`: PASS.
- Executed demonstration: synthetic two-peak baseline, fit, residual and export pipeline.

`python -m build`: PASS for wheel and source archive. A separate fresh virtual
environment installed the wheel and ran the documented demo successfully.
Dependency audit of the validation environment after updating its installation
tool: no known vulnerabilities reported by pip-audit on 2026-10-23. GitHub
Actions configuration is not itself evidence that a remote run passed; consult
the actual Actions runs for the published commit.

The original measurement datasets and full legacy experiment pipelines are NOT_RUN. The retained legacy AsLS function was differentially checked on a nonflat synthetic fixture.

Synthetic fitted centers: 350.073450, 720.168367 cm^-1; generated centers: 350 and 720 cm^-1. RMSE: 0.154795 a.u. Lag-1 residual correlation: 0.000976. These are software-check observations, not measured experimental results.

An initial Linux CI run exposed a platform-dependent false PASS for a flat spectrum. The implementation now explicitly rejects affine inputs as unidentifiable for Gaussian peaks, with flat inputs at multiple intensity scales and a linear-background regression fixture. The acceptance rule was fixed rather than weakening the test.
