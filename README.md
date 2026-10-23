# Raman Spectroscopy Analysis

[![checks](https://github.com/Harihi86/raman-spectroscopy-analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/Harihi86/raman-spectroscopy-analysis/actions/workflows/ci.yml)

[English](README.md) · [简体中文](README.zh-CN.md) · [中文使用指南](docs/usage.zh-CN.md)

**Version 0.2** — a reproducible command-line tool for individual Raman spectra,
manifest-based batches and peak tracking across pressure acquisitions. It extends
the owner-authored SHproject source with configurable, independently tested APIs.
No experiment profile, fixed pressure table or sample-specific window is imported
by the supported pipeline.

## What you can do

- Read explicit columns from CSV, TSV, delimited TXT and XLSX; keep original inputs intact.
- Estimate an AsLS or reference-region polynomial baseline and fit 1–8 Gaussian peaks.
- Process spectra with shared defaults and per-spectrum overrides; save each failure independently.
- Track accepted fitted peaks within independently named series, in acquisition order.
- Preserve repeated/decreasing pressures, missing metadata, peak births/disappearances and gaps.
- Use pressure-linear predictions and global assignment; retain ambiguous identities without forcing labels.
- Open an offline bilingual HTML result index, numerical tables, figures and detailed JSON reports.
- Recompute tracking settings from saved fits without rerunning the spectral optimizer.
- Use Chinese help and console output; CSV summaries include a UTF-8 BOM for spreadsheet compatibility.

The public examples are **synthetic**. Track IDs are algorithmic continuity
labels, not chemical or vibrational-mode identifications.

## Install

Python 3.11 or newer. Install into a virtual environment from this repository:

```bash
python -m venv .venv
# macOS / Linux
source .venv/bin/activate
python -m pip install .
```

On Windows PowerShell, activation is optional:

```powershell
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install .
.venv\Scripts\python.exe -m raman_analysis.cli --lang zh --help
```

Runtime dependencies are declared in `pyproject.toml`. The current validated
runtime versions are recorded in `requirements-demo.txt`; other supported
platforms resolve the declared version ranges.

## Try a complete batch

```bash
raman-analyze --lang zh demo-batch --output runs/first-batch
raman-analyze demo-batch --scenario crossing --output runs/crossing
```

Open `runs/first-batch/analysis/report.html` locally. The first example fits
11 generated spectra in two independent sweeps, including repeated and decreasing
pressures. The crossing example fits six generated spectra and checks continuity
when the frequency ordering of two trajectories reverses. Demonstration fitting
priors are explicit and informed by the synthetic construction; this is not a
benchmark of blind peak discovery.

![Generated synthetic pressure trajectories](docs/crossing-tracks.png)

## Analyze your own batch

A complete [portable synthetic batch](examples/batch) includes its manifest, CSV inputs and truth table.

Create a manifest. Each spectrum has a stable ID, an input path and optional
pressure/series metadata. Paths resolve relative to the manifest file.

```json
{
  "schema_version": 1,
  "data_kind": "measurement",
  "defaults": {
    "input": {"x": "shift_cm1", "y": "intensity", "units": "cm^-1"},
    "baseline": {"method": "none"},
    "fit": {"peaks": [{"center": 350, "sigma": 20}, {"center": 720, "sigma": 30}]}
  },
  "tracking": {
    "prediction": "linear",
    "max_distance_cm1": 20,
    "ambiguity_margin_cm1": 1,
    "max_gap": 1
  },
  "spectra": [
    {"id": "s01", "series": "sample-a", "pressure_gpa": 0, "input": {"path": "data/s01.csv"}},
    {"id": "s02", "series": "sample-a", "pressure_gpa": 1, "input": {"path": "data/s02.csv"}}
  ]
}
```

These numerical seeds and tracking thresholds are illustrative. Choose them for
the actual data and resolution. A peak model is an explicit scientific input;
the program does not select every material's lineshape or peak count automatically.

```bash
raman-analyze validate manifest.json --check-files
raman-analyze batch manifest.json --output runs/batch-01
raman-analyze track runs/batch-01/batch-report.json --settings tracking.json --output runs/retracked
```

Use a new output directory for each run. `--no-figures` keeps the numerical and
HTML outputs while disabling figure generation. Add `--json` for machine-readable
console summaries; progress is written to stderr.

[Manifest and output reference](docs/configuration.md) · [Tracking method](docs/tracking.md)
· [中文常见问题](docs/troubleshooting.zh-CN.md).

## Results and quality states

`report.html` links per-spectrum raw/corrected data, fit parameters, residuals and
diagnostics. `summary.csv` preserves acquisition IDs and pressures;
`tracking/assignments.csv` preserves every accepted, ambiguous or excluded peak
and explicit gaps. Empty pressures are not replaced with invented values.

| State | Interpretation |
|---|---|
| PASS | Configured numerical/assignment criteria passed |
| FAIL | A numerical fit did not pass |
| INCONCLUSIVE | Partial results, unresolved fit diagnostics or assignment ambiguity |
| BLOCKED | An input/configuration prevented analysis |
| NOT_RUN | Tracking was explicitly disabled or a run is not yet evaluated |

A bad file does not discard successful spectra. Batch status never reports PASS
when a requested fit or tracking result remains unresolved. An interruption leaves
`run.json` and flushed `spectra.jsonl` records; an incomplete batch is not accepted
for retracking. Automatic resume is not implemented.

## Single spectra and uncertainty

The original supported single-spectrum interface remains available:

```bash
raman-analyze demo --output runs/single-example
raman-analyze run examples/synthetic/config.json --output runs/single-fit
```

AsLS requires approximately uniform sampling. Fits use positive Gaussians with a
linear residual background. Covariance, FWHM and area errors are local estimates
conditional on a fixed baseline and selected model. Baseline/calibration/model
uncertainty is not included. Failed or unidentifiable fits do not update tracks.
See [spectral methodology](docs/methodology.md).

## Verification and scope

```bash
python -m pip install ".[dev]"
python -m pytest -q
python -m ruff check src tests
python -m ruff format --check src tests
python -m build
```

Tests cover analytical spectra, cross-platform fit failures, manifest validation,
partial batches, interruption records, retracking, independent series, crossings,
missing/duplicate pressures and an exhaustive small assignment oracle. CI runs on
Linux/Windows with Python 3.11/3.12. [Observed validation](docs/validation.md).

Current scope: explicit Gaussian models, header-based tabular inputs and local
serial batches of up to 10,000 entries. Non-Gaussian models, automatic peak/model
selection, calibrated tracking probabilities, bootstrap uncertainty, instrument
formats and authorized measured-spectrum validation remain future work.

## Source, roadmap and references

[Preserved SHproject processing source](legacy/SHproject) is experimental legacy
code, with experiment profiles replaced by public templates. It is not imported
by the new batch/tracking API. Raw experimental spreadsheets, historical figures,
run logs and private profiles are excluded. [Source audit](docs/legacy-audit.md).

[Roadmap](ROADMAP.md) · [Changelog](CHANGELOG.md) · [Contributing](CONTRIBUTING.md)
· [Publication scope](docs/publication.md).

References: [Eilers–Boelens AsLS method report](https://prod-dcd-datasets-public-files-eu-west-1.s3.eu-west-1.amazonaws.com/dd7c1919-302c-4ba0-8f88-8aa61e86bb9d),
[SciPy least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html),
[SciPy linear_sum_assignment](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html),
[jsonschema validation](https://python-jsonschema.readthedocs.io/en/stable/validate/).

MIT licensed original software. The software license does not grant rights to
unreleased measurements or external literature. This repository is independent
of the φ³ theoretical-physics project.
