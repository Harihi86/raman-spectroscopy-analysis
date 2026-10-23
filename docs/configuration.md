# Manifest and output reference

[中文使用指南](usage.zh-CN.md). Bundled machine-readable definitions:
[batch manifest schema](../src/raman_analysis/schemas/batch-manifest.schema.json) and
[tracking settings schema](../src/raman_analysis/schemas/tracking-settings.schema.json).
The runtime uses JSON Schema Draft 2020-12 and rejects unknown fields and nonfinite values.

| Field | Meaning |
|---|---|
| `schema_version` | Required; currently 1 |
| `data_kind` | synthetic / measurement / unspecified |
| `defaults` | Shared single-spectrum input options, baseline, fit and metadata |
| `spectra` | 1–10,000 acquisitions, in acquisition order |
| `spectra[].id` | Unique, case-insensitive, filesystem-safe ASCII ID, up to 64 characters |
| `spectra[].label` | Optional display name, including Unicode |
| `spectra[].series` | Independent sequence; defaults to `default`; never inferred from a filename |
| `spectra[].pressure_gpa` | Optional finite pressure, in GPa; missing/null is retained |
| `spectra[].input.path` | Required per acquisition; relative to the manifest directory |
| `spectra[].input` | Explicit columns, error column, delimiter, missing policy, sheet and cm^-1 units |
| `spectra[].analysis` | Per-acquisition overrides of defaults; lists replace entire lists |
| `tracking` | Algorithm settings; see tracking method |
| `outputs.per_spectrum_figures` | Whether each fit writes PNG/PDF; default true |
| `outputs.tracking_figures` | Whether to plot each tracked series; default true |
| `metadata` | User-defined JSON metadata; no executable configuration |

When the baseline method changes, the replacement baseline object is used in full.
When the method stays the same, explicitly specified fields are merged. Each
spectrum requires exactly one resolved `fit.peaks` list or `fit.auto` object. Switching these modes removes the previous mode settings. Input columns may be names or
nonnegative zero-based indices; `header` defaults to row 0; `header: null` preserves a headerless numeric first row and requires integer column indices.

`batch-report.json` uses `format=raman-analysis-batch`, `schema_version=1`, a program
version and a lifecycle state. It contains per-spectrum outcomes, fit peaks,
provenance, acquisition IDs and tracking settings. Per-spectrum raw intensities,
baselines, fit curves and residuals are saved separately. `spectra.jsonl` is flushed
once per completed acquisition, avoiding quadratic rewriting of an entire batch
report during long runs. `run.json` is atomically replaced as progress changes.

`tracking/tracking.json` contains every assignment, gap, candidate identity and
trajectory summary. `assignments.csv` serializes candidate IDs as JSON in a single
cell. A gap has no measured fitted center. Prediction values live in a separate
column and are not substituted for missing observations.

`track` reads a complete batch report and optionally a replacement tracking
settings object. It records the source report's SHA-256 and does not call the
spectral optimizer. A supplied settings file replaces the saved settings as a
whole; omitted fields then use the documented defaults.

Console exit codes: 0 for configured PASS/explicitly disabled tracking, 1 for
completed results requiring review, 2 for invalid configuration/execution input,
and 130 for an interrupted operation. All JSON and HTML use UTF-8; tabular batch
summaries use UTF-8 with BOM. Runtime reports can contain user data and local
paths in error diagnostics; review them before sharing.

Automatic settings and count-score formulas are documented in [automatic.md](automatic.md). AsLS `resample: true` explicitly supports irregular axes by estimating only the baseline on a same-length uniform grid; the fit uses the original samples.
