# Contributing

Use Python 3.11 or newer and an isolated environment:

```bash
python -m pip install ".[dev]"
python -m pytest -q
python -m ruff check src tests tools
python -m ruff format --check src tests tools
python -m build
```

Changes should include an observable example or regression test for the intended
behavior. Prefer analytical fixtures, seeded synthetic spectra and independent
oracles to tests that merely repeat implementation details. Keep the schema,
English/Chinese documentation and supported CLI consistent.

Report an issue with the package/Python version, a small synthetic reproduction,
configuration, observed outcome and expected behavior. Do not attach private
experimental data, credentials or unreleasable operational logs. Runtime `runs/`
outputs are intentionally ignored. Include uncertainty and failed cases when
proposing scientific interpretation or tracking changes.

The legacy modules are preserved source, not the supported API. Migrate a legacy
behavior only with matched assumptions and regression evidence. Dataset rights,
calibration and measured-spectrum validation are separate from software tests.

For completed local batches, `python tools/check_saved_batch.py path/to/batch-report.json` independently checks saved curve/likelihood identities and hashes. Its PASS is artifact consistency, not a physical peak-count accuracy label.
