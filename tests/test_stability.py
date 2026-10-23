import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from raman_analysis.batch import analyze_batch
from raman_analysis.stability import compare_saved_batches


def batches(tmp_path, *, altered=False):
    x = np.linspace(0, 100, 201)
    y = 4 + np.random.default_rng(9).normal(0, 0.05, len(x))
    pd.DataFrame({"shift_cm1": x, "intensity": y}).to_csv(tmp_path / "input.csv", index=False)
    manifest = {
        "schema_version": 1,
        "defaults": {
            "baseline": {"method": "none"},
            "fit": {"auto": {"max_peaks": 2, "shapes": ["gaussian", "lorentzian", "voigt"]}},
        },
        "tracking": {"enabled": False},
        "spectra": [{"id": "s1", "input": {"path": "input.csv"}}],
    }
    paths = []
    for index in range(2):
        if index:
            manifest["defaults"]["baseline"] = {
                "method": "polynomial",
                "degree": 1,
                "regions": [[0, 100]],
            }
        if index and altered:
            manifest["defaults"]["fit"]["auto"]["prominence_snr"] = 2
        path = tmp_path / f"manifest-{index}.json"
        path.write_text(json.dumps(manifest))
        report = analyze_batch(path, tmp_path / f"out-{index}", no_figures=True)
        assert report["state"] == "COMPLETE"
        paths.append(tmp_path / f"out-{index}/batch-report.json")
    return paths


def test_stability_reports_count_agreement_and_auditable_sources(tmp_path):
    paths = batches(tmp_path)
    result = compare_saved_batches(paths, tmp_path / "stability")
    assert result["spectra"] == 1 and result["variant_count"] == 2
    assert result["stability_counts"] == {"SUPPORTED_STABLE": 1}
    assert result["physical_peak_count_accuracy"] == "NOT_ASSESSED"
    assert len(result["variants"][0]["report_sha256"]) == 64
    assert (tmp_path / "stability/stability.csv").read_bytes().startswith(b"\xef\xbb\xbf")
    assert "物理峰数准确" in (tmp_path / "stability/report.html").read_text(encoding="utf-8")
    with pytest.raises(FileExistsError):
        compare_saved_batches(paths, tmp_path / "stability")


def test_stability_rejects_threshold_drift_and_duplicate_reports(tmp_path):
    paths = batches(tmp_path, altered=True)
    with pytest.raises(ValueError, match="settings or fit budgets differ"):
        compare_saved_batches(paths, tmp_path / "stability")
    with pytest.raises(ValueError, match="Repeated reports"):
        compare_saved_batches([paths[0], paths[0]], tmp_path / "stability")


def test_independent_saved_checker_recalculates_every_family(tmp_path):
    paths = batches(tmp_path)
    tool = Path(__file__).resolve().parents[1] / "tools/check_saved_batch.py"
    spec = importlib.util.spec_from_file_location("saved_batch_checker", tool)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    for path in paths:
        result = checker.check_batch(path)
        assert result["artifact_consistency"] == "PASS"
        assert result["automatic_score_checks"] == 1 and result["family_score_checks"] == 3
    csv = paths[0].parent / "spectra/s1/shape-comparison.csv"
    csv.write_text(csv.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="hash mismatch"):
        compare_saved_batches(paths, tmp_path / "stability")


def test_family_checker_covers_nonzero_voigt_and_known_errors(tmp_path):
    from raman_analysis.profiles import voigt

    x = np.linspace(0, 200, 401)
    errors = np.linspace(0.02, 0.05, len(x))
    y = 2 + voigt(x, 12, 80, 5, 4) + np.random.default_rng(49).normal(0, errors)
    pd.DataFrame({"shift_cm1": x, "intensity": y, "errors": errors}).to_csv(
        tmp_path / "input.csv", index=False
    )
    manifest = {
        "schema_version": 1,
        "defaults": {
            "fit": {"auto": {"max_peaks": 2, "shapes": ["gaussian", "lorentzian", "voigt"]}}
        },
        "tracking": {"enabled": False},
        "spectra": [{"id": "v", "input": {"path": "input.csv", "sigma": "errors"}}],
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest))
    result = analyze_batch(path, tmp_path / "batch", no_figures=True)
    assert result["spectra"][0]["model_selection"]["selected_shape"] == "voigt"
    tool = Path(__file__).resolve().parents[1] / "tools/check_saved_batch.py"
    spec = importlib.util.spec_from_file_location("saved_batch_checker", tool)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    assert checker.check_batch(tmp_path / "batch/batch-report.json")["family_score_checks"] == 3


@pytest.mark.parametrize("case", ["count_cap", "family_disagreement", "blocked"])
def test_stability_never_promotes_agreement_or_drops_blocked_records(tmp_path, case):
    from raman_analysis.profiles import gaussian, lorentzian

    x = np.linspace(0, 200, 401)
    shape = lorentzian if case == "family_disagreement" else gaussian
    y = 3 + shape(x, 12, 80, 7) + np.random.default_rng(401).normal(0, 0.03, len(x))
    pd.DataFrame({"shift_cm1": x, "intensity": y}).to_csv(tmp_path / "input.csv", index=False)
    shapes = ["gaussian", "lorentzian", "voigt"] if case == "family_disagreement" else ["gaussian"]
    reports = []
    for i in range(2):
        baseline = (
            {"method": "none"}
            if not i
            else {"method": "polynomial", "degree": 1, "regions": [[0, 200]]}
        )
        manifest = {
            "schema_version": 1,
            "defaults": {
                "baseline": baseline,
                "fit": {
                    "auto": {
                        "max_peaks": 3 if case == "family_disagreement" else 1,
                        "shapes": shapes,
                    }
                },
            },
            "spectra": [
                {"id": "s", "input": {"path": "missing.csv" if case == "blocked" else "input.csv"}}
            ],
            "tracking": {"enabled": False},
        }
        path = tmp_path / f"m{i}.json"
        path.write_text(json.dumps(manifest))
        analyze_batch(path, tmp_path / f"b{i}", no_figures=True)
        reports.append(tmp_path / f"b{i}/batch-report.json")
    result = compare_saved_batches(reports, tmp_path / "comparison")
    expected = {
        "count_cap": "AGREEMENT_ONLY",
        "family_disagreement": "SENSITIVE",
        "blocked": "UNRESOLVED",
    }[case]
    assert result["stability_counts"] == {expected: 1}
    assert result["status"] == "INCONCLUSIVE"
    assert len(result["records"][0]["variants"]) == 2
