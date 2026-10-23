import json

import numpy as np
import pandas as pd
import pytest

from raman_analysis.automatic import AutoSettings, auto_fit, information_criteria
from raman_analysis.baseline import estimate_baseline
from raman_analysis.fitting import gaussian
from raman_analysis.io import load_spectrum
from raman_analysis.manifest import merge_settings, validate_manifest
from raman_analysis.pipeline import analyze


def signal(count, *, noise=0.1):
    x = np.linspace(0, 200, 401)
    centers = {0: [], 1: [80], 2: [60, 135], 3: [40, 100, 160]}[count]
    y = 3 + 0.005 * x + sum(gaussian(x, 12 - i * 2, c, 6 + i * 2) for i, c in enumerate(centers))
    return x, y + np.random.default_rng(414 + count).normal(0, noise, len(x)), centers


@pytest.mark.parametrize("count", [0, 1, 2, 3])
def test_automatic_count_and_centers_without_truth_informed_seeds(count):
    x, y, centers = signal(count)
    original = y.copy()
    result = auto_fit(x, y, auto={"max_peaks": 4})
    decision = result["model_selection"]
    assert decision["selected_peak_count"] == count
    assert count in decision["plausible_peak_counts"]
    np.testing.assert_allclose([p["center_cm1"] for p in result["peaks"]], centers, atol=0.15)
    np.testing.assert_array_equal(y, original)
    assert [c["peak_count"] for c in decision["candidates"]] == list(range(5))
    assert decision["smoothing_scope"].startswith("initial proposals only")
    json.dumps(decision, allow_nan=False)


def test_overlapping_pair_is_recovered_even_with_one_initial_local_maximum():
    x = np.linspace(0, 200, 401)
    y = (
        3
        + gaussian(x, 12, 85, 9)
        + gaussian(x, 10, 104, 10)
        + np.random.default_rng(31).normal(0, 0.08, len(x))
    )
    result = auto_fit(x, y, auto={"max_peaks": 3})
    assert result["model_selection"]["selected_peak_count"] == 2
    np.testing.assert_allclose([p["center_cm1"] for p in result["peaks"]], [85, 104], atol=0.6)


def test_criterion_scores_match_likelihood_ratios_and_parameter_penalties():
    r = np.array([1.0, -2.0, 0.5, 3.0, -1.0, 2.0, 0.75, -0.5, 1.2, 2.0, 1.5, 2.5])
    a = information_criteria(r, 2)
    b = information_criteria(r / 2, 5)
    assert a["parameters_including_noise"] == 3
    expected = len(r) * np.log(0.25) + 3 * np.log(len(r))
    assert b["bic"] - a["bic"] == pytest.approx(expected)
    errors = np.linspace(0.5, 2, len(r))
    known = information_criteria(r, 2, errors)
    assert known["parameters_including_noise"] == 2
    assert known["bic"] == pytest.approx(sum((r / errors) ** 2) + 2 * np.log(len(r)))
    assert known["aicc"] == pytest.approx(sum((r / errors) ** 2) + 4 + 12 / (len(r) - 3))
    assert information_criteria(r[:8], 7)["aicc"] is None


def test_known_errors_are_used_for_likelihood_and_deterministic_selection():
    x, _, _ = signal(2)
    errors = np.linspace(0.05, 0.2, len(x))
    y = (
        3
        + 0.005 * x
        + gaussian(x, 12, 60, 6)
        + gaussian(x, 10, 135, 8)
        + np.random.default_rng(9).normal(0, errors)
    )
    result = auto_fit(x, y, errors=errors, auto={"max_peaks": 3, "criterion": "bic"})
    selected = result["model_selection"]["selected_peak_count"]
    assert selected == 2
    candidate = result["model_selection"]["candidates"][selected]
    expected = information_criteria(result["residuals"], 3 * selected + 2, errors)
    assert candidate["bic"] == expected["bic"]
    repeated = auto_fit(x, y, errors=errors, auto={"max_peaks": 3, "criterion": "bic"})
    assert result["model_selection"] == repeated["model_selection"]


def test_near_equal_peak_counts_remain_explicitly_inconclusive():
    x, y, _ = signal(1)
    result = auto_fit(x, y, auto={"max_peaks": 2, "delta_threshold": 1e9})
    assert result["status"] == "INCONCLUSIVE"
    assert result["model_selection"]["plausible_peak_counts"] == [0, 1, 2]
    assert result["model_selection"]["count_status"] == "INCONCLUSIVE"


def test_hitting_count_limit_is_not_claimed_as_a_confident_count():
    x, y, _ = signal(2)
    result = auto_fit(x, y, auto={"max_peaks": 2, "delta_threshold": 0})
    assert len(result["peaks"]) == 2 and result["status"] == "INCONCLUSIVE"
    assert any("max_peaks" in d for d in result["diagnostics"])


def test_failed_candidate_search_does_not_fall_back_to_a_signal_maximum():
    x, y, _ = signal(2)
    result = auto_fit(x, y, auto={"max_peaks": 2}, max_nfev=1)
    assert result["status"] == "INCONCLUSIVE"
    assert result["model_selection"]["search_complete"] is False
    assert result["peaks"] == []
    assert all(c["fit_status"] == "FAIL" for c in result["model_selection"]["candidates"][1:])


@pytest.mark.parametrize("level", [0.0, 1.0, 1e6])
def test_affine_input_reports_zero_supported_peaks(level):
    x = np.linspace(0, 10, 100)
    result = auto_fit(x, level * (2 + 0.2 * x))
    assert result["status"] == "PASS"
    assert result["model_selection"]["selected_peak_count"] == 0
    assert all(c["fit_status"] == "NOT_RUN" for c in result["model_selection"]["candidates"][1:])


@pytest.mark.parametrize(
    "settings",
    [
        {"max_peaks": 0},
        {"max_peaks": True},
        {"criterion": "aic"},
        {"multistarts": 0},
        {"prominence_snr": float("nan")},
        {"delta_threshold": True},
        {"smooth_sigma_points": -1},
    ],
)
def test_invalid_automatic_settings_are_rejected(settings):
    with pytest.raises(ValueError):
        AutoSettings(**settings)


def test_invalid_sampling_and_sigma_bounds_are_rejected():
    x, y, _ = signal(1)
    with pytest.raises(ValueError, match="increasing"):
        auto_fit(x[::-1], y)
    with pytest.raises(ValueError, match="sigma bounds"):
        auto_fit(x, y, auto={"min_sigma_cm1": 20, "max_sigma_cm1": 10})


@pytest.mark.parametrize("extension", ["csv", "xlsx"])
def test_explicit_headerless_input_keeps_the_first_numeric_row(tmp_path, extension):
    x, y, _ = signal(1)
    data = pd.DataFrame({"x": x, "y": y})
    path = tmp_path / f"raw.{extension}"
    if extension == "xlsx":
        data.to_excel(path, index=False, header=False)
    else:
        data.to_csv(path, index=False, header=False)
    original = path.read_bytes()
    loaded = load_spectrum(path, x=0, y=1, header=None)
    assert loaded.provenance["header"] is None
    np.testing.assert_array_equal(loaded.x, x)
    np.testing.assert_allclose(loaded.y, y)
    assert path.read_bytes() == original
    with pytest.raises(ValueError, match="integer column"):
        load_spectrum(path, header=None)


def test_opt_in_irregular_baseline_interpolation_preserves_original_axis():
    x = np.linspace(0, 100, 300) ** 1.03
    true = 2 + 0.01 * x
    y = true + gaussian(x, 5, 50, 4)
    original = x.copy()
    with pytest.raises(ValueError, match="uniform"):
        estimate_baseline(x, y, {"method": "asls"})
    baseline = estimate_baseline(x, y, {"method": "asls", "lam": 1e7, "resample": True})
    assert np.max(np.abs(baseline - true)) < 0.05
    np.testing.assert_array_equal(x, original)


def test_auto_pipeline_exports_count_comparison_and_readable_zero_peak_table(tmp_path):
    x, y, _ = signal(0)
    source = tmp_path / "noise.csv"
    pd.DataFrame({"shift_cm1": x, "intensity": y}).to_csv(source, index=False)
    result = analyze(
        {"input": {"path": "noise.csv"}, "fit": {"auto": {"max_peaks": 2}}},
        tmp_path / "out",
        base=tmp_path,
        save_figures=False,
    )
    assert result["model_selection"]["selected_peak_count"] == 0
    assert pd.read_csv(tmp_path / "out/peaks.csv").empty
    assert len(pd.read_csv(tmp_path / "out/model-selection.csv")) == 3
    assert "自动峰数判断" in (tmp_path / "out/model-selection.html").read_text(encoding="utf-8")
    assert (
        json.loads((tmp_path / "out/model-selection.json").read_text())["selected_peak_count"] == 0
    )


def test_manual_and_automatic_modes_replace_each_other_without_stale_seeds():
    manual = {"fit": {"peaks": [{"center": 50}], "max_nfev": 100}}
    automatic = {"fit": {"auto": {"max_peaks": 3}}}
    assert "peaks" not in merge_settings(manual, automatic)["fit"]
    assert "auto" not in merge_settings(automatic, manual)["fit"]
    manifest = {
        "schema_version": 1,
        "defaults": {"fit": {"auto": {}}},
        "spectra": [{"id": "s1", "input": {"path": "x.csv"}}],
    }
    assert validate_manifest(manifest)[0]["analysis"]["fit"] == {"auto": {}}
    manifest["defaults"]["fit"]["peaks"] = [{"center": 50}]
    with pytest.raises(ValueError):
        validate_manifest(manifest)


@pytest.mark.parametrize("weighted", [False, True])
def test_analytic_model_jacobian_matches_central_finite_differences(weighted):
    from raman_analysis.fitting import model_jacobian

    x = np.linspace(0, 200, 141)
    params = np.array([12.0, 65.0, 7.0, 9.0, 130.0, 12.0, 3.0, 0.01])
    errors = np.linspace(0.05, 0.2, len(x)) if weighted else None

    def model(p):
        result = gaussian(x, *p[:3]) + gaussian(x, *p[3:6]) + p[-2] + p[-1] * (x - x.mean())
        return result if errors is None else result / errors

    numerical = np.empty((len(x), len(params)))
    for i in range(len(params)):
        delta = np.zeros(len(params))
        delta[i] = 1e-5 * max(abs(params[i]), 1)
        numerical[:, i] = (model(params + delta) - model(params - delta)) / (2 * delta[i])
    np.testing.assert_allclose(
        model_jacobian(x, params, 2, errors), numerical, atol=1e-7, rtol=1e-5
    )


def test_automatic_batch_tracks_crossing_without_generated_seeds(tmp_path):
    import importlib.util
    from pathlib import Path

    from raman_analysis.examples import demo_batch

    result = demo_batch(
        tmp_path / "automatic", scenario="crossing", automatic=True, no_figures=True
    )
    assert result["status"] == "PASS" and result["tracks"] == 2
    manifest = json.loads((tmp_path / "automatic/manifest.json").read_text(encoding="utf-8"))
    assert "auto" in manifest["defaults"]["fit"]
    assert all("analysis" not in s for s in manifest["spectra"])
    assert all(s["model_selection"]["selected_peak_count"] == 2 for s in result["spectra"])
    tool_path = Path(__file__).resolve().parents[1] / "tools/check_saved_batch.py"
    spec = importlib.util.spec_from_file_location("saved_batch_checker", tool_path)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    report_path = tmp_path / "automatic/analysis/batch-report.json"
    evidence = checker.check_batch(report_path)
    assert evidence["automatic_score_checks"] == 6
    assert evidence["physical_peak_count_accuracy"] == "NOT_ASSESSED"
    curve = report_path.parent / result["spectra"][0]["output"] / "fit.csv"
    curve.write_text(curve.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="hash mismatch"):
        checker.check_batch(report_path)


def test_aicc_candidates_choose_best_weighted_start_against_independent_oracle(monkeypatch):
    from raman_analysis import automatic

    x, _, _ = signal(2)
    errors = np.linspace(0.05, 0.2, len(x))
    y = (
        3
        + 0.005 * x
        + gaussian(x, 12, 60, 6)
        + gaussian(x, 10, 135, 8)
        + np.random.default_rng(9).normal(0, errors)
    )
    observed = {}
    original = automatic.fit_peaks

    def capture(x, y, peaks, **kwargs):
        fitted = original(x, y, peaks, **kwargs)
        if fitted["status"] != "FAIL":
            k = 3 * len(peaks) + 2
            score = float(
                np.sum((fitted["residuals"] / errors) ** 2)
                + 2 * k
                + 2 * k * (k + 1) / (len(x) - k - 1)
            )
            observed.setdefault(len(peaks), []).append(score)
        return fitted

    monkeypatch.setattr(automatic, "fit_peaks", capture)
    result = automatic.auto_fit(x, y, errors=errors, auto={"criterion": "aicc", "max_peaks": 3})
    candidates = result["model_selection"]["candidates"]
    for candidate in candidates[1:]:
        assert candidate["aicc"] == pytest.approx(min(observed[candidate["peak_count"]]))
    expected = min(candidates, key=lambda c: c["aicc"])["peak_count"]
    assert result["model_selection"]["selected_peak_count"] == expected


def test_narrow_valid_sigma_interval_keeps_initialization_inside_bounds():
    x = np.linspace(0, 200, 401)
    y = 3 + gaussian(x, 12, 80, 6) + np.random.default_rng(44).normal(0, 0.1, len(x))
    result = auto_fit(x, y, auto={"max_peaks": 2, "min_sigma_cm1": 5.99, "max_sigma_cm1": 6.01})
    assert result["model_selection"]["selected_peak_count"] == 1
    assert result["model_selection"]["search_complete"] is True
    assert 5.99 <= result["peaks"][0]["sigma_cm1"] <= 6.01


def test_non_object_automatic_configuration_is_not_treated_as_defaults():
    x, y, _ = signal(1)
    with pytest.raises(ValueError, match="must be an object"):
        auto_fit(x, y, auto=[])
