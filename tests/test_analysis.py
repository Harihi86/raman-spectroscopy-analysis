import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from raman_analysis.baseline import estimate_baseline
from raman_analysis.cli import demo
from raman_analysis.fitting import FWHM_FACTOR, fit_peaks, gaussian
from raman_analysis.io import load_spectrum


def write_data(tmp_path, x=None, y=None, suffix="csv"):
    x = np.linspace(100, 1100, 501) if x is None else x
    y = gaussian(np.asarray(x), 20, 500, 25) if y is None else y
    frame = pd.DataFrame({"shift_cm1": x, "intensity": y})
    path = tmp_path / f"input.{suffix}"
    if suffix == "xlsx":
        frame.to_excel(path, index=False)
    else:
        frame.to_csv(path, index=False)
    return path


@pytest.mark.parametrize("extension", ["csv", "xlsx"])
def test_explicit_import_and_raw_preservation(tmp_path, extension):
    path = write_data(tmp_path, suffix=extension)
    before = path.read_bytes()
    spectrum = load_spectrum(path)
    assert len(spectrum.x) == 501
    assert path.read_bytes() == before
    assert not spectrum.y.flags.writeable


def test_missing_values_require_explicit_drop(tmp_path):
    path = write_data(tmp_path, y=np.r_[np.nan, np.ones(500)])
    with pytest.raises(ValueError, match="missing"):
        load_spectrum(path)
    assert load_spectrum(path, missing="drop").provenance["dropped_rows"] == 1


def test_duplicate_axis_and_wrong_units_rejected(tmp_path):
    path = write_data(tmp_path, x=np.r_[100, np.linspace(100, 1100, 500)])
    with pytest.raises(ValueError, match="Duplicate"):
        load_spectrum(path)
    with pytest.raises(ValueError, match="unit"):
        load_spectrum(path, units="nm")


def test_descending_axis_sorted_and_recorded(tmp_path):
    data = load_spectrum(write_data(tmp_path, x=np.linspace(1100, 100, 501)))
    assert np.all(np.diff(data.x) > 0)
    assert not data.provenance["input_already_sorted"]


def test_asls_flat_input_and_irregular_sampling():
    x = np.linspace(0, 10, 100)
    baseline = estimate_baseline(x, np.zeros(100), {"method": "asls"})
    np.testing.assert_allclose(baseline, 0)
    with pytest.raises(ValueError, match="uniform"):
        estimate_baseline(x**2, np.ones(100), {"method": "asls"})


def test_polynomial_reference_regions():
    x = np.linspace(0, 100, 501)
    y = 2 + 0.02 * x + gaussian(x, 8, 50, 2)
    baseline = estimate_baseline(
        x, y, {"method": "polynomial", "degree": 1, "regions": [[0, 20], [80, 100]]}
    )
    np.testing.assert_allclose(baseline, 2 + 0.02 * x, atol=1e-9)


@pytest.mark.parametrize("centers", [[230], [230, 750], [230, 500, 850]])
def test_arbitrary_peak_counts_and_analytical_widths(centers):
    x = np.linspace(0, 1100, 1101)
    y = (
        sum(gaussian(x, 10 + i, center, 18 + 2 * i) for i, center in enumerate(centers))
        + 1
        + x * 0.001
    )
    peaks = [{"center": c - 2, "sigma": 20} for c in centers]
    result = fit_peaks(x, y, peaks)
    assert result["status"] == "PASS"
    np.testing.assert_allclose([p["center_cm1"] for p in result["peaks"]], centers, atol=1e-5)
    np.testing.assert_allclose(
        [p["fwhm_cm1"] for p in result["peaks"]],
        [(18 + 2 * i) * FWHM_FACTOR for i in range(len(centers))],
        atol=1e-4,
    )
    assert result["rss"] < 1e-10


def test_absolute_measurement_sigma_and_covariance():
    x = np.linspace(0, 100, 501)
    y = gaussian(x, 10, 50, 5) + np.random.default_rng(7).normal(0, 0.1, len(x))
    result = fit_peaks(x, y, [{"center": 49, "sigma": 4}], errors=np.full(len(x), 0.1))
    assert result["status"] == "PASS"
    assert abs(result["peaks"][0]["center_cm1"] - 50) < 0.1
    assert result["peaks"][0]["area_se"] > 0
    assert result["covariance_mode"] == "absolute measurement sigma"


def test_nonconvergence_does_not_fall_back_to_maximum():
    x = np.linspace(0, 100, 101)
    result = fit_peaks(x, gaussian(x, 10, 50, 5), [{"center": 30, "sigma": 20}], max_nfev=1)
    assert result["status"] == "FAIL"
    assert "optimizer did not converge" in result["diagnostics"]


def test_demo_recovers_known_truth_and_exports(tmp_path):
    result = demo(tmp_path / "demo")
    assert result["data_kind"] == "synthetic"
    assert result["status"] == "PASS"
    np.testing.assert_allclose([p["center_cm1"] for p in result["peaks"]], [350, 720], atol=0.5)
    for name in [
        "report.json",
        "processed.csv",
        "fit.csv",
        "peaks.csv",
        "spectrum.png",
        "spectrum.pdf",
        "artifacts.json",
    ]:
        assert (tmp_path / "demo" / "analysis" / name).stat().st_size > 0


def test_preserved_legacy_baseline_agrees_on_nonflat_fixture():
    legacy = Path(__file__).parents[1] / "legacy/SHproject/core_processing.py"
    spec = importlib.util.spec_from_file_location("legacy_processing", legacy)
    module = importlib.util.module_from_spec(spec)
    import sys

    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    x = np.linspace(0, 100, 501)
    y = 1 + x / 100 + gaussian(x, 5, 50, 8)
    expected = module.als_baseline(y, lam=1e5, p=0.01, niter=20)
    actual = estimate_baseline(x, y, {"method": "asls", "lam": 1e5, "p": 0.01, "iterations": 20})
    np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-8)


@pytest.mark.parametrize(
    "bad",
    [
        {"method": "asls", "p": 0},
        {"method": "asls", "lam": -1},
        {"method": "asls", "iterations": 0},
        {"method": "asls", "typo": 1},
    ],
)
def test_invalid_baseline_parameters_fail(bad):
    with pytest.raises(ValueError):
        estimate_baseline(np.arange(20), np.ones(20), bad)


@pytest.mark.parametrize("level", [0, 1, -12, 1e8])
def test_flat_spectrum_covariance_is_inconclusive(level):
    x = np.linspace(0, 100, 201)
    result = fit_peaks(x, np.full_like(x, level), [{"center": 50, "sigma": 10}])
    assert result["status"] == "INCONCLUSIVE"
    assert result["covariance"] is None


def test_affine_input_does_not_establish_a_peak():
    x = np.linspace(0, 100, 201)
    result = fit_peaks(x, 3 + x / 10, [{"center": 50, "sigma": 10}])
    assert result["status"] == "INCONCLUSIVE"
    assert result["covariance"] is None


def test_bound_solution_is_inconclusive():
    x = np.linspace(0, 100, 201)
    result = fit_peaks(
        x, gaussian(x, 10, 70, 5), [{"center": 30, "sigma": 10, "center_bounds": [20, 40]}]
    )
    assert result["status"] == "INCONCLUSIVE"
    assert any("bound" in d for d in result["diagnostics"])


def test_correlated_model_mismatch_remains_inconclusive():
    x = np.linspace(0, 100, 501)
    y = gaussian(x, 10, 50, 5) + 0.2 * np.sin(x / 4)
    result = fit_peaks(x, y, [{"center": 49, "sigma": 4}])
    assert result["status"] == "INCONCLUSIVE"
    assert any("correlation" in d for d in result["diagnostics"])
