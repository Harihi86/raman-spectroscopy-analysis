import json

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import quad
from scipy.special import voigt_profile

from raman_analysis.automatic import AutoSettings, auto_fit
from raman_analysis.baseline import baseline_with_diagnostics
from raman_analysis.fitting import fit_peaks, model_jacobian
from raman_analysis.pipeline import analyze
from raman_analysis.profiles import (
    FWHM_FACTOR,
    PARAMETERS,
    gaussian,
    lorentzian,
    peak_fwhm,
    peak_summary,
    profile,
    voigt,
)


@pytest.mark.parametrize(
    "shape,params",
    [
        ("gaussian", [7, 80, 5]),
        ("lorentzian", [7, 80, 5]),
        ("voigt", [7, 80, 5, 3]),
        ("voigt", [7, 80, 0.04, 8]),
        ("voigt", [7, 80, 8, 0.04]),
        ("voigt", [7, 80, 0.001, 10]),
        ("voigt", [7, 80, 0.001, 1000]),
    ],
)
@pytest.mark.parametrize("weighted", [False, True])
def test_profile_jacobians_against_central_differences(shape, params, weighted):
    x = np.linspace(0, 200, 421)
    params = np.array([*params, 2, 0.01], dtype=float)
    errors = np.linspace(0.1, 0.8, len(x)) if weighted else None

    def model(p):
        curve = profile(x, p[:-2], shape) + p[-2] + p[-1] * (x - x.mean())
        return curve if errors is None else curve / errors

    numerical = np.empty((len(x), len(params)))
    for i in range(len(params)):
        delta = np.zeros(len(params))
        delta[i] = 1e-5 * max(abs(params[i]), 1)
        numerical[:, i] = (model(params + delta) - model(params - delta)) / (2 * delta[i])
    np.testing.assert_allclose(
        model_jacobian(x, params, 1, errors, shape), numerical, atol=3e-7, rtol=3e-5
    )


@pytest.mark.parametrize(
    "shape,params",
    [
        ("gaussian", [7, 80, 5]),
        ("lorentzian", [7, 80, 5]),
        ("voigt", [7, 80, 5, 3]),
    ],
)
def test_height_area_fwhm_and_covariance_against_independent_oracles(shape, params):
    p = np.array(params, dtype=float)
    covariance = np.eye(len(p)) * 0.01
    covariance[0, 2] = covariance[2, 0] = -0.003
    row = peak_summary(p, shape, covariance)
    assert profile(np.array([p[1]]), p, shape)[0] == pytest.approx(p[0])
    integral = quad(lambda t: float(profile(np.array([t]), p, shape)[0]), -np.inf, np.inf)[0]
    assert row["area"] == pytest.approx(integral, rel=1e-8)
    assert profile(np.array([p[1] + row["fwhm_cm1"] / 2]), p, shape)[0] == pytest.approx(p[0] / 2)
    area_derivative, width_derivative = [], []
    for i in range(len(p)):
        delta = np.zeros(len(p))
        delta[i] = 1e-5 * max(abs(p[i]), 1)
        plus, minus = peak_summary(p + delta, shape), peak_summary(p - delta, shape)
        area_derivative.append((plus["area"] - minus["area"]) / (2 * delta[i]))
        width_derivative.append((plus["fwhm_cm1"] - minus["fwhm_cm1"]) / (2 * delta[i]))
    for key, gradient in [("area_se", area_derivative), ("fwhm_se_cm1", width_derivative)]:
        gradient = np.array(gradient)
        assert row[key] == pytest.approx(np.sqrt(gradient @ covariance @ gradient), rel=1e-7)
    assert (row["sigma_cm1"] is None) == (shape == "lorentzian")
    assert (row["gamma_cm1"] is None) == (shape == "gaussian")


def test_exact_voigt_matches_scipy_and_pure_family_limits():
    x = np.linspace(-80, 80, 301)
    expected = 7 * voigt_profile(x, 5, 3) / voigt_profile(0, 5, 3)
    np.testing.assert_allclose(voigt(x, 7, 0, 5, 3), expected)
    np.testing.assert_allclose(voigt(x, 7, 0, 5, 0), gaussian(x, 7, 0, 5), atol=1e-14)
    np.testing.assert_allclose(voigt(x, 7, 0, 0, 3), lorentzian(x, 7, 0, 3), atol=1e-14)
    assert peak_fwhm([7, 0, 5, 0], "voigt") == pytest.approx(FWHM_FACTOR * 5)
    assert peak_fwhm([7, 0, 0, 3], "voigt") == pytest.approx(6)


@pytest.mark.parametrize(
    "shape,seed",
    [
        ("lorentzian", {"center": 78, "gamma": 6}),
        ("voigt", {"center": 78, "sigma": 6, "gamma": 2}),
    ],
)
def test_fit_recovers_non_gaussian_parameters_and_full_parameter_penalty(shape, seed):
    x = np.linspace(0, 200, 401)
    truth = [12, 80, 5] if shape == "lorentzian" else [12, 80, 5, 3]
    y = 2 + 0.01 * x + profile(x, truth, shape) + np.random.default_rng(68).normal(0, 0.03, len(x))
    fitted = fit_peaks(x, y, [seed], shape=shape, errors=np.full(len(x), 0.03))
    assert fitted["status"] == "PASS"
    assert fitted["mean_parameter_count"] == len(PARAMETERS[shape]) + 2
    assert fitted["peaks"][0]["center_cm1"] == pytest.approx(80, abs=0.03)
    result = auto_fit(x, y, auto={"max_peaks": 2, "shapes": [shape]})
    assert result["model_selection"]["selected_peak_count"] == 1
    best = result["model_selection"]["candidates"][1]
    k = len(PARAMETERS[shape]) + 3  # two background terms + estimated variance
    assert best["parameters_including_noise"] == k
    assert best["bic"] == pytest.approx(
        len(x) * np.log(result["rss"] / len(x)) + k * np.log(len(x))
    )


def test_family_comparison_preserves_disagreement_instead_of_favorable_pass():
    x = np.linspace(0, 200, 501)
    y = 2 + lorentzian(x, 10, 80, 7) + np.random.default_rng(401).normal(0, 0.03, len(x))
    result = auto_fit(x, y, auto={"max_peaks": 3, "shapes": list(PARAMETERS)})
    decision = result["model_selection"]
    assert decision["selected_shape"] == "lorentzian"
    assert decision["selected_peak_count"] == 1
    assert len(decision["candidates"]) == 10  # one shared zero model, 3 * 3 positive models
    assert decision["shape_count_stability"] == "SENSITIVE"
    assert result["status"] == "INCONCLUSIVE"
    assert any("across the compared peak families" in d for d in result["diagnostics"])
    assert len(result["comparison_fits"]) == 3
    json.dumps(decision, allow_nan=False)


@pytest.mark.parametrize(
    "shapes", [[], "voigt", ["pseudo_voigt"], ["gaussian", "gaussian"], [True]]
)
def test_invalid_family_lists_are_rejected(shapes):
    with pytest.raises(ValueError):
        AutoSettings(shapes=shapes)


def test_arpls_avoids_asls_noise_offset_and_reports_nonconvergence():
    x = np.linspace(0, 200, 801)
    truth = 5 + 0.01 * x + 0.0001 * (x - 100) ** 2
    y = truth + gaussian(x, 12, 80, 5) + np.random.default_rng(501).normal(0, 0.3, len(x))
    original = y.copy()
    ar, info = baseline_with_diagnostics(x, y, {"method": "arpls", "lam": 1e6})
    als, _ = baseline_with_diagnostics(x, y, {"method": "asls", "lam": 1e6})
    outside_peak = abs(x - 80) > 25
    assert info["converged"] is True
    assert np.sqrt(np.mean((ar[outside_peak] - truth[outside_peak]) ** 2)) < 0.1
    assert abs(np.mean(ar[outside_peak] - truth[outside_peak])) < abs(
        np.mean(als[outside_peak] - truth[outside_peak])
    )
    _, short = baseline_with_diagnostics(x, y, {"method": "arpls", "lam": 1e6, "iterations": 1})
    assert short["converged"] is False and short["diagnostics"]
    np.testing.assert_array_equal(y, original)


@pytest.mark.parametrize("level", [0, 1e6])
def test_arpls_affine_input_and_explicit_irregular_grid(level):
    x = np.linspace(0, 20, 200) ** 1.03
    y = level * (2 + 0.3 * x)
    baseline, info = baseline_with_diagnostics(x, y, {"method": "arpls", "resample": True})
    np.testing.assert_allclose(baseline, y, atol=max(1, level) * 1e-9)
    assert info["converged"] is True and info["resampled"] is True
    with pytest.raises(ValueError, match="uniform"):
        baseline_with_diagnostics(x, y, {"method": "arpls"})


def test_pipeline_exports_all_family_curves_and_baseline_warnings(tmp_path):
    x = np.linspace(0, 200, 401)
    y = 4 + lorentzian(x, 12, 80, 7) + np.random.default_rng(604).normal(0, 0.1, len(x))
    pd.DataFrame({"shift_cm1": x, "intensity": y}).to_csv(tmp_path / "input.csv", index=False)
    result = analyze(
        {
            "input": {"path": "input.csv"},
            "baseline": {"method": "arpls", "iterations": 1},
            "fit": {"auto": {"max_peaks": 2, "shapes": list(PARAMETERS)}},
        },
        tmp_path / "out",
        base=tmp_path,
        save_figures=False,
    )
    assert result["baseline_diagnostics"]["converged"] is False
    assert result["status"] == "INCONCLUSIVE"
    assert result["model_selection"]["count_status"] == "INCONCLUSIVE"
    curves = pd.read_csv(tmp_path / "out/shape-comparison.csv")
    assert all(f"{shape}_residual" in curves for shape in PARAMETERS)
    assert "Voigt" not in curves.columns
    html = (tmp_path / "out/model-selection.html").read_text(encoding="utf-8")
    assert "lorentzian" in html and "voigt" in html and "Across families" in html


def test_two_voigt_components_use_full_width_blocks_and_input_order():
    x = np.linspace(0, 200, 501)
    errors = np.full(len(x), 0.03)
    y = 3 + 0.01 * x + voigt(x, 12, 65, 4, 2) + voigt(x, 9, 135, 6, 4)
    y += np.random.default_rng(849).normal(0, errors)
    result = fit_peaks(
        x,
        y,
        [{"center": 137, "sigma": 5, "gamma": 3}, {"center": 63, "sigma": 5, "gamma": 3}],
        shape="voigt",
        errors=errors,
    )
    assert result["status"] == "PASS"
    assert result["mean_parameter_count"] == 10
    assert np.shape(result["covariance"]) == (10, 10)
    np.testing.assert_allclose([p["center_cm1"] for p in result["peaks"]], [65, 135], atol=0.04)
    assert result["covariance_parameter_order"].startswith(
        "input peak order: amplitude, center, sigma, gamma"
    )
    np.testing.assert_allclose(
        sum(result["components"]) + result["residual_background"], result["fitted"]
    )
