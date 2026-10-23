"""Independently check saved numerical outputs; this does not label true peaks."""

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import voigt_profile
from scipy.stats import cauchy, norm


def require(condition, message):
    if not condition:
        raise ValueError(message)


def likelihood_scores(residual, p, errors):
    n = len(residual)
    rss = float(residual @ residual)
    if errors is not None:
        # Gaussian densities provide an independent oracle for the fitter's
        # likelihood constants and parameter/noise penalties.
        deviance = -2 * float(norm.logpdf(residual, scale=errors).sum()) - float(
            np.log(2 * np.pi * errors**2).sum()
        )
        k = p
    else:
        variance = max(rss / n, np.finfo(float).tiny)
        deviance = -2 * float(norm.logpdf(residual, scale=np.sqrt(variance)).sum()) - n * (
            np.log(2 * np.pi) + 1
        )
        if rss / n < np.finfo(float).tiny:
            deviance = n * np.log(np.finfo(float).tiny)
        k = p + 1
    return {
        "bic": deviance + k * np.log(n),
        "aicc": deviance + 2 * k + 2 * k * (k + 1) / (n - k - 1) if n > k + 1 else None,
    }


def check_score(residual, candidate, count, shape, errors):
    p = (4 if shape == "voigt" else 3) * count + 2
    for criterion, expected in likelihood_scores(residual, p, errors).items():
        if expected is None:
            require(
                candidate[criterion] is None, "AICc must be unavailable with insufficient samples"
            )
        else:
            require(
                np.isclose(candidate[criterion], expected, rtol=1e-9, atol=1e-6),
                f"Selected {criterion} disagrees with saved samples",
            )
    if "mean_parameters" in candidate:
        require(candidate["mean_parameters"] == p, "Incorrect mean-parameter penalty")


def reconstructed_peaks(x, peaks, shape):
    result = np.zeros(len(x))
    for peak in peaks:
        offset = x - peak["center_cm1"]
        if shape == "gaussian":
            value = norm.pdf(offset, scale=peak["sigma_cm1"]) / norm.pdf(0, scale=peak["sigma_cm1"])
        elif shape == "lorentzian":
            value = cauchy.pdf(offset, scale=peak["gamma_cm1"]) / cauchy.pdf(
                0, scale=peak["gamma_cm1"]
            )
        else:
            sigma, gamma = peak["sigma_cm1"], peak["gamma_cm1"]
            value = voigt_profile(offset, sigma, gamma) / voigt_profile(0, sigma, gamma)
        result += peak["amplitude"] * value
    return result


def check_batch(report_path):
    report_path = Path(report_path)
    root = report_path.parent
    raw = report_path.read_bytes()
    report = json.loads(raw)
    require(
        report.get("format") == "raman-analysis-batch" and report.get("state") == "COMPLETE",
        "A complete batch report is required",
    )
    for manifest in root.rglob("artifacts.json"):
        for name, expected in json.loads(manifest.read_text(encoding="utf-8")).items():
            require(
                hashlib.sha256((manifest.parent / name).read_bytes()).hexdigest() == expected,
                f"Artifact hash mismatch: {name}",
            )
    checked = 0
    automatic = 0
    family_checks = 0
    for record in report["spectra"]:
        if record["status"] == "BLOCKED":
            continue
        folder = root / record["output"]
        fitted = pd.read_csv(folder / "fit.csv")
        processed = pd.read_csv(folder / "processed.csv")
        result = json.loads((folder / "report.json").read_text(encoding="utf-8"))
        require(
            record["status"] == result["status"], "Saved fit status disagrees with batch record"
        )
        index = np.searchsorted(processed["shift_cm1"].to_numpy(), fitted["shift_cm1"].to_numpy())
        require(np.all(index < len(processed)), "Fit axis is outside the saved input axis")
        np.testing.assert_allclose(
            processed["shift_cm1"].to_numpy()[index], fitted["shift_cm1"], rtol=1e-12
        )
        corrected = processed["corrected_intensity"].to_numpy()[index]
        scale = max(float(np.max(np.abs(processed["raw_intensity"]))), 1)
        np.testing.assert_allclose(
            processed["raw_intensity"] - processed["estimated_baseline"],
            processed["corrected_intensity"],
            atol=scale * 1e-12,
        )
        residual = fitted["residual"].to_numpy()
        np.testing.assert_allclose(
            corrected - fitted["fitted_corrected"], residual, atol=scale * 1e-12
        )
        components = [c for c in fitted if re.fullmatch(r"(?:gaussian|lorentzian|voigt)_\d+", c)]
        require(
            len(components) == len(result["peaks"]),
            "Component count does not match saved peak parameters",
        )
        curves = fitted[components].sum(axis=1) + fitted["residual_background"]
        np.testing.assert_allclose(curves, fitted["fitted_corrected"], atol=scale * 1e-12)
        require(np.isfinite(fitted.to_numpy()).all(), "Nonfinite numerical output")
        rss = float(residual @ residual)
        require(
            np.isclose(rss, result["rss"], rtol=1e-9, atol=scale**2 * 1e-20),
            "Saved RSS disagrees with residual samples",
        )
        checked += 1
        if "model_selection" not in result:
            continue
        decision = json.loads((folder / "model-selection.json").read_text(encoding="utf-8"))
        count = len(result["peaks"])
        require(
            decision["selected_peak_count"] == count, "Selected count does not match saved peaks"
        )
        require(
            decision["count_status"] == record["model_selection"]["count_status"],
            "Count status disagrees with batch record",
        )
        selected_shape = decision.get("selected_shape", "gaussian")
        errors = fitted["measurement_sigma"].to_numpy() if "measurement_sigma" in fitted else None
        candidate = next(
            c
            for c in decision["candidates"]
            if c["peak_count"] == count and c.get("shape", "gaussian") == selected_shape
        )
        check_score(residual, candidate, count, selected_shape, errors)
        if "shape_comparison" in decision:
            comparison = pd.read_csv(folder / "shape-comparison.csv")
            np.testing.assert_allclose(comparison["shift_cm1"], fitted["shift_cm1"], rtol=1e-12)
            for family in decision["shape_comparison"]:
                shape = family["shape"]
                family_residual = comparison[f"{shape}_residual"].to_numpy()
                family_curve = comparison[f"{shape}_fitted_corrected"].to_numpy()
                np.testing.assert_allclose(
                    corrected - family_curve, family_residual, atol=scale * 1e-12
                )
                count_in_family = family["selected_peak_count"]
                columns = [c for c in comparison if re.fullmatch(rf"{shape}_component_\d+", c)]
                require(len(columns) == count_in_family, "Family component count mismatch")
                np.testing.assert_allclose(
                    comparison[columns].sum(axis=1) + comparison[f"{shape}_residual_background"],
                    family_curve,
                    atol=scale * 1e-12,
                )
                model = reconstructed_peaks(
                    comparison["shift_cm1"].to_numpy(), family["peaks"], shape
                )
                np.testing.assert_allclose(
                    model + comparison[f"{shape}_residual_background"],
                    family_curve,
                    atol=scale * 1e-10,
                )
                require(
                    np.isfinite(family_curve).all() and np.isfinite(family_residual).all(),
                    "Nonfinite family curves",
                )
                require(
                    np.isclose(
                        family_residual @ family_residual,
                        family["rss"],
                        rtol=1e-9,
                        atol=scale**2 * 1e-20,
                    ),
                    "Family RSS mismatch",
                )
                check_score(family_residual, family, count_in_family, shape, errors)
                criterion = decision["criterion"]
                scores = [
                    c[criterion]
                    for c in decision["candidates"]
                    if c.get("shape") in {shape, "background"} and c[criterion] is not None
                ]
                require(
                    family[criterion] <= min(scores) + 1e-7, "Family best is not its smallest score"
                )
                family_checks += 1
        criterion = decision["criterion"]
        require(
            candidate[criterion]
            <= min(c[criterion] for c in decision["candidates"] if c[criterion] is not None) + 1e-7,
            "Selected count is not the smallest recorded score",
        )
        automatic += 1
    return {
        "artifact_consistency": "PASS",
        "source_report_sha256": hashlib.sha256(raw).hexdigest(),
        "records": len(report["spectra"]),
        "numerical_spectra_checked": checked,
        "automatic_score_checks": automatic,
        "family_score_checks": family_checks,
        "fit_outcomes": dict(Counter(s["status"] for s in report["spectra"])),
        "physical_peak_count_accuracy": "NOT_ASSESSED",
        "scope": "Saved artifact hashes, curve identities, residual RSS and selected and per-family Gaussian-error likelihood scores (Gaussian, Lorentzian, exact Voigt peaks); no reference peak labels or scientific-validity certification",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        evidence = check_batch(args.report)
    except (ValueError, OSError, KeyError, AssertionError) as exc:
        parser.exit(2, f"Saved-output validation failed: {exc}\n")
    content = json.dumps(evidence, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(content)
    print(content, end="")
