"""Automatic Gaussian proposals and explicit, conditional peak-count selection.

Detection uses smoothed data only for initialization. Every candidate model is
fitted to the same original samples. Information criteria are rankings under an
independent Gaussian error model, not probabilities or physical mode labels.
"""

from dataclasses import asdict, dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

from .fitting import FWHM_FACTOR, fit_peaks


@dataclass(frozen=True)
class AutoSettings:
    max_peaks: int = 5
    criterion: str = "bic"
    delta_threshold: float = 6.0
    prominence_snr: float = 4.0
    smooth_sigma_points: float = 1.5
    min_sigma_cm1: float | None = None
    max_sigma_cm1: float | None = None
    min_separation_cm1: float | None = None
    multistarts: int = 3

    def __post_init__(self):
        if type(self.max_peaks) is not int or not 1 <= self.max_peaks <= 8:
            raise ValueError("auto.max_peaks must be an integer in [1, 8]")
        if type(self.multistarts) is not int or not 1 <= self.multistarts <= 12:
            raise ValueError("auto.multistarts must be an integer in [1, 12]")
        if self.criterion not in {"bic", "aicc"}:
            raise ValueError("auto.criterion must be bic or aicc")
        for name in (
            "delta_threshold",
            "prominence_snr",
            "smooth_sigma_points",
            "min_sigma_cm1",
            "max_sigma_cm1",
            "min_separation_cm1",
        ):
            value = getattr(self, name)
            if value is None and name.endswith("cm1"):
                continue
            if type(value) not in (int, float) or not np.isfinite(value) or value < 0:
                raise ValueError(f"auto.{name} must be a nonnegative finite number")
            if name in {"prominence_snr", "min_sigma_cm1", "max_sigma_cm1"} and value == 0:
                raise ValueError(f"auto.{name} must be positive")


def information_criteria(residuals, parameter_count, errors=None):
    """Drop common likelihood constants; include estimated variance in k."""
    r = np.asarray(residuals, dtype=float)
    n = len(r)
    k = parameter_count + (errors is None)
    if errors is None:
        mse = max(float(np.mean(r**2)), np.finfo(float).tiny)
        deviance = n * np.log(mse)
    else:
        deviance = float(np.dot(r / errors, r / errors))
    return {
        "parameters_including_noise": int(k),
        "bic": float(deviance + k * np.log(n)),
        "aicc": float(deviance + 2 * k + 2 * k * (k + 1) / (n - k - 1)) if n > k + 1 else None,
    }


def _background(x, y, errors, correlation_limit):
    centered = x - np.mean(x)
    design = np.column_stack([np.ones(len(x)), centered])
    weights = np.ones(len(x)) if errors is None else 1 / errors
    parameters = np.linalg.lstsq(design * weights[:, None], y * weights, rcond=None)[0]
    fitted = design @ parameters
    residuals = y - fitted
    rss = float(residuals @ residuals)
    dof = len(x) - 2
    lag1, diagnostics = None, []
    if np.std(residuals) > 1e-8 * max(float(np.max(np.abs(y))), 1):
        if np.std(residuals[:-1]) > 0 and np.std(residuals[1:]) > 0:
            lag1 = float(np.corrcoef(residuals[:-1], residuals[1:])[0, 1])
            if correlation_limit is not None and abs(lag1) > correlation_limit:
                diagnostics.append(
                    "lag-1 residual correlation exceeds the configured heuristic limit"
                )
    chi = float(np.sum((residuals / errors) ** 2) / dof) if errors is not None else None
    if chi is not None and chi > 4:
        diagnostics.append(
            "residuals exceed declared measurement noise (reduced chi-square > 4 heuristic)"
        )
    return {
        "status": "INCONCLUSIVE" if diagnostics else "PASS",
        "diagnostics": diagnostics,
        "peaks": [],
        "rss": rss,
        "rmse": float(np.sqrt(rss / len(x))),
        "residual_lag1_correlation": lag1,
        "residual_correlation_limit": correlation_limit,
        "reduced_chi_square": chi,
        "degrees_of_freedom": dof,
        "covariance": None,
        "covariance_mode": "no peak parameters",
        "uncertainty_scope": "conditional on a fixed baseline and background model",
        "optimizer_message": "weighted linear background fit; zero peaks",
        "evaluations": 1,
        "fitted": fitted,
        "residuals": residuals,
        "components": [],
        "residual_background": fitted,
    }


def _clip_sigma(value, minimum, maximum, margin):
    low, high = minimum * (1 + margin), maximum * (1 - margin)
    if low >= high:
        low = minimum + (maximum - minimum) / 3
        high = maximum - (maximum - minimum) / 3
    if not minimum < low <= high < maximum:
        raise ValueError("Sigma bounds need a representable interior")
    return float(np.clip(value, low, high))


def _proposals(x, signal, noise, config, minimum, maximum):
    smooth = (
        gaussian_filter1d(signal, config.smooth_sigma_points)
        if config.smooth_sigma_points
        else signal.copy()
    )
    indices, properties = find_peaks(
        smooth, prominence=config.prominence_snr * noise, width=(None, None)
    )
    order = sorted(range(len(indices)), key=lambda i: properties["prominences"][i], reverse=True)
    step = float(np.median(np.diff(x)))
    separation = config.min_separation_cm1 if config.min_separation_cm1 is not None else 2 * step
    proposals = []
    sample_indices = np.arange(len(x))
    for i in order:
        center = float(x[indices[i]])
        if any(abs(center - p["center"]) < separation for p in proposals):
            continue
        width = float(
            np.interp(properties["right_ips"][i], sample_indices, x)
            - np.interp(properties["left_ips"][i], sample_indices, x)
        )
        proposals.append(
            {
                "center": center,
                "sigma": _clip_sigma(width / FWHM_FACTOR, minimum, maximum, 0.1),
                "amplitude": max(float(properties["prominences"][i]), noise),
                "prominence": float(properties["prominences"][i]),
            }
        )
    return proposals[:32], smooth, separation


def _seeds(count, start, proposals, previous, x, smooth, noise, minimum, maximum):
    span = float(np.ptp(x))
    if previous and start == 1:
        seeds = [
            {"center": p["center_cm1"], "sigma": p["sigma_cm1"], "amplitude": p["amplitude"]}
            for p in previous["peaks"]
        ]
        residual = gaussian_filter1d(previous["residuals"], 1.5)
        # Grow a previous model at an unexplained positive residual, never use it
        # as a reported peak or skip optimization of the original samples.
        ordering = np.argsort(residual)[::-1]
        index = next(
            (
                int(i)
                for i in ordering
                if 0 < i < len(x) - 1 and all(abs(x[i] - s["center"]) > minimum for s in seeds)
            ),
            len(x) // 2,
        )
        seeds.append(
            {
                "center": float(x[index]),
                "sigma": span / (5 * count),
                "amplitude": max(float(residual[index]), noise),
            }
        )
    elif previous and start >= 2:
        seeds = [
            {"center": p["center_cm1"], "sigma": p["sigma_cm1"], "amplitude": p["amplitude"]}
            for p in previous["peaks"]
        ]
        broadest = max(range(len(seeds)), key=lambda i: seeds[i]["sigma"])
        peak = seeds.pop(broadest)
        shift = peak["sigma"] * (0.45 + 0.1 * (start - 2))
        seeds.extend(
            [
                {**peak, "center": peak["center"] - shift, "amplitude": peak["amplitude"] / 2},
                {**peak, "center": peak["center"] + shift, "amplitude": peak["amplitude"] / 2},
            ]
        )
    else:
        proposals = proposals[start:] + proposals[:start]
        seeds = [{k: p[k] for k in ("center", "sigma", "amplitude")} for p in proposals[:count]]
        for center in np.linspace(x[0], x[-1], 2 * count + 3)[1:-1]:
            if len(seeds) >= count:
                break
            if all(abs(center - p["center"]) > minimum for p in seeds):
                seeds.append(
                    {
                        "center": float(center),
                        "sigma": span / (5 * count),
                        "amplitude": max(float(np.interp(center, x, smooth)), noise),
                    }
                )
    if len(seeds) != count:
        return None
    normalized = []
    for i, seed in enumerate(sorted(seeds, key=lambda p: p["center"])):
        center = float(np.clip(seed["center"], x[0] + span * 1e-5, x[-1] - span * 1e-5))
        if normalized and center <= normalized[-1]["center"] + span * 1e-6:
            center = normalized[-1]["center"] + span * 1e-6
        if center >= x[-1]:
            return None
        normalized.append(
            {
                "center": center,
                "sigma": _clip_sigma(seed["sigma"], minimum, maximum, 0.01),
                "amplitude": max(float(seed["amplitude"]), noise),
                "sigma_bounds": [minimum, maximum],
            }
        )
    return normalized


def auto_fit(x, y, *, auto=None, errors=None, max_nfev=3000, residual_correlation_limit=0.35):
    if auto is not None and not isinstance(auto, dict):
        raise ValueError("fit.auto must be an object")
    config = AutoSettings(**(auto or {}))
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if (
        x.ndim != 1
        or x.shape != y.shape
        or len(x) < 8
        or not np.isfinite([x, y]).all()
        or np.any(np.diff(x) <= 0)
    ):
        raise ValueError(
            "Automatic fitting needs at least eight finite, strictly increasing samples"
        )
    if type(max_nfev) is not int or not 1 <= max_nfev <= 100000:
        raise ValueError("max_nfev must be an integer in [1, 100000]")
    if residual_correlation_limit is not None and not 0 < residual_correlation_limit < 1:
        raise ValueError("residual_correlation_limit must be in (0, 1) or null")
    if errors is not None:
        errors = np.asarray(errors, dtype=float)
        if errors.shape != x.shape or not np.isfinite(errors).all() or np.any(errors <= 0):
            raise ValueError("Measurement errors must be positive and match the data")
    span, step = float(np.ptp(x)), float(np.median(np.diff(x)))
    minimum = config.min_sigma_cm1 or step
    maximum = config.max_sigma_cm1 or span / 3
    if not 0 < minimum < maximum:
        raise ValueError("Automatic sigma bounds must satisfy 0 < minimum < maximum")
    background = _background(x, y, errors, residual_correlation_limit)
    signal = background["residuals"]
    difference = np.diff(y)
    noise = float(
        np.median(np.abs(difference - np.median(difference))) / (0.6744897501960817 * np.sqrt(2))
    )
    numerical_floor = max(
        float(np.max(np.abs(y))) * np.finfo(float).eps * 32, np.sqrt(np.finfo(float).tiny)
    )
    noise = max(noise, numerical_floor)
    proposals, smooth, separation = _proposals(x, signal, noise, config, minimum, maximum)
    affine = bool(
        np.max(np.abs(signal)) <= 1e-12 * max(float(np.max(np.abs(y))), np.finfo(float).tiny)
    )
    fits = {0: background}
    candidates = [
        {
            "peak_count": 0,
            "fit_status": background["status"],
            "attempts": [],
            **information_criteria(signal, 2, errors),
        }
    ]
    incomplete = False
    previous = None
    for count in range(1, config.max_peaks + 1):
        attempts = []
        if not affine and len(x) > 3 * count + 4:
            for start in range(config.multistarts):
                seeds = _seeds(
                    count, start, proposals, previous, x, smooth, noise, minimum, maximum
                )
                if seeds is None:
                    continue
                try:
                    fitted = fit_peaks(
                        x,
                        y,
                        seeds,
                        errors=errors,
                        max_nfev=max_nfev,
                        residual_correlation_limit=residual_correlation_limit,
                    )
                except (ValueError, np.linalg.LinAlgError) as exc:
                    attempts.append(
                        {
                            "start": start,
                            "status": "BLOCKED",
                            "diagnostics": [str(exc)],
                            "initial_peaks": seeds,
                        }
                    )
                    continue
                attempts.append(
                    {
                        "start": start,
                        "status": fitted["status"],
                        "rss": fitted["rss"],
                        "diagnostics": fitted["diagnostics"],
                        "initial_peaks": seeds,
                    }
                )
                if fitted["status"] != "FAIL" and np.isfinite(fitted["rss"]):
                    score = information_criteria(fitted["residuals"], 3 * count + 2, errors)[
                        config.criterion
                    ]
                    previous_score = (
                        information_criteria(fits[count]["residuals"], 3 * count + 2, errors)[
                            config.criterion
                        ]
                        if count in fits
                        else None
                    )
                    if score is not None and (previous_score is None or score < previous_score):
                        fits[count] = fitted
        if count in fits:
            fitted = fits[count]
            candidates.append(
                {
                    "peak_count": count,
                    "fit_status": fitted["status"],
                    "attempts": attempts,
                    **information_criteria(fitted["residuals"], 3 * count + 2, errors),
                }
            )
            previous = fitted
        else:
            candidates.append(
                {
                    "peak_count": count,
                    "fit_status": "NOT_RUN" if affine else "FAIL",
                    "bic": None,
                    "aicc": None,
                    "attempts": attempts,
                    "reason": "affine signal has no peak contrast"
                    if affine
                    else "no converged fit or insufficient degrees of freedom",
                }
            )
            incomplete |= not affine
    scored = [c for c in candidates if c[config.criterion] is not None]
    best = min(scored, key=lambda c: (c[config.criterion], c["peak_count"]))
    selected = best["peak_count"]
    plausible = []
    for candidate in candidates:
        score = candidate[config.criterion]
        candidate["delta"] = float(score - best[config.criterion]) if score is not None else None
        if score is not None and candidate["delta"] <= config.delta_threshold:
            plausible.append(candidate["peak_count"])
    result = fits[selected]
    diagnostics = list(result["diagnostics"])
    reasons = []
    if len(plausible) > 1:
        reasons.append(
            "multiple peak counts are within the configured information-criterion margin"
        )
    if selected == config.max_peaks:
        reasons.append("best model reaches max_peaks; larger counts have not been assessed")
    if incomplete:
        reasons.append(
            "candidate search is incomplete; one or more peak counts have no converged fit"
        )
    if any(p["amplitude"] < config.prominence_snr * noise for p in result["peaks"]):
        reasons.append(
            "a fitted component is below the configured noise-relative amplitude threshold"
        )
    centers = sorted(p["center_cm1"] for p in result["peaks"])
    if any(b - a < separation for a, b in zip(centers, centers[1:])):
        reasons.append("selected components are closer than the configured minimum separation")
    if any(
        min(p["center_cm1"] - x[0], x[-1] - p["center_cm1"]) < p["sigma_cm1"]
        for p in result["peaks"]
    ):
        reasons.append("a selected component is truncated by the fit window")
    diagnostics.extend(reasons)
    result["diagnostics"] = diagnostics
    result["status"] = "INCONCLUSIVE" if diagnostics else "PASS"
    result["model_selection"] = {
        "method": "automatic Gaussian count comparison",
        "settings": asdict(config),
        "criterion": config.criterion,
        "selected_peak_count": selected,
        "plausible_peak_counts": plausible,
        "count_status": result["status"],
        "candidates": candidates,
        "detected_proposals": proposals,
        "noise_estimate": noise,
        "noise_estimate_method": "MAD of first differences / (0.67449 sqrt(2)); heuristic, not supplied measurement sigma",
        "resolved_sigma_bounds_cm1": [minimum, maximum],
        "resolved_min_separation_cm1": separation,
        "smoothing_scope": "initial proposals only; all models fit original samples",
        "scope": "conditional rankings under independent Gaussian errors and a fixed baseline; not physical mode identification or calibrated probabilities",
        "search_complete": not incomplete,
    }
    return result
