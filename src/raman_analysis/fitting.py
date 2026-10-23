"""Bounded symmetric peak fits with residual and local covariance diagnostics."""

import numpy as np
from scipy.optimize import least_squares

from .profiles import FWHM_FACTOR, PARAMETERS, gaussian, peak_summary, profile, profile_jacobian

__all__ = ["FWHM_FACTOR", "gaussian", "model_jacobian", "fit_peaks"]


def model_jacobian(x, params, peak_count, errors=None, shape="gaussian"):
    """Analytic derivatives in input-peak order, followed by intercept/slope."""
    size = len(PARAMETERS[shape])
    jac = np.empty((len(x), size * peak_count + 2))
    for i in range(peak_count):
        jac[:, size * i : size * (i + 1)] = profile_jacobian(
            x, params[size * i : size * (i + 1)], shape
        )
    jac[:, -2] = 1
    jac[:, -1] = x - np.mean(x)
    return jac if errors is None else jac / errors[:, None]


def fit_peaks(
    x,
    y,
    peaks,
    *,
    shape="gaussian",
    errors=None,
    max_nfev=10000,
    residual_correlation_limit=0.35,
):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or x.shape != y.shape or len(x) < 8 or not np.isfinite([x, y]).all():
        raise ValueError("Fit arrays must have matching shapes and finite values")
    if not np.all(np.diff(x) > 0):
        raise ValueError("Fit axis must be strictly increasing")
    if not 1 <= len(peaks) <= 8:
        raise ValueError("Specify 1..8 peaks explicitly")
    if shape not in PARAMETERS:
        raise ValueError("shape must be gaussian, lorentzian or voigt")
    size = len(PARAMETERS[shape])
    if type(max_nfev) is not int or not 1 <= max_nfev <= 100000:
        raise ValueError("max_nfev must be in [1, 100000]")
    if residual_correlation_limit is not None and not 0 < residual_correlation_limit < 1:
        raise ValueError("residual_correlation_limit must be in (0, 1) or null")
    if len(x) <= size * len(peaks) + 2:
        raise ValueError("Insufficient residual degrees of freedom")
    if errors is not None:
        errors = np.asarray(errors, dtype=float)
        if errors.shape != x.shape or not np.isfinite(errors).all() or np.any(errors <= 0):
            raise ValueError("Measurement errors must be positive and match the data")
    reference = float(np.mean(x))
    x0 = x - reference
    # A numerical Jacobian can appear full rank at a vanishing peak on some
    # BLAS/optimizer combinations. Exact affine input supplies no peak evidence.
    affine = np.polynomial.Polynomial.fit(x0, y, 1)(x0)
    intensity_scale = max(float(np.max(np.abs(y))), np.finfo(float).tiny)
    no_peak_contrast = bool(np.max(np.abs(y - affine)) <= 1e-12 * intensity_scale)
    step, span = float(np.min(np.diff(x))), float(np.ptp(x))
    initial, lower, upper = [], [], []
    seeds = []
    for peak in peaks:
        allowed = {*PARAMETERS[shape], "center_bounds"}
        allowed.update(f"{name}_bounds" for name in PARAMETERS[shape][2:])
        if set(peak) - allowed:
            raise ValueError(f"Unknown or incompatible {shape} peak configuration field")
        center = float(peak["center"])
        amplitude = float(peak.get("amplitude", max(float(np.ptp(y)), 1e-3)))
        cb = peak.get("center_bounds", [float(x[0]), float(x[-1])])
        if len(cb) != 2 or not np.isfinite([center, amplitude, *cb]).all():
            raise ValueError("Peak parameters and bounds must be finite")
        if not x[0] <= cb[0] < center < cb[1] <= x[-1] or amplitude <= 0:
            raise ValueError("Initial peaks must lie strictly inside positive physical bounds")
        seeds.append(center)
        initial.extend([amplitude, center])
        lower.extend([0, cb[0]])
        upper.extend([np.inf, cb[1]])
        for name in PARAMETERS[shape][2:]:
            width = float(peak.get(name, span / (12 * len(peaks))))
            wb = peak.get(f"{name}_bounds", [step / 4, span])
            if len(wb) != 2 or not np.isfinite([width, *wb]).all():
                raise ValueError("Peak parameters and bounds must be finite")
            if not 0 < wb[0] < width < wb[1]:
                raise ValueError("Initial peaks must lie strictly inside positive physical bounds")
            initial.append(width)
            lower.append(wb[0])
            upper.append(wb[1])
    if len(set(seeds)) != len(seeds):
        raise ValueError("Initial peak centers must be distinct")
    initial.extend([float(np.percentile(y, 10)), 0.0])
    lower.extend([-np.inf, -np.inf])
    upper.extend([np.inf, np.inf])

    def components(params):
        return [profile(x, params[size * i : size * (i + 1)], shape) for i in range(len(peaks))]

    def model(params):
        return np.sum(components(params), axis=0) + params[-2] + params[-1] * x0

    def residual(params):
        r = model(params) - y
        return r if errors is None else r / errors

    result = least_squares(
        residual,
        initial,
        jac=lambda params: model_jacobian(x, params, len(peaks), errors, shape),
        bounds=(lower, upper),
        max_nfev=max_nfev,
        x_scale="jac",
        loss="linear",
    )
    fitted = model(result.x)
    residuals = y - fitted
    rank = int(np.linalg.matrix_rank(result.jac))
    condition = float(np.linalg.cond(result.jac))
    identifiable = (
        rank == len(initial)
        and np.isfinite(condition)
        and condition < 1e10
        and not no_peak_contrast
    )
    dof = len(x) - len(initial)
    covariance = None
    if identifiable:
        _, singular, vt = np.linalg.svd(result.jac, full_matrices=False)
        covariance = (vt.T / singular**2) @ vt
        if errors is None:
            covariance *= float(np.dot(result.fun, result.fun)) / dof
    diagnostics = []
    if not result.success:
        diagnostics.append("optimizer did not converge")
    if rank != len(initial) or not np.isfinite(condition) or condition >= 1e10:
        diagnostics.append("local covariance is rank deficient or poorly conditioned")
    if no_peak_contrast:
        diagnostics.append("flat or affine input contains no identifiable peak")
    if np.any(result.active_mask):
        diagnostics.append("one or more fit parameters lie on a bound")
    residual_rmse = float(np.sqrt(np.mean(residuals**2)))
    lag1 = None
    if residual_rmse > 1e-8 * max(float(np.max(np.abs(y))), 1):
        if np.std(residuals[:-1]) > 0 and np.std(residuals[1:]) > 0:
            lag1 = float(np.corrcoef(residuals[:-1], residuals[1:])[0, 1])
            if residual_correlation_limit is not None and abs(lag1) > residual_correlation_limit:
                diagnostics.append(
                    "lag-1 residual correlation exceeds the configured heuristic limit"
                )
    reduced_chi_square = float(np.dot(result.fun, result.fun) / dof) if errors is not None else None
    if reduced_chi_square is not None and reduced_chi_square > 4:
        diagnostics.append(
            "residuals exceed declared measurement noise (reduced chi-square > 4 heuristic)"
        )
    status = "FAIL" if not result.success else ("INCONCLUSIVE" if diagnostics else "PASS")
    table = []
    for i in range(len(peaks)):
        block = (
            None
            if covariance is None
            else covariance[size * i : size * (i + 1), size * i : size * (i + 1)]
        )
        table.append(peak_summary(result.x[size * i : size * (i + 1)], shape, block))
    table.sort(key=lambda row: row["center_cm1"])
    rss = float(np.dot(residuals, residuals))
    return {
        "status": status,
        "shape": shape,
        "mean_parameter_count": len(initial),
        "diagnostics": diagnostics,
        "peaks": table,
        "rss": rss,
        "rmse": float(np.sqrt(rss / len(x))),
        "residual_lag1_correlation": lag1,
        "residual_correlation_limit": residual_correlation_limit,
        "reduced_chi_square": reduced_chi_square,
        "degrees_of_freedom": dof,
        "jacobian_rank": rank,
        "input_is_affine_within_tolerance": no_peak_contrast,
        "affine_tolerance_relative": 1e-12,
        "jacobian_condition": condition if np.isfinite(condition) else None,
        "covariance": covariance.tolist() if covariance is not None else None,
        "covariance_parameter_order": "input peak order: "
        + ", ".join(PARAMETERS[shape])
        + "; background intercept, slope",
        "covariance_mode": "absolute measurement sigma"
        if errors is not None
        else "residual-scaled iid noise",
        "uncertainty_scope": "local linear approximation conditional on fixed baseline and model; excludes baseline and calibration uncertainty",
        "optimizer_message": str(result.message),
        "evaluations": int(result.nfev),
        "fitted": fitted,
        "residuals": residuals,
        "components": components(result.x),
        "residual_background": result.x[-2] + result.x[-1] * x0,
    }
