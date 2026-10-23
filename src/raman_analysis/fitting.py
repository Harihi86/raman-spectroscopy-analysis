"""Bounded N-Gaussian fits with residual and local covariance diagnostics."""

import math

import numpy as np
from scipy.optimize import least_squares

FWHM_FACTOR = 2 * math.sqrt(2 * math.log(2))


def gaussian(x, amplitude, center, sigma):
    return amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def fit_peaks(x, y, peaks, *, errors=None, max_nfev=10000, residual_correlation_limit=0.35):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or x.shape != y.shape or len(x) < 8 or not np.isfinite([x, y]).all():
        raise ValueError("Fit arrays must have matching shapes and finite values")
    if not np.all(np.diff(x) > 0):
        raise ValueError("Fit axis must be strictly increasing")
    if not 1 <= len(peaks) <= 8:
        raise ValueError("Specify 1..8 peaks explicitly")
    if type(max_nfev) is not int or not 1 <= max_nfev <= 100000:
        raise ValueError("max_nfev must be in [1, 100000]")
    if residual_correlation_limit is not None and not 0 < residual_correlation_limit < 1:
        raise ValueError("residual_correlation_limit must be in (0, 1) or null")
    if len(x) <= 3 * len(peaks) + 2:
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
        if set(peak) - {"center", "sigma", "amplitude", "center_bounds", "sigma_bounds"}:
            raise ValueError("Unknown peak configuration field")
        center = float(peak["center"])
        sigma = float(peak.get("sigma", span / (12 * len(peaks))))
        amplitude = float(peak.get("amplitude", max(float(np.ptp(y)), 1e-3)))
        cb = peak.get("center_bounds", [float(x[0]), float(x[-1])])
        sb = peak.get("sigma_bounds", [step / 4, span])
        if (
            len(cb) != 2
            or len(sb) != 2
            or not np.isfinite([center, sigma, amplitude, *cb, *sb]).all()
        ):
            raise ValueError("Peak parameters and bounds must be finite")
        if (
            not x[0] <= cb[0] < center < cb[1] <= x[-1]
            or not 0 < sb[0] < sigma < sb[1]
            or amplitude <= 0
        ):
            raise ValueError("Initial peaks must lie strictly inside positive physical bounds")
        seeds.append(center)
        initial.extend([amplitude, center, sigma])
        lower.extend([0, cb[0], sb[0]])
        upper.extend([np.inf, cb[1], sb[1]])
    if len(set(seeds)) != len(seeds):
        raise ValueError("Initial peak centers must be distinct")
    initial.extend([float(np.percentile(y, 10)), 0.0])
    lower.extend([-np.inf, -np.inf])
    upper.extend([np.inf, np.inf])

    def components(params):
        return [gaussian(x, *params[3 * i : 3 * i + 3]) for i in range(len(peaks))]

    def model(params):
        return np.sum(components(params), axis=0) + params[-2] + params[-1] * x0

    def residual(params):
        r = model(params) - y
        return r if errors is None else r / errors

    result = least_squares(
        residual, initial, bounds=(lower, upper), max_nfev=max_nfev, x_scale="jac", loss="linear"
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
        diagnostics.append("flat or affine input contains no identifiable Gaussian peak")
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
        amp, center, sigma = [float(v) for v in result.x[3 * i : 3 * i + 3]]
        area = amp * sigma * math.sqrt(2 * math.pi)
        row = {
            "amplitude": amp,
            "center_cm1": center,
            "sigma_cm1": sigma,
            "fwhm_cm1": FWHM_FACTOR * sigma,
            "area": area,
            "amplitude_se": None,
            "center_se_cm1": None,
            "sigma_se_cm1": None,
            "fwhm_se_cm1": None,
            "area_se": None,
        }
        if covariance is not None:
            block = covariance[3 * i : 3 * i + 3, 3 * i : 3 * i + 3]
            se = np.sqrt(np.maximum(np.diag(block), 0))
            gradient = np.array([sigma, 0, amp]) * math.sqrt(2 * math.pi)
            row.update(
                amplitude_se=float(se[0]),
                center_se_cm1=float(se[1]),
                sigma_se_cm1=float(se[2]),
                fwhm_se_cm1=float(FWHM_FACTOR * se[2]),
                area_se=float(math.sqrt(max(float(gradient @ block @ gradient), 0))),
            )
        table.append(row)
    table.sort(key=lambda row: row["center_cm1"])
    rss = float(np.dot(residuals, residuals))
    return {
        "status": status,
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
        "covariance_parameter_order": "input peak order: amplitude, center, sigma; background intercept, slope",
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
