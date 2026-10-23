"""Penalized and explicit reference-region baselines with solver diagnostics."""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.special import expit


def estimate_baseline(x, y, config):
    """Compatibility interface returning just the estimated curve."""
    return baseline_with_diagnostics(x, y, config)[0]


def baseline_with_diagnostics(x, y, config):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if (
        x.ndim != 1
        or x.shape != y.shape
        or len(x) < 3
        or not np.isfinite([x, y]).all()
        or np.any(np.diff(x) <= 0)
    ):
        raise ValueError("Baseline arrays need matching finite increasing samples")
    config = dict(config)
    method = config.pop("method", "none")
    details = {"method": method, "resampled": False, "input_points": len(x), "diagnostics": []}
    if method == "none":
        if config:
            raise ValueError("none baseline takes no parameters")
        return np.zeros_like(y), details
    if method in {"asls", "arpls"}:
        allowed = {"lam", "iterations", "resample", "p" if method == "asls" else "tolerance"}
        if set(config) - allowed:
            raise ValueError(f"Unknown {method} parameter")
        resample = config.pop("resample", False)
        if type(resample) is not bool:
            raise ValueError("Penalized baseline resample must be true or false")
        spacing = np.diff(x)
        if not np.allclose(spacing, np.median(spacing), rtol=1e-3, atol=1e-9):
            if resample:
                # Only the baseline estimator uses this explicit interpolation;
                # original spectral samples still enter every fit and export.
                grid = np.linspace(x[0], x[-1], len(x))
                baseline, details = baseline_with_diagnostics(
                    grid, np.interp(grid, x, y), {"method": method, **config}
                )
                details.update(
                    resampled=True,
                    interpolation="linear; uniform grid for baseline only; fit samples unchanged",
                )
                return np.interp(x, grid, baseline), details
            raise ValueError(
                "Penalized baseline uses an index-space penalty and requires uniform sampling"
            )
        lam, p, iterations = (
            config.get("lam", 1e6),
            config.get("p", 0.01),
            config.get("iterations", 30 if method == "asls" else 100),
        )
        if not np.isfinite(lam) or lam <= 0 or not np.isfinite(p) or not 0 < p < 0.5:
            raise ValueError("Positive finite lam and 0 < p < 0.5 are required")
        if type(iterations) is not int or not 1 <= iterations <= 500:
            raise ValueError("iterations must be an integer in [1, 500]")
        n = len(y)
        difference = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n), format="csc")
        penalty = lam * (difference.T @ difference)
        weights = np.ones(n)
        if method == "arpls":
            tolerance = config.get("tolerance", 1e-3)
            if (
                type(tolerance) not in (int, float)
                or not np.isfinite(tolerance)
                or not 0 < tolerance < 1
            ):
                raise ValueError("arPLS tolerance must be in (0, 1)")
            # The second difference annihilates affine trends. Removing one
            # before solving improves conditioning without changing the model.
            trend = np.polynomial.Polynomial.fit(x, y, 1)(x)
            contrast = y - trend
            scale = max(float(np.max(np.abs(y))), np.finfo(float).tiny)
            if np.max(np.abs(contrast)) <= scale * 1e-12:
                details.update(
                    converged=True, iterations=0, termination="affine input", tolerance=tolerance
                )
                return trend, details
            converged, change, reason = False, None, "iteration limit reached"
            for iteration in range(1, iterations + 1):
                baseline = spsolve(
                    sparse.diags(weights, format="csc") + penalty, weights * contrast
                )
                residual = contrast - baseline
                negative = residual[residual < 0]
                spread = float(np.std(negative)) if len(negative) >= 2 else 0
                if spread <= np.finfo(float).eps * scale * 32:
                    reason = "negative residual spread is degenerate"
                    break
                # Baek et al., DOI 10.1039/C4AN01061B. Population SD (ddof=0).
                updated = expit(-2 * (residual - (2 * spread - float(np.mean(negative)))) / spread)
                # A machine-precision floor only prevents singular arithmetic.
                updated = np.maximum(updated, np.finfo(float).eps)
                change = float(np.linalg.norm(updated - weights) / np.linalg.norm(weights))
                if change < tolerance:
                    converged, reason = True, "relative weight change below tolerance"
                    break
                weights = updated
            baseline = baseline + trend
            details.update(
                converged=converged,
                iterations=iteration,
                max_iterations=iterations,
                relative_weight_change=change,
                tolerance=tolerance,
                termination=reason,
                negative_residual_sd_ddof=0,
            )
            if not converged:
                details["diagnostics"].append(
                    "arPLS baseline iteration did not converge: " + reason
                )
        else:
            for _ in range(iterations):
                baseline = spsolve(sparse.diags(weights, format="csc") + penalty, weights * y)
                # Equality keeps a positive weight, including on exactly flat spectra.
                weights = np.where(y > baseline, p, 1 - p)
            details.update(
                iterations=iterations,
                converged=None,
                termination="fixed iteration budget; convergence not asserted",
            )
        if not np.all(np.isfinite(baseline)):
            raise ValueError("Penalized baseline did not produce a finite curve")
        return baseline, details
    if method == "polynomial":
        if set(config) - {"degree", "regions"}:
            raise ValueError("Unknown polynomial parameter")
        degree, regions = config.get("degree", 1), config.get("regions", [])
        if type(degree) is not int or not 0 <= degree <= 5:
            raise ValueError("Polynomial degree must be in [0, 5]")
        mask = np.zeros(len(x), dtype=bool)
        for low, high in regions:
            if not np.isfinite([low, high]).all() or low >= high:
                raise ValueError("Baseline regions must be finite increasing intervals")
            mask |= (x >= low) & (x <= high)
        if mask.sum() <= degree + 1:
            raise ValueError("Provide enough explicit baseline reference points")
        polynomial = np.polynomial.Polynomial.fit(x[mask], y[mask], degree)
        return polynomial(x), details
    raise ValueError("baseline.method must be none, asls, arpls or polynomial")
