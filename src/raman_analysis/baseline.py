"""Asymmetric least squares and explicit reference-region polynomial baselines."""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def estimate_baseline(x, y, config):
    config = dict(config)
    method = config.pop("method", "none")
    if method == "none":
        if config:
            raise ValueError("none baseline takes no parameters")
        return np.zeros_like(y)
    if method == "asls":
        if set(config) - {"lam", "p", "iterations"}:
            raise ValueError("Unknown AsLS parameter")
        spacing = np.diff(x)
        if not np.allclose(spacing, np.median(spacing), rtol=1e-3, atol=1e-9):
            raise ValueError("AsLS uses an index-space penalty and requires uniform sampling")
        lam, p, iterations = (
            config.get("lam", 1e6),
            config.get("p", 0.01),
            config.get("iterations", 30),
        )
        if not np.isfinite(lam) or lam <= 0 or not np.isfinite(p) or not 0 < p < 0.5:
            raise ValueError("Positive finite lam and 0 < p < 0.5 are required")
        if type(iterations) is not int or not 1 <= iterations <= 500:
            raise ValueError("iterations must be an integer in [1, 500]")
        n = len(y)
        difference = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n), format="csc")
        penalty = lam * (difference.T @ difference)
        weights = np.ones(n)
        for _ in range(iterations):
            baseline = spsolve(sparse.diags(weights, format="csc") + penalty, weights * y)
            # Equality keeps a positive weight, including on exactly flat spectra.
            weights = np.where(y > baseline, p, 1 - p)
        if not np.all(np.isfinite(baseline)):
            raise ValueError("AsLS did not produce a finite baseline")
        return baseline
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
        return polynomial(x)
    raise ValueError("baseline.method must be none, asls or polynomial")
