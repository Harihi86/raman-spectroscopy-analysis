"""Height-normalized symmetric peak families and their analytic derivatives.

Voigt is the exact normal/Cauchy convolution, evaluated with the Faddeeva
function. sigma is a Gaussian SD; gamma is a Lorentzian HWHM, never an SD.
"""

import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import voigt_profile, wofz

FWHM_FACTOR = 2 * math.sqrt(2 * math.log(2))
PARAMETERS = {
    "gaussian": ("amplitude", "center", "sigma"),
    "lorentzian": ("amplitude", "center", "gamma"),
    "voigt": ("amplitude", "center", "sigma", "gamma"),
}


def gaussian(x, amplitude, center, sigma):
    return amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def lorentzian(x, amplitude, center, gamma):
    return amplitude / (1 + ((x - center) / gamma) ** 2)


def voigt(x, amplitude, center, sigma, gamma):
    return amplitude * voigt_profile(x - center, sigma, gamma) / voigt_profile(0, sigma, gamma)


def profile(x, params, shape):
    return {"gaussian": gaussian, "lorentzian": lorentzian, "voigt": voigt}[shape](x, *params)


def _voigt_partials(offset, sigma, gamma):
    """Unit-area derivatives, avoiding cancellation near the Cauchy limit.

    Differentiate the large-z expansion directly instead of subtracting nearly
    equal terms in -2*z*w(z) + 2i/sqrt(pi). Values still use SciPy's ufunc.
    """
    offset = np.atleast_1d(np.asarray(offset, dtype=float))
    q = offset + 1j * gamma
    large = np.abs(q) > 40 * sigma
    dx, ds, dg = [np.empty(len(offset)) for _ in range(3)]
    if np.any(~large):
        z = q[~large] / (math.sqrt(2) * sigma)
        w = wofz(z)
        dw = -2 * z * w + 2j / math.sqrt(math.pi)
        dx[~large] = dw.real / (2 * sigma**2 * math.sqrt(math.pi))
        ds[~large] = -(w.real + (z * dw).real) / (sigma**2 * math.sqrt(2 * math.pi))
        dg[~large] = -dw.imag / (2 * sigma**2 * math.sqrt(math.pi))
    if np.any(large):
        inverse = 1 / q[large]
        ratio = (sigma * inverse) ** 2
        dx_series = inverse**2 * np.polynomial.polynomial.polyval(
            ratio, [1, 3, 15, 105, 945, 10395, 135135]
        )
        ds_series = (
            sigma
            * inverse**3
            * np.polynomial.polynomial.polyval(ratio, [2, 12, 90, 840, 9450, 124740])
        )
        dx[large] = (-1j * dx_series).real / math.pi
        ds[large] = (1j * ds_series).real / math.pi
        dg[large] = dx_series.real / math.pi
    return dx, ds, dg


def profile_jacobian(x, params, shape):
    amplitude, center, width = params[:3]
    offset = x - center
    if shape == "gaussian":
        z = offset / width
        value = np.exp(-0.5 * z**2)
        return np.column_stack(
            [value, amplitude * value * z / width, amplitude * value * z**2 / width]
        )
    if shape == "lorentzian":
        z = offset / width
        value = 1 / (1 + z**2)
        return np.column_stack(
            [value, 2 * amplitude * z * value**2 / width, 2 * amplitude * z**2 * value**2 / width]
        )
    sigma, gamma = params[2:4]
    v0 = voigt_profile(0, sigma, gamma)
    h = voigt_profile(offset, sigma, gamma) / v0
    dx, ds_raw, dg_raw = _voigt_partials(offset, sigma, gamma)
    _, ds0, dg0 = _voigt_partials([0], sigma, gamma)
    dc = -dx / v0
    ds = (ds_raw - h * ds0[0]) / v0
    dg = (dg_raw - h * dg0[0]) / v0
    return np.column_stack([h, amplitude * dc, amplitude * ds, amplitude * dg])


def peak_fwhm(params, shape):
    if shape == "gaussian":
        return FWHM_FACTOR * params[2]
    if shape == "lorentzian":
        return 2 * params[2]
    sigma, gamma = params[2:4]
    # Solve the exact half-height equation, not an empirical FWHM approximation.
    half_height = voigt_profile(0, sigma, gamma) / 2
    upper = FWHM_FACTOR * sigma + 2 * gamma
    return 2 * brentq(
        lambda distance: voigt_profile(distance, sigma, gamma) - half_height,
        0,
        upper,
        xtol=max(upper * 1e-13, np.finfo(float).tiny),
    )


def peak_summary(params, shape, covariance=None):
    """Full-line area and local covariance propagation, including cross terms."""
    amplitude, center, width = [float(v) for v in params[:3]]
    fwhm = float(peak_fwhm(params, shape))
    if shape == "gaussian":
        area = amplitude * width * math.sqrt(2 * math.pi)
        area_gradient = np.array([width, 0, amplitude]) * math.sqrt(2 * math.pi)
        width_gradient = np.array([0, 0, FWHM_FACTOR])
    elif shape == "lorentzian":
        area = math.pi * amplitude * width
        area_gradient = math.pi * np.array([width, 0, amplitude])
        width_gradient = np.array([0, 0, 2])
    else:
        sigma, gamma = params[2:4]
        v0 = voigt_profile(0, sigma, gamma)
        area = amplitude / v0
        _, ds0, dg0 = _voigt_partials([0], sigma, gamma)
        dv_sigma, dv_gamma = ds0[0], dg0[0]
        area_gradient = np.array(
            [1 / v0, 0, -amplitude * dv_sigma / v0**2, -amplitude * dv_gamma / v0**2]
        )
        # Implicit differentiation of H(FWHM/2, sigma, gamma) = 1/2.
        edge_jac = profile_jacobian(np.array([center + fwhm / 2]), params, shape)[0]
        width_gradient = np.array(
            [0, 0, 2 * edge_jac[2] / edge_jac[1], 2 * edge_jac[3] / edge_jac[1]]
        )
    row = {
        "shape": shape,
        "amplitude": amplitude,
        "center_cm1": center,
        "sigma_cm1": width if shape != "lorentzian" else None,
        "gamma_cm1": float(params[-1]) if shape != "gaussian" else None,
        "fwhm_cm1": fwhm,
        "area": float(area),
        "amplitude_se": None,
        "center_se_cm1": None,
        "sigma_se_cm1": None,
        "gamma_se_cm1": None,
        "fwhm_se_cm1": None,
        "area_se": None,
    }
    if covariance is not None:
        se = np.sqrt(np.maximum(np.diag(covariance), 0))
        row.update(
            amplitude_se=float(se[0]),
            center_se_cm1=float(se[1]),
            sigma_se_cm1=float(se[2]) if shape != "lorentzian" else None,
            gamma_se_cm1=float(se[-1]) if shape != "gaussian" else None,
            fwhm_se_cm1=float(np.sqrt(max(float(width_gradient @ covariance @ width_gradient), 0))),
            area_se=float(np.sqrt(max(float(area_gradient @ covariance @ area_gradient), 0))),
        )
    return row
