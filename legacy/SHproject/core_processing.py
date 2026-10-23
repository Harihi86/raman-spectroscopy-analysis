#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import copy
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, least_squares
from scipy.signal import find_peaks, savgol_filter
from scipy.sparse.linalg import spsolve


@dataclass
class Spectrum:
    label: str
    pressure: float
    x: np.ndarray
    y_raw: np.ndarray
    y_rate_norm: np.ndarray
    baseline: np.ndarray
    y_corrected: np.ndarray
    zero_reference_offset: float = 0.0


def deep_copy_config(cfg: Dict[str, Any]) -> Dict[str, Any]:
    return copy.deepcopy(cfg)


def deep_update(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    result = deep_copy_config(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(result.get(k), dict):
            result[k] = deep_update(result[k], v)
        else:
            result[k] = copy.deepcopy(v)
    return result


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)

try:
    from tqdm.auto import tqdm as _tqdm  # type: ignore
except Exception:
    _tqdm = None


class _SimpleProgress:
    def __init__(self, total: Optional[int] = None, desc: str = '', enabled: bool = True, leave: bool = True) -> None:
        self.total = total if total is not None else 0
        self.desc = desc
        self.enabled = enabled
        self.leave = leave
        self.n = 0
        self.postfix = ''
        self._last_render = 0.0
        if self.enabled:
            self._render(force=True)

    def set_postfix_str(self, text: str) -> None:
        self.postfix = text
        if self.enabled:
            self._render(force=True)

    def update(self, n: int = 1) -> None:
        self.n += n
        if self.enabled:
            self._render(force=False)

    def close(self) -> None:
        if not self.enabled:
            return
        self._render(force=True)
        if self.leave:
            sys.stdout.write('\n')
        else:
            sys.stdout.write('\r' + ' ' * 140 + '\r')
        sys.stdout.flush()

    def _render(self, force: bool = False) -> None:
        now = time.time()
        if (not force) and (now - self._last_render < 0.08):
            return
        self._last_render = now
        total = self.total if self.total and self.total > 0 else None
        if total is None:
            msg = f"{self.desc}: {self.n}"
        else:
            frac = min(max(self.n / total, 0.0), 1.0)
            width = 24
            filled = int(round(width * frac))
            bar = '#' * filled + '-' * (width - filled)
            pct = int(round(frac * 100.0))
            msg = f"{self.desc}: [{bar}] {self.n}/{total} {pct:3d}%"
        if self.postfix:
            msg += f" | {self.postfix}"
        sys.stdout.write('\r' + msg[:180])
        sys.stdout.flush()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()


def make_progress_bar(total: Optional[int], desc: str, enabled: bool = True, leave: bool = True):
    if not enabled:
        return _SimpleProgress(total=total, desc=desc, enabled=False, leave=leave)
    if _tqdm is not None:
        return _tqdm(total=total, desc=desc, leave=leave, dynamic_ncols=True)
    return _SimpleProgress(total=total, desc=desc, enabled=True, leave=leave)


def progress_iter(iterable: Iterable[Any], desc: str, total: Optional[int] = None, enabled: bool = True, leave: bool = False):
    if not enabled:
        for item in iterable:
            yield item
        return
    if _tqdm is not None:
        for item in _tqdm(iterable, total=total, desc=desc, leave=leave, dynamic_ncols=True):
            yield item
        return
    pbar = _SimpleProgress(total=total, desc=desc, enabled=True, leave=leave)
    try:
        for item in iterable:
            yield item
            pbar.update(1)
    finally:
        pbar.close()


def resolve_output_dir_from_cfg(cfg: Dict[str, Any]) -> Path:
    raw = str(cfg.get('output_dir', cfg.get('name', 'run_output'))).strip()
    path = Path(raw)
    if path.is_absolute():
        return path
    parts = path.parts
    if parts and parts[0].lower() == 'outputs':
        return path
    group = str(cfg.get('output_group', '')).strip().replace('\\', '/')
    if group:
        return Path('outputs') / Path(group) / path
    return Path('outputs') / path


def trim_xy(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = min(len(x), len(y))
    x = np.asarray(x[:n], dtype=float)
    y = np.asarray(y[:n], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) == 0:
        return x, y
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    ux, inv = np.unique(x, return_inverse=True)
    if len(ux) != len(x):
        y_avg = np.zeros_like(ux, dtype=float)
        counts = np.zeros_like(ux, dtype=float)
        for i, j in enumerate(inv):
            y_avg[j] += y[i]
            counts[j] += 1.0
        y = y_avg / np.maximum(counts, 1.0)
        x = ux
    return x, y


def in_window(x: np.ndarray, window: Optional[Tuple[float, float]]) -> np.ndarray:
    if window is None:
        return np.ones_like(x, dtype=bool)
    lo, hi = window
    return (x >= lo) & (x <= hi)


def build_exclusion_mask(x: np.ndarray, exclude_ranges: Optional[Sequence[Tuple[float, float]]]) -> np.ndarray:
    mask = np.ones_like(x, dtype=bool)
    if not exclude_ranges:
        return mask
    for lo, hi in exclude_ranges:
        mask &= ((x < lo) | (x > hi))
    return mask


def apply_nan_masks(x: np.ndarray, y: np.ndarray, mask_ranges: Optional[Sequence[Tuple[float, float]]]) -> np.ndarray:
    out = np.asarray(y, dtype=float).copy()
    if not mask_ranges:
        return out
    for lo, hi in _normalize_mask_ranges(mask_ranges):
        bad = (x >= lo) & (x <= hi)
        out[bad] = np.nan
    return out


def _normalize_mask_ranges(mask_ranges: Optional[Sequence[Tuple[float, float]]]) -> List[Tuple[float, float]]:
    if not mask_ranges:
        return []
    out: List[Tuple[float, float]] = []
    for item in mask_ranges:
        if item is None:
            continue
        if len(item) < 2:
            continue
        lo = float(item[0])
        hi = float(item[1])
        if not (np.isfinite(lo) and np.isfinite(hi)):
            continue
        if hi < lo:
            lo, hi = hi, lo
        out.append((lo, hi))
    if not out:
        return []
    out.sort(key=lambda t: (t[0], t[1]))
    merged: List[Tuple[float, float]] = [out[0]]
    for lo, hi in out[1:]:
        plo, phi = merged[-1]
        if lo <= phi:
            merged[-1] = (plo, max(phi, hi))
        else:
            merged.append((lo, hi))
    return merged


def _mask_rule_matches(rule: Dict[str, Any], spectrum_index: Optional[int], pressure: Optional[float], label: Optional[str]) -> bool:
    idx1 = None if spectrum_index is None else int(spectrum_index) + 1
    index_min = rule.get('index_min')
    index_max = rule.get('index_max')
    if index_min is not None and (idx1 is None or idx1 < int(index_min)):
        return False
    if index_max is not None and (idx1 is None or idx1 > int(index_max)):
        return False

    pressure_min = rule.get('pressure_min')
    pressure_max = rule.get('pressure_max')
    if pressure_min is not None and (pressure is None or float(pressure) < float(pressure_min)):
        return False
    if pressure_max is not None and (pressure is None or float(pressure) > float(pressure_max)):
        return False

    labels = rule.get('labels')
    if labels is not None:
        labels_set = {str(v) for v in labels}
        if label is None or str(label) not in labels_set:
            return False

    label_contains = rule.get('label_contains')
    if label_contains is not None:
        needle = str(label_contains)
        if label is None or needle not in str(label):
            return False

    return True


def get_target_mask_ranges(
    cfg: Dict[str, Any],
    target: str,
    spectrum_index: Optional[int] = None,
    pressure: Optional[float] = None,
    label: Optional[str] = None,
) -> List[Tuple[float, float]]:
    target = str(target).lower()
    ranges: List[Tuple[float, float]] = []
    ranges.extend(cfg.get('mask_ranges', []))
    ranges.extend(cfg.get(f'{target}_mask_ranges', []))

    for rule in cfg.get('mask_rules', []):
        rule_targets = {str(t).lower() for t in rule.get('targets', ['all'])}
        if 'all' not in rule_targets and target not in rule_targets:
            continue
        if not _mask_rule_matches(rule, spectrum_index, pressure, label):
            continue
        ranges.extend(rule.get('ranges', []))

    return _normalize_mask_ranges(ranges)


def fill_masked_regions_linear(x: np.ndarray, y: np.ndarray, mask_ranges: Optional[Sequence[Tuple[float, float]]]) -> np.ndarray:
    out = np.asarray(y, dtype=float).copy()
    ranges = _normalize_mask_ranges(mask_ranges)
    if not ranges:
        return out

    x = np.asarray(x, dtype=float)
    finite = np.isfinite(x) & np.isfinite(out)
    if finite.sum() < 2:
        return out

    bad = np.zeros_like(out, dtype=bool)
    for lo, hi in ranges:
        bad |= ((x >= lo) & (x <= hi))
    if not np.any(bad):
        return out

    good = finite & (~bad)
    if good.sum() < 2:
        return out

    out[bad] = np.interp(x[bad], x[good], out[good])
    return out


def _split_contiguous_segments(xf: np.ndarray, yf: np.ndarray) -> List[Tuple[np.ndarray, np.ndarray]]:
    if len(xf) == 0:
        return []
    if len(xf) == 1:
        return [(xf, yf)]
    diffs = np.diff(xf)
    pos = diffs[diffs > 0]
    if pos.size == 0:
        return [(xf, yf)]
    dx_ref = float(np.nanmedian(pos))
    if not np.isfinite(dx_ref) or dx_ref <= 0:
        return [(xf, yf)]
    gap_threshold = max(4.0, 3.0 * dx_ref)
    breaks = np.where(diffs > gap_threshold)[0]
    if breaks.size == 0:
        return [(xf, yf)]
    segments: List[Tuple[np.ndarray, np.ndarray]] = []
    start = 0
    for br in breaks:
        stop = int(br) + 1
        segments.append((xf[start:stop], yf[start:stop]))
        start = stop
    segments.append((xf[start:], yf[start:]))
    return segments


def _segment_prominence_score(xseg: np.ndarray, yseg: np.ndarray) -> float:
    if len(xseg) < 3:
        return -np.inf
    m0, b0 = _linear_guess(xseg, yseg)
    resid = yseg - (m0 * xseg + b0)
    if not np.isfinite(resid).any():
        return -np.inf
    return float(np.nanmax(resid))


def _select_primary_fit_segment(xf: np.ndarray, yf: np.ndarray, min_points: int) -> Tuple[np.ndarray, np.ndarray]:
    segments = _split_contiguous_segments(xf, yf)
    if len(segments) <= 1:
        return xf, yf

    min_keep = max(8, min_points // 2)
    candidates = [seg for seg in segments if len(seg[0]) >= min_keep]
    if not candidates:
        candidates = segments

    scored = []
    for xseg, yseg in candidates:
        scored.append((_segment_prominence_score(xseg, yseg), len(xseg), xseg, yseg))
    scored.sort(key=lambda t: (t[0], t[1]))
    _, _, xbest, ybest = scored[-1]
    return xbest, ybest


def _fit_center_bounds(xf: np.ndarray, fit_window: Tuple[float, float]) -> Tuple[float, float]:
    if len(xf) == 0:
        return float(fit_window[0]), float(fit_window[1])
    xlo = float(np.nanmin(xf))
    xhi = float(np.nanmax(xf))
    span = xhi - xlo
    if not np.isfinite(span) or span <= 0:
        return xlo, xhi
    edge_margin = min(12.0, max(2.0, 0.08 * span))
    clo = max(float(fit_window[0]), xlo + edge_margin)
    chi = min(float(fit_window[1]), xhi - edge_margin)
    if chi <= clo:
        clo, chi = xlo, xhi
    return float(clo), float(chi)


def _safe_savgol(y: np.ndarray, window_length: int, polyorder: int) -> np.ndarray:
    arr = np.asarray(y, dtype=float)
    n = len(arr)
    if n < 5:
        return arr.copy()
    win = int(max(window_length, 5))
    if win % 2 == 0:
        win += 1
    if win > n:
        win = n if (n % 2 == 1) else (n - 1)
    if win < 5:
        return arr.copy()
    po = int(max(1, polyorder))
    if po >= win:
        po = win - 2
    if po < 1:
        return arr.copy()
    try:
        return np.asarray(savgol_filter(arr, window_length=win, polyorder=po, mode='interp'), dtype=float)
    except Exception:
        return arr.copy()


def _prepare_peak_search_residual(
    xf: np.ndarray,
    yf: np.ndarray,
    search_window: Optional[Tuple[float, float]],
    local_bg_fraction: float,
    smooth_window: int,
    smooth_polyorder: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    x_use = np.asarray(xf, dtype=float)
    y_use = np.asarray(yf, dtype=float)
    if search_window is not None:
        sw_mask = in_window(x_use, search_window)
        if np.count_nonzero(sw_mask) >= 5:
            x_use = x_use[sw_mask]
            y_use = y_use[sw_mask]
    bg, m0, b0 = _edge_linear_background(x_use, y_use, fraction=local_bg_fraction)
    resid = np.asarray(y_use - bg, dtype=float)
    resid_smooth = _safe_savgol(resid, window_length=smooth_window, polyorder=smooth_polyorder)
    return x_use, y_use, resid_smooth, float(m0), float(b0)


def detect_prominent_peak_candidate(
    xf: np.ndarray,
    yf: np.ndarray,
    search_window: Optional[Tuple[float, float]] = None,
    edge_margin: float = 8.0,
    min_prominence_ratio: float = 0.08,
    smooth_window: int = 11,
    smooth_polyorder: int = 2,
    local_bg_fraction: float = 0.12,
) -> Optional[Dict[str, float]]:
    if len(xf) < 5:
        return None
    x_use, y_use, resid_smooth, m0, b0 = _prepare_peak_search_residual(
        xf, yf, search_window=search_window, local_bg_fraction=local_bg_fraction,
        smooth_window=smooth_window, smooth_polyorder=smooth_polyorder,
    )
    if len(x_use) < 5 or not np.isfinite(resid_smooth).any():
        return None

    edge_margin = float(max(edge_margin, 0.0))
    interior = np.ones_like(x_use, dtype=bool)
    if len(x_use) >= 3 and edge_margin > 0:
        interior &= (x_use >= float(x_use[0]) + edge_margin) & (x_use <= float(x_use[-1]) - edge_margin)
        if np.count_nonzero(interior) < 5:
            interior = np.ones_like(x_use, dtype=bool)

    x_int = x_use[interior]
    y_int = resid_smooth[interior]
    if len(x_int) < 5 or not np.isfinite(y_int).any():
        return None

    scale = float(np.nanpercentile(y_int, 95.0) - np.nanpercentile(y_int, 20.0))
    if not np.isfinite(scale) or scale <= 0:
        scale = float(np.nanmax(y_int) - np.nanmin(y_int)) if np.isfinite(np.nanmax(y_int) - np.nanmin(y_int)) else 0.0
    prominence_abs = max(1e-9, float(min_prominence_ratio) * max(scale, 1e-9))

    peaks, props = find_peaks(y_int, prominence=prominence_abs)
    if len(peaks) == 0:
        return None

    prominences = np.asarray(props.get('prominences', np.zeros(len(peaks))), dtype=float)
    heights = y_int[peaks]
    scores = prominences + 0.15 * np.maximum(heights, 0.0)
    best = int(np.nanargmax(scores))
    center = float(x_int[peaks[best]])
    height = float(heights[best])
    prominence = float(prominences[best]) if len(prominences) > best else np.nan
    return {
        'center': center,
        'height': height,
        'prominence': prominence,
        'bg_slope': float(m0),
        'bg_intercept': float(b0),
        'search_x_min': float(np.nanmin(x_use)),
        'search_x_max': float(np.nanmax(x_use)),
    }


def detect_read_mode(df: pd.DataFrame) -> str:
    ncols = df.shape[1]
    numeric_counts = [pd.to_numeric(df.iloc[:, i], errors="coerce").notna().sum() for i in range(ncols)]

    shared_like = ncols >= 2 and numeric_counts[0] > 10 and any(c > 10 for c in numeric_counts[1:])

    xy_like = ncols >= 2 and ncols % 2 == 0
    if xy_like:
        checks = []
        for i in range(0, ncols, 2):
            xvals = pd.to_numeric(df.iloc[:, i], errors="coerce").dropna().to_numpy()
            checks.append(len(xvals) > 10)
        xy_like = all(checks)

    if xy_like:
        return "xy_pairs"
    if shared_like:
        return "shared_x"
    raise ValueError("Could not auto-detect Excel layout. Set read_mode manually.")


def read_excel_spectra(
    excel_path: str,
    sheet_name: Any = 0,
    read_mode: str = "auto",
    labels: Optional[Sequence[str]] = None,
) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    df = pd.read_excel(excel_path, sheet_name=sheet_name, header=None)
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")

    if read_mode == "auto":
        read_mode = detect_read_mode(df)

    spectra: List[Tuple[str, np.ndarray, np.ndarray]] = []

    if read_mode == "shared_x":
        x = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
        for j in range(1, df.shape[1]):
            y = pd.to_numeric(df.iloc[:, j], errors="coerce").to_numpy(dtype=float)
            xv, yv = trim_xy(x, y)
            if len(xv) == 0:
                continue
            label = labels[j - 1] if labels is not None and j - 1 < len(labels) else f"S{len(spectra) + 1:02d}"
            spectra.append((label, xv, yv))
    elif read_mode == "xy_pairs":
        for j in range(0, df.shape[1], 2):
            x = pd.to_numeric(df.iloc[:, j], errors="coerce").to_numpy(dtype=float)
            y = pd.to_numeric(df.iloc[:, j + 1], errors="coerce").to_numpy(dtype=float)
            xv, yv = trim_xy(x, y)
            if len(xv) == 0:
                continue
            idx = j // 2
            label = labels[idx] if labels is not None and idx < len(labels) else f"S{len(spectra) + 1:02d}"
            spectra.append((label, xv, yv))
    else:
        raise ValueError("read_mode must be one of: auto, shared_x, xy_pairs")

    if not spectra:
        raise ValueError("No valid spectra were read from Excel.")
    return spectra


def build_rate_factors(rate_cfg: Dict[str, Any], nspectra: int) -> Optional[np.ndarray]:
    if not rate_cfg.get("enabled", False):
        return None

    direct = rate_cfg.get("direct_factors", [])
    if direct:
        if len(direct) != nspectra:
            raise ValueError("Length of direct_factors must match number of spectra.")
        factors = np.asarray(direct, dtype=float)
    else:
        exposure_s = rate_cfg.get("exposure_s", [])
        n_frames = rate_cfg.get("n_frames", [])
        if len(exposure_s) != nspectra or len(n_frames) != nspectra:
            raise ValueError("When rate normalization is enabled, exposure_s and n_frames must both match number of spectra.")
        factors = np.asarray(exposure_s, dtype=float) * np.asarray(n_frames, dtype=float)

    if np.any(factors <= 0):
        raise ValueError("All rate normalization factors must be > 0.")
    return factors


def als_baseline(y: np.ndarray, lam: float = 1e7, p: float = 0.01, niter: int = 30) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return np.zeros_like(y)
    D = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n), format="csc")
    w = np.ones(n, dtype=float)
    for _ in range(niter):
        W = sparse.spdiags(w, 0, n, n, format="csc")
        Z = W + lam * (D.T @ D)
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1.0 - p) * (y < z)
    return np.asarray(z, dtype=float)


def make_common_grid(window: Tuple[float, float], step: float) -> np.ndarray:
    lo, hi = window
    n = int(round((hi - lo) / step)) + 1
    return np.linspace(lo, hi, n)


def interpolate_spectrum(x: np.ndarray, y: np.ndarray, xnew: np.ndarray) -> np.ndarray:
    f = interp1d(x, y, kind="linear", bounds_error=False, fill_value=np.nan)
    return f(xnew)


def minmax_normalize_ignore_nan(y: np.ndarray) -> np.ndarray:
    out = np.asarray(y, dtype=float).copy()
    finite = np.isfinite(out)
    if finite.sum() == 0:
        return np.full_like(out, np.nan, dtype=float)
    ymin = np.nanmin(out)
    ymax = np.nanmax(out)
    if not np.isfinite(ymin) or not np.isfinite(ymax) or ymax <= ymin:
        out[finite] = 0.0
        return out
    out[finite] = (out[finite] - ymin) / (ymax - ymin)
    return out


def max_normalize_ignore_nan(y: np.ndarray) -> np.ndarray:
    out = np.asarray(y, dtype=float).copy()
    finite = np.isfinite(out)
    if finite.sum() == 0:
        return np.full_like(out, np.nan, dtype=float)
    ymax = np.nanmax(out)
    if not np.isfinite(ymax) or ymax == 0:
        out[finite] = 0.0
        return out
    out[finite] = out[finite] / ymax
    return out


def normalize_vector(y: np.ndarray, mode: str) -> np.ndarray:
    mode = (mode or "none").lower()
    if mode == "none":
        return np.asarray(y, dtype=float)
    if mode in {"minmax", "per_spectrum_minmax"}:
        return minmax_normalize_ignore_nan(y)
    if mode in {"max", "peakmax", "per_spectrum_max"}:
        return max_normalize_ignore_nan(y)
    raise ValueError(f"Unsupported normalization mode: {mode}")


def _pressure_rule_matches(pressure: float, rule: Dict[str, Any]) -> bool:
    pmin = rule.get("pressure_min", None)
    pmax = rule.get("pressure_max", None)
    if pmin is not None and pressure < float(pmin):
        return False
    if pmax is not None and pressure > float(pmax):
        return False
    return True


def _select_heatmap_norm_window(
    pressure: float,
    default_window: Tuple[float, float],
    norm_window_rules: Optional[Sequence[Dict[str, Any]]],
) -> Tuple[float, float]:
    if not norm_window_rules:
        return tuple(default_window)
    for rule in norm_window_rules:
        if _pressure_rule_matches(float(pressure), rule):
            win = rule.get("window", None)
            if win is not None:
                return (float(win[0]), float(win[1]))
    return tuple(default_window)


def compute_zero_reference_offset(
    x: np.ndarray,
    y: np.ndarray,
    window: Optional[Tuple[float, float]],
    statistic: str = "median",
    exclude_ranges: Optional[Sequence[Tuple[float, float]]] = None,
) -> float:
    mask = in_window(x, window) & build_exclusion_mask(x, exclude_ranges) & np.isfinite(y)
    if np.count_nonzero(mask) == 0:
        return 0.0
    vals = np.asarray(y[mask], dtype=float)
    stat = str(statistic).strip().lower()
    if stat == 'mean':
        out = float(np.nanmean(vals))
    elif stat in {'min', 'minimum'}:
        out = float(np.nanmin(vals))
    elif stat in {'p10', 'percentile10'}:
        out = float(np.nanpercentile(vals, 10.0))
    elif stat in {'p25', 'percentile25'}:
        out = float(np.nanpercentile(vals, 25.0))
    else:
        out = float(np.nanmedian(vals))
    return out if np.isfinite(out) else 0.0


def maybe_apply_zero_reference(
    x: np.ndarray,
    ycorr: np.ndarray,
    cfg: Optional[Dict[str, Any]],
) -> Tuple[np.ndarray, float]:
    if not cfg or not bool(cfg.get('enabled', False)):
        return np.asarray(ycorr, dtype=float), 0.0
    window_raw = cfg.get('window')
    window = tuple(window_raw) if isinstance(window_raw, (list, tuple)) and len(window_raw) == 2 else None
    exclude_ranges = _normalize_mask_ranges(cfg.get('exclude_ranges', []))
    offset = compute_zero_reference_offset(
        np.asarray(x, dtype=float),
        np.asarray(ycorr, dtype=float),
        window=window,
        statistic=str(cfg.get('statistic', 'median')),
        exclude_ranges=exclude_ranges,
    )
    return np.asarray(ycorr, dtype=float) - offset, float(offset)


def _exclude_heatmap_rules_from_norm(
    norm_mask: np.ndarray,
    xgrid: np.ndarray,
    pressure: float,
    fill_rules: Optional[Sequence[Dict[str, Any]]],
) -> np.ndarray:
    out = np.asarray(norm_mask, dtype=bool).copy()
    if not fill_rules:
        return out
    for rule in fill_rules:
        if not _pressure_rule_matches(float(pressure), rule):
            continue
        if not bool(rule.get("exclude_from_normalization", False)):
            continue
        x_min = float(rule.get("x_min", np.nanmin(xgrid)))
        x_max = float(rule.get("x_max", np.nanmax(xgrid)))
        out &= ~((xgrid >= x_min) & (xgrid <= x_max))
    return out


def _apply_heatmap_fill_rules(
    y_norm: np.ndarray,
    xgrid: np.ndarray,
    pressure: float,
    fill_rules: Optional[Sequence[Dict[str, Any]]],
) -> np.ndarray:
    out = np.asarray(y_norm, dtype=float).copy()
    if not fill_rules:
        return out
    for rule in fill_rules:
        if not _pressure_rule_matches(float(pressure), rule):
            continue
        x_min = float(rule.get("x_min", np.nanmin(xgrid)))
        x_max = float(rule.get("x_max", np.nanmax(xgrid)))
        fill_value = float(rule.get("fill_value", 0.0))
        mask = (xgrid >= x_min) & (xgrid <= x_max)
        out[mask] = fill_value
    return out


def gaussian(x: np.ndarray, amp: float, cen: float, sigma: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - cen) / sigma) ** 2)


def single_gaussian_with_linear_bg(x: np.ndarray, amp: float, cen: float, sigma: float, m: float, b: float) -> np.ndarray:
    return gaussian(x, amp, cen, sigma) + m * x + b


def double_gaussian_with_linear_bg(
    x: np.ndarray,
    amp1: float,
    cen1: float,
    sigma1: float,
    amp2: float,
    cen2: float,
    sigma2: float,
    m: float,
    b: float,
) -> np.ndarray:
    return gaussian(x, amp1, cen1, sigma1) + gaussian(x, amp2, cen2, sigma2) + m * x + b


def fit_stats(y: np.ndarray, yhat: np.ndarray, k: int) -> Dict[str, float]:
    resid = np.asarray(y, dtype=float) - np.asarray(yhat, dtype=float)
    n = len(resid)
    rss = float(np.nansum(resid ** 2))
    rss = max(rss, 1e-30)
    aic = float(n * np.log(rss / max(n, 1)) + 2 * k)
    bic = float(n * np.log(rss / max(n, 1)) + k * np.log(max(n, 1)))
    return {"rss": rss, "aic": aic, "bic": bic}


def prepare_fit_region(
    x: np.ndarray,
    y: np.ndarray,
    fit_window: Tuple[float, float],
    exclude_ranges: Optional[Sequence[Tuple[float, float]]],
) -> Tuple[np.ndarray, np.ndarray]:
    mask = in_window(x, fit_window) & build_exclusion_mask(x, exclude_ranges) & np.isfinite(y)
    xf = np.asarray(x[mask], dtype=float)
    yf = np.asarray(y[mask], dtype=float)
    return xf, yf


def _linear_guess(xf: np.ndarray, yf: np.ndarray) -> Tuple[float, float]:
    x1, x2 = float(xf[0]), float(xf[-1])
    y1, y2 = float(yf[0]), float(yf[-1])
    m0 = 0.0 if x2 == x1 else (y2 - y1) / (x2 - x1)
    b0 = y1 - m0 * x1
    return float(m0), float(b0)




def _edge_linear_background(xf: np.ndarray, yf: np.ndarray, fraction: float = 0.12) -> Tuple[np.ndarray, float, float]:
    if len(xf) < 4:
        m0, b0 = _linear_guess(xf, yf)
        return m0 * xf + b0, float(m0), float(b0)
    frac = float(np.clip(fraction, 0.05, 0.35))
    n_edge = max(6, int(math.ceil(len(xf) * frac)))
    n_edge = min(n_edge, max(2, len(xf) // 2))
    xs = np.concatenate([xf[:n_edge], xf[-n_edge:]])
    ys = np.concatenate([yf[:n_edge], yf[-n_edge:]])
    A = np.column_stack([xs, np.ones_like(xs)])
    try:
        m, b = np.linalg.lstsq(A, ys, rcond=None)[0]
    except Exception:
        m, b = _linear_guess(xf, yf)
    bg = m * xf + b
    return np.asarray(bg, dtype=float), float(m), float(b)


def _covariance_from_least_squares(jac: Optional[np.ndarray], residuals: np.ndarray, n_params: int) -> Optional[np.ndarray]:
    if jac is None:
        return None
    jac = np.asarray(jac, dtype=float)
    if jac.ndim != 2 or jac.shape[0] <= n_params:
        return None
    try:
        jtj = jac.T @ jac
        inv = np.linalg.pinv(jtj)
        dof = max(jac.shape[0] - n_params, 1)
        rss = float(np.nansum(np.asarray(residuals, dtype=float) ** 2))
        s_sq = rss / dof
        return inv * s_sq
    except Exception:
        return None


def _parse_double_pair(
    value: Optional[Sequence[float] | Sequence[Sequence[float]]],
    default: Tuple[float, float],
) -> Tuple[float, float]:
    if value is None:
        return float(default[0]), float(default[1])
    if len(value) < 2:
        return float(default[0]), float(default[1])
    try:
        return float(value[0]), float(value[1])
    except Exception:
        return float(default[0]), float(default[1])


def _parse_double_bounds(
    value: Optional[Sequence[Sequence[float]]],
    default1: Tuple[float, float],
    default2: Tuple[float, float],
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    if value is None or len(value) < 2:
        return (float(default1[0]), float(default1[1])), (float(default2[0]), float(default2[1]))
    try:
        b1 = (float(value[0][0]), float(value[0][1]))
        b2 = (float(value[1][0]), float(value[1][1]))
        return b1, b2
    except Exception:
        return (float(default1[0]), float(default1[1])), (float(default2[0]), float(default2[1]))



def fit_single_gaussian(
    x: np.ndarray,
    y: np.ndarray,
    fit_window: Tuple[float, float],
    exclude_ranges: Optional[Sequence[Tuple[float, float]]] = None,
    min_points: int = 20,
    search_window: Optional[Tuple[float, float]] = None,
    center_bounds: Optional[Sequence[float]] = None,
    sigma_bounds: Optional[Sequence[float]] = None,
    strict_local_peak: bool = False,
    edge_margin: float = 8.0,
    smooth_window: int = 11,
    smooth_polyorder: int = 2,
    min_prominence_ratio: float = 0.08,
    local_bg_fraction: float = 0.12,
    center_half_width: float = 18.0,
    peak_match_tolerance: float = 14.0,
    reject_boundary_fit: bool = False,
) -> Dict[str, Any]:
    xf_all, yf_all = prepare_fit_region(x, y, fit_window, exclude_ranges)
    xf, yf = _select_primary_fit_segment(xf_all, yf_all, min_points)
    result: Dict[str, Any] = {
        'success': False,
        'mode_used': 'single',
        'center': np.nan,
        'center_se': np.nan,
        'center1': np.nan,
        'center1_se': np.nan,
        'center2': np.nan,
        'center2_se': np.nan,
        'effective_center': np.nan,
        'effective_center_se': np.nan,
        'midpoint_center': np.nan,
        'fwhm': np.nan,
        'fwhm1': np.nan,
        'fwhm2': np.nan,
        'sigma': np.nan,
        'sigma1': np.nan,
        'sigma2': np.nan,
        'amplitude': np.nan,
        'amplitude1': np.nan,
        'amplitude2': np.nan,
        'area': np.nan,
        'area1': np.nan,
        'area2': np.nan,
        'bg_slope': np.nan,
        'bg_intercept': np.nan,
        'n_points': int(len(xf_all)),
        'fit_x': xf,
        'fit_y': yf,
        'fit_yhat': None,
        'fit_y_bg': None,
        'peak1_y': None,
        'peak2_y': None,
        'rss': np.nan,
        'aic': np.nan,
        'bic': np.nan,
        'message': '',
    }
    if len(xf) < min_points:
        result['message'] = 'too_few_points'
        return result

    cen_lo, cen_hi = _fit_center_bounds(xf, fit_window)
    xspan = max(float(xf[-1]) - float(xf[0]), 1.0)
    m0, b0 = _linear_guess(xf, yf)
    resid0 = yf - (m0 * xf + b0)
    if not np.isfinite(resid0).any():
        resid0 = yf - np.nanmedian(yf)

    candidate = None
    if strict_local_peak:
        candidate = detect_prominent_peak_candidate(
            xf,
            yf,
            search_window=search_window,
            edge_margin=float(edge_margin),
            min_prominence_ratio=float(min_prominence_ratio),
            smooth_window=int(smooth_window),
            smooth_polyorder=int(smooth_polyorder),
            local_bg_fraction=float(local_bg_fraction),
        )
        if candidate is None:
            result['message'] = 'no_prominent_local_peak'
            return result
        m0 = float(candidate['bg_slope'])
        b0 = float(candidate['bg_intercept'])
        resid0 = yf - (m0 * xf + b0)
        cen0 = float(candidate['center'])
        amp0 = float(max(candidate['height'], np.nanpercentile(resid0, 90), 1e-6))
        if center_bounds is None:
            half = float(max(center_half_width, 4.0))
            cen_lo = max(cen_lo, cen0 - half)
            cen_hi = min(cen_hi, cen0 + half)
            if cen_hi <= cen_lo:
                cen_lo, cen_hi = _fit_center_bounds(xf, fit_window)
    else:
        i_max = int(np.nanargmax(resid0))
        cen0 = float(np.clip(xf[i_max], cen_lo, cen_hi))
        amp0 = float(max(resid0[i_max], np.nanpercentile(resid0, 90), 1e-6))

    if center_bounds is not None and len(center_bounds) >= 2:
        cen_lo = max(cen_lo, float(center_bounds[0]))
        cen_hi = min(cen_hi, float(center_bounds[1]))
        if cen_hi <= cen_lo:
            result['message'] = 'invalid_center_bounds'
            return result
        cen0 = float(np.clip(cen0, cen_lo, cen_hi))

    sigma0 = max(xspan / 12.0, 2.0)
    sigma_lo = 0.5
    sigma_hi = max(xspan, 1.0)
    if sigma_bounds is not None and len(sigma_bounds) >= 2:
        sigma_lo = max(sigma_lo, float(sigma_bounds[0]))
        sigma_hi = min(sigma_hi, float(sigma_bounds[1]))
        if sigma_hi <= sigma_lo:
            result['message'] = 'invalid_sigma_bounds'
            return result
    sigma0 = float(np.clip(sigma0, sigma_lo, sigma_hi))

    p0 = [amp0, cen0, sigma0, m0, b0]
    bounds = (
        [0.0, cen_lo, sigma_lo, -np.inf, -np.inf],
        [np.inf, cen_hi, sigma_hi, np.inf, np.inf],
    )

    try:
        popt, pcov = curve_fit(single_gaussian_with_linear_bg, xf, yf, p0=p0, bounds=bounds, maxfev=50000)
        amp, cen, sigma, m, b = [float(v) for v in popt]
        perr = np.sqrt(np.diag(pcov)) if pcov is not None else np.full(5, np.nan)
        cen_se = float(perr[1]) if len(perr) > 1 else np.nan

        if strict_local_peak and candidate is not None:
            tol = float(max(peak_match_tolerance, 2.0))
            if abs(cen - float(candidate['center'])) > tol:
                result['message'] = f'fitted_center_drifted:{cen:.3f}'
                return result
        if reject_boundary_fit:
            boundary_margin = float(max(edge_margin, 2.0))
            if (cen >= cen_hi - boundary_margin) or (cen <= cen_lo + boundary_margin):
                result['message'] = f'boundary_fit_rejected:{cen:.3f}'
                return result

        bg = m * xf + b
        peak1 = gaussian(xf, amp, cen, sigma)
        yhat = bg + peak1

        stats = fit_stats(yf, yhat, k=5)
        area = amp * sigma * math.sqrt(2.0 * math.pi)
        fwhm = 2.0 * math.sqrt(2.0 * math.log(2.0)) * sigma

        result.update({
            'success': True,
            'center': cen,
            'center_se': cen_se,
            'center1': cen,
            'center1_se': cen_se,
            'effective_center': cen,
            'effective_center_se': cen_se,
            'midpoint_center': cen,
            'fwhm': fwhm,
            'fwhm1': fwhm,
            'sigma': sigma,
            'sigma1': sigma,
            'amplitude': amp,
            'amplitude1': amp,
            'area': area,
            'area1': area,
            'bg_slope': m,
            'bg_intercept': b,
            'fit_yhat': yhat,
            'fit_y_bg': bg,
            'peak1_y': peak1,
            'peak2_y': None,
            **stats,
        })
    except Exception as exc:
        result['message'] = f'single_fit_failed: {exc}'
    return result

def fit_double_gaussian(
    x: np.ndarray,
    y: np.ndarray,
    fit_window: Tuple[float, float],
    exclude_ranges: Optional[Sequence[Tuple[float, float]]] = None,
    min_points: int = 25,
    min_separation: float = 3.0,
    effective_center_method: str = 'area_weighted',
    strategy: str = 'simultaneous_linear',
    init_centers: Optional[Sequence[float]] = None,
    center_bounds: Optional[Sequence[Sequence[float]]] = None,
    init_sigmas: Optional[Sequence[float]] = None,
    sigma_bounds: Optional[Sequence[Sequence[float]]] = None,
    local_bg_fraction: float = 0.12,
    robust_loss: str = 'linear',
    robust_f_scale: float = 1.0,
) -> Dict[str, Any]:
    xf_all, yf_all = prepare_fit_region(x, y, fit_window, exclude_ranges)
    xf, yf = _select_primary_fit_segment(xf_all, yf_all, min_points)
    result: Dict[str, Any] = {
        'success': False,
        'mode_used': 'double',
        'center': np.nan,
        'center_se': np.nan,
        'center1': np.nan,
        'center1_se': np.nan,
        'center2': np.nan,
        'center2_se': np.nan,
        'effective_center': np.nan,
        'effective_center_se': np.nan,
        'midpoint_center': np.nan,
        'fwhm': np.nan,
        'fwhm1': np.nan,
        'fwhm2': np.nan,
        'sigma': np.nan,
        'sigma1': np.nan,
        'sigma2': np.nan,
        'amplitude': np.nan,
        'amplitude1': np.nan,
        'amplitude2': np.nan,
        'area': np.nan,
        'area1': np.nan,
        'area2': np.nan,
        'bg_slope': np.nan,
        'bg_intercept': np.nan,
        'peak_separation': np.nan,
        'n_points': int(len(xf_all)),
        'fit_x': xf,
        'fit_y': yf,
        'fit_yhat': None,
        'fit_y_bg': None,
        'peak1_y': None,
        'peak2_y': None,
        'rss': np.nan,
        'aic': np.nan,
        'bic': np.nan,
        'message': '',
    }
    if len(xf) < min_points:
        result['message'] = 'too_few_points'
        return result

    xlo = float(xf[0])
    xhi = float(xf[-1])
    span = max(xhi - xlo, 1.0)
    m0, b0 = _linear_guess(xf, yf)
    resid0 = yf - (m0 * xf + b0)
    if not np.isfinite(resid0).any():
        resid0 = yf - np.nanmedian(yf)
    i_main = int(np.nanargmax(resid0))
    cen_main = float(xf[i_main])
    amp_main = float(max(resid0[i_main], np.nanpercentile(resid0, 90), 1e-6))
    cen_lo, cen_hi = _fit_center_bounds(xf, fit_window)

    cen1_heur = max(cen_lo, cen_main - 0.08 * span)
    cen2_heur = min(cen_hi, cen_main + 0.08 * span)
    if cen2_heur <= cen1_heur + min_separation:
        cen1_heur = max(cen_lo, cen_main - 0.15 * span)
        cen2_heur = min(cen_hi, cen_main + 0.15 * span)
    if cen2_heur <= cen1_heur + min_separation:
        cen1_heur = cen_lo
        cen2_heur = cen_hi

    cen1_0, cen2_0 = _parse_double_pair(init_centers, (cen1_heur, cen2_heur))
    if cen2_0 < cen1_0:
        cen1_0, cen2_0 = cen2_0, cen1_0

    sigma1_def = max(span / 18.0, 2.0)
    sigma2_def = max(span / 18.0, 2.0)
    sigma1_0, sigma2_0 = _parse_double_pair(init_sigmas, (sigma1_def, sigma2_def))

    (cen1_lb, cen1_ub), (cen2_lb, cen2_ub) = _parse_double_bounds(center_bounds, (cen_lo, cen_hi), (cen_lo, cen_hi))
    (sig1_lb, sig1_ub), (sig2_lb, sig2_ub) = _parse_double_bounds(sigma_bounds, (0.5, span), (0.5, span))

    cen1_lb = max(cen_lo, min(cen1_lb, cen_hi))
    cen1_ub = max(cen1_lb + 1e-6, min(cen1_ub, cen_hi))
    cen2_lb = max(cen_lo, min(cen2_lb, cen_hi))
    cen2_ub = max(cen2_lb + 1e-6, min(cen2_ub, cen_hi))
    sig1_lb = max(0.5, sig1_lb)
    sig2_lb = max(0.5, sig2_lb)
    sig1_ub = max(sig1_lb + 1e-6, sig1_ub)
    sig2_ub = max(sig2_lb + 1e-6, sig2_ub)

    cen1_0 = float(np.clip(cen1_0, cen1_lb, cen1_ub))
    cen2_0 = float(np.clip(cen2_0, cen2_lb, cen2_ub))
    sigma1_0 = float(np.clip(sigma1_0, sig1_lb, sig1_ub))
    sigma2_0 = float(np.clip(sigma2_0, sig2_lb, sig2_ub))

    strategy = str(strategy).lower().strip()

    try:
        if strategy == 'local_bg':
            bg, m, b = _edge_linear_background(xf, yf, fraction=local_bg_fraction)
            y_sub = yf - bg
            amp1_0 = float(max(np.interp(cen1_0, xf, y_sub), np.nanpercentile(y_sub, 90), 1e-6))
            amp2_0 = float(max(np.interp(cen2_0, xf, y_sub), 0.3 * np.nanpercentile(y_sub, 90), 1e-6))

            p0 = np.array([amp1_0, cen1_0, sigma1_0, amp2_0, cen2_0, sigma2_0], dtype=float)
            lb = np.array([0.0, cen1_lb, sig1_lb, 0.0, cen2_lb, sig2_lb], dtype=float)
            ub = np.array([np.inf, cen1_ub, sig1_ub, np.inf, cen2_ub, sig2_ub], dtype=float)

            def residuals_local(params: np.ndarray) -> np.ndarray:
                a1, c1, s1, a2, c2, s2 = params
                return gaussian(xf, a1, c1, s1) + gaussian(xf, a2, c2, s2) - y_sub

            lsq = least_squares(
                residuals_local,
                p0,
                bounds=(lb, ub),
                loss=str(robust_loss),
                f_scale=float(max(robust_f_scale, 1e-6)),
                max_nfev=100000,
            )
            amp1, cen1, sigma1, amp2, cen2, sigma2 = [float(v) for v in lsq.x]
            pcov = _covariance_from_least_squares(lsq.jac, residuals_local(lsq.x), 6)
            perr = np.sqrt(np.diag(pcov)) if pcov is not None else np.full(6, np.nan)
            peak1 = gaussian(xf, amp1, cen1, sigma1)
            peak2 = gaussian(xf, amp2, cen2, sigma2)
            yhat = bg + peak1 + peak2
        else:
            amp1_0 = 0.6 * amp_main
            amp2_0 = 0.4 * amp_main
            p0 = [amp1_0, cen1_0, sigma1_0, amp2_0, cen2_0, sigma2_0, m0, b0]
            bounds = (
                [0.0, cen1_lb, sig1_lb, 0.0, cen2_lb, sig2_lb, -np.inf, -np.inf],
                [np.inf, cen1_ub, sig1_ub, np.inf, cen2_ub, sig2_ub, np.inf, np.inf],
            )
            popt, pcov = curve_fit(double_gaussian_with_linear_bg, xf, yf, p0=p0, bounds=bounds, maxfev=100000)
            amp1, cen1, sigma1, amp2, cen2, sigma2, m, b = [float(v) for v in popt]
            perr = np.sqrt(np.diag(pcov)) if pcov is not None else np.full(8, np.nan)
            peak1 = gaussian(xf, amp1, cen1, sigma1)
            peak2 = gaussian(xf, amp2, cen2, sigma2)
            bg = m * xf + b
            yhat = bg + peak1 + peak2

        if cen2 < cen1:
            amp1, amp2 = amp2, amp1
            cen1, cen2 = cen2, cen1
            sigma1, sigma2 = sigma2, sigma1
            peak1, peak2 = peak2, peak1
            if strategy == 'local_bg':
                perr = np.array([perr[3], perr[4], perr[5], perr[0], perr[1], perr[2]])
            else:
                perr = np.array([perr[3], perr[4], perr[5], perr[0], perr[1], perr[2], perr[6], perr[7]])

        separation = cen2 - cen1
        if separation < min_separation:
            result['message'] = f'double_fit_rejected_small_separation:{separation:.3f}'
            return result

        area1 = amp1 * sigma1 * math.sqrt(2.0 * math.pi)
        area2 = amp2 * sigma2 * math.sqrt(2.0 * math.pi)
        midpoint = 0.5 * (cen1 + cen2)
        method = str(effective_center_method).lower().strip()
        if method == 'midpoint':
            eff = midpoint
        elif method == 'center1':
            eff = cen1
        elif method == 'center2':
            eff = cen2
        else:
            denom = area1 + area2
            eff = midpoint if denom <= 0 else (area1 * cen1 + area2 * cen2) / denom

        k = 6 if strategy == 'local_bg' else 8
        stats = fit_stats(yf, yhat, k=k)
        fwhm1 = 2.0 * math.sqrt(2.0 * math.log(2.0)) * sigma1
        fwhm2 = 2.0 * math.sqrt(2.0 * math.log(2.0)) * sigma2

        if strategy == 'local_bg':
            cen1_se = float(perr[1]) if len(perr) > 1 else np.nan
            cen2_se = float(perr[4]) if len(perr) > 4 else np.nan
        else:
            cen1_se = float(perr[1]) if len(perr) > 1 else np.nan
            cen2_se = float(perr[4]) if len(perr) > 4 else np.nan
        eff_se = np.sqrt(cen1_se**2 + cen2_se**2) / 2.0 if np.isfinite(cen1_se) and np.isfinite(cen2_se) else np.nan

        result.update({
            'success': True,
            'center': eff,
            'center_se': eff_se,
            'center1': cen1,
            'center1_se': cen1_se,
            'center2': cen2,
            'center2_se': cen2_se,
            'effective_center': eff,
            'effective_center_se': eff_se,
            'midpoint_center': midpoint,
            'fwhm1': fwhm1,
            'fwhm2': fwhm2,
            'sigma1': sigma1,
            'sigma2': sigma2,
            'amplitude1': amp1,
            'amplitude2': amp2,
            'area': area1 + area2,
            'area1': area1,
            'area2': area2,
            'bg_slope': float(m),
            'bg_intercept': float(b),
            'peak_separation': separation,
            'fit_yhat': yhat,
            'fit_y_bg': bg,
            'peak1_y': peak1,
            'peak2_y': peak2,
            **stats,
        })
    except Exception as exc:
        result['message'] = f'double_fit_failed: {exc}'
    return result


def fallback_peak_max(

    x: np.ndarray,
    y: np.ndarray,
    fit_window: Tuple[float, float],
    exclude_ranges: Optional[Sequence[Tuple[float, float]]] = None,
) -> Dict[str, Any]:
    xf_all, yf_all = prepare_fit_region(x, y, fit_window, exclude_ranges)
    xf, yf = _select_primary_fit_segment(xf_all, yf_all, min_points=12)
    result: Dict[str, Any] = {
        'success': False,
        'mode_used': 'max',
        'center': np.nan,
        'center_se': np.nan,
        'center1': np.nan,
        'center1_se': np.nan,
        'center2': np.nan,
        'center2_se': np.nan,
        'effective_center': np.nan,
        'effective_center_se': np.nan,
        'midpoint_center': np.nan,
        'fwhm': np.nan,
        'fwhm1': np.nan,
        'fwhm2': np.nan,
        'sigma': np.nan,
        'sigma1': np.nan,
        'sigma2': np.nan,
        'amplitude': np.nan,
        'amplitude1': np.nan,
        'amplitude2': np.nan,
        'area': np.nan,
        'area1': np.nan,
        'area2': np.nan,
        'bg_slope': np.nan,
        'bg_intercept': np.nan,
        'peak_separation': np.nan,
        'n_points': int(len(xf_all)),
        'fit_x': xf,
        'fit_y': yf,
        'fit_yhat': None,
        'fit_y_bg': None,
        'peak1_y': None,
        'peak2_y': None,
        'rss': np.nan,
        'aic': np.nan,
        'bic': np.nan,
        'message': '',
    }
    if len(xf) == 0:
        result['message'] = 'empty_fit_region'
        return result
    m0, b0 = _linear_guess(xf, yf)
    resid = yf - (m0 * xf + b0)
    if not np.isfinite(resid).any():
        resid = yf
    idx = int(np.nanargmax(resid))
    center = float(xf[idx])
    amp = float(resid[idx])
    result.update({
        'success': True,
        'center': center,
        'center1': center,
        'effective_center': center,
        'midpoint_center': center,
        'amplitude': amp,
        'amplitude1': amp,
    })
    return result


def choose_fit_result(
    single_res: Dict[str, Any],
    double_res: Dict[str, Any],
    auto_cfg: Dict[str, Any],
) -> Dict[str, Any]:
    criterion = auto_cfg.get("criterion", "aic").lower()
    delta_threshold = float(auto_cfg.get("delta_threshold", 6.0))
    require_both = bool(auto_cfg.get("require_both_success", False))

    if single_res.get("success") and not double_res.get("success"):
        return single_res
    if double_res.get("success") and not single_res.get("success"):
        return double_res
    if not single_res.get("success") and not double_res.get("success"):
        return single_res

    if require_both and not (single_res.get("success") and double_res.get("success")):
        return single_res if single_res.get("success") else double_res

    if criterion == "bic":
        ds = float(single_res["bic"] - double_res["bic"])
    else:
        ds = float(single_res["aic"] - double_res["aic"])

    return double_res if ds >= delta_threshold else single_res




def _peak_rule_matches(rule: Dict[str, Any], spec_idx: int, pressure: float, label: str) -> bool:
    idx1 = int(spec_idx) + 1
    if 'index_min' in rule and idx1 < int(rule['index_min']):
        return False
    if 'index_max' in rule and idx1 > int(rule['index_max']):
        return False
    if 'pressure_min' in rule and float(pressure) < float(rule['pressure_min']):
        return False
    if 'pressure_max' in rule and float(pressure) > float(rule['pressure_max']):
        return False
    labels = rule.get('labels')
    if labels is not None and str(label) not in [str(v) for v in labels]:
        return False
    return True


def resolve_peak_fit_cfg(base_cfg: Dict[str, Any], spec_idx: int, pressure: float, label: str) -> Dict[str, Any]:
    cfg = deep_copy_config(base_cfg)
    rules = list(cfg.pop('fit_rules', []))
    for rule in rules:
        if _peak_rule_matches(rule, spec_idx=spec_idx, pressure=pressure, label=label):
            cfg = deep_update(cfg, rule.get('update', {}))
    return cfg


def run_peak_fit(spec: Spectrum, fit_cfg: Dict[str, Any]) -> Dict[str, Any]:
    mode = str(fit_cfg.get("mode", "single")).lower()
    fit_window = tuple(fit_cfg["fit_window"])
    exclude_ranges = fit_cfg.get("exclude_ranges", [])
    min_points = int(fit_cfg.get("min_points", 20))
    fallback = bool(fit_cfg.get("fallback_to_max_if_fit_fails", True))

    if mode == "none":
        return {"success": False, "mode_used": "none", "label": spec.label, "pressure": spec.pressure}

    if mode == "single":
        res = fit_single_gaussian(
            spec.x,
            spec.y_corrected,
            fit_window,
            exclude_ranges,
            min_points=min_points,
            search_window=tuple(fit_cfg.get("single_search_window")) if fit_cfg.get("single_search_window") is not None else None,
            center_bounds=fit_cfg.get("single_center_bounds"),
            sigma_bounds=fit_cfg.get("single_sigma_bounds"),
            strict_local_peak=bool(fit_cfg.get("single_strict_local_peak", False)),
            edge_margin=float(fit_cfg.get("single_edge_margin", 8.0)),
            smooth_window=int(fit_cfg.get("single_savgol_window", 11)),
            smooth_polyorder=int(fit_cfg.get("single_savgol_polyorder", 2)),
            min_prominence_ratio=float(fit_cfg.get("single_min_prominence_ratio", 0.08)),
            local_bg_fraction=float(fit_cfg.get("single_local_bg_fraction", 0.12)),
            center_half_width=float(fit_cfg.get("single_center_half_width", 18.0)),
            peak_match_tolerance=float(fit_cfg.get("single_peak_match_tolerance", 14.0)),
            reject_boundary_fit=bool(fit_cfg.get("single_reject_boundary_fit", False)),
        )
    elif mode == "double":
        res = fit_double_gaussian(
            spec.x,
            spec.y_corrected,
            fit_window,
            exclude_ranges,
            min_points=max(min_points, 25),
            min_separation=float(fit_cfg.get("double_min_separation", 3.0)),
            effective_center_method=str(fit_cfg.get("effective_center_method", "area_weighted")),
            strategy=str(fit_cfg.get("double_strategy", "simultaneous_linear")),
            init_centers=fit_cfg.get("double_init_centers"),
            center_bounds=fit_cfg.get("double_center_bounds"),
            init_sigmas=fit_cfg.get("double_init_sigmas"),
            sigma_bounds=fit_cfg.get("double_sigma_bounds"),
            local_bg_fraction=float(fit_cfg.get("double_local_bg_fraction", 0.12)),
            robust_loss=str(fit_cfg.get("double_loss", "linear")),
            robust_f_scale=float(fit_cfg.get("double_f_scale", 1.0)),
        )
    elif mode == "auto":
        single_res = fit_single_gaussian(
            spec.x,
            spec.y_corrected,
            fit_window,
            exclude_ranges,
            min_points=min_points,
            search_window=tuple(fit_cfg.get("single_search_window")) if fit_cfg.get("single_search_window") is not None else None,
            center_bounds=fit_cfg.get("single_center_bounds"),
            sigma_bounds=fit_cfg.get("single_sigma_bounds"),
            strict_local_peak=bool(fit_cfg.get("single_strict_local_peak", False)),
            edge_margin=float(fit_cfg.get("single_edge_margin", 8.0)),
            smooth_window=int(fit_cfg.get("single_savgol_window", 11)),
            smooth_polyorder=int(fit_cfg.get("single_savgol_polyorder", 2)),
            min_prominence_ratio=float(fit_cfg.get("single_min_prominence_ratio", 0.08)),
            local_bg_fraction=float(fit_cfg.get("single_local_bg_fraction", 0.12)),
            center_half_width=float(fit_cfg.get("single_center_half_width", 18.0)),
            peak_match_tolerance=float(fit_cfg.get("single_peak_match_tolerance", 14.0)),
            reject_boundary_fit=bool(fit_cfg.get("single_reject_boundary_fit", False)),
        )
        double_res = fit_double_gaussian(
            spec.x,
            spec.y_corrected,
            fit_window,
            exclude_ranges,
            min_points=max(min_points, 25),
            min_separation=float(fit_cfg.get("double_min_separation", 3.0)),
            effective_center_method=str(fit_cfg.get("effective_center_method", "area_weighted")),
            strategy=str(fit_cfg.get("double_strategy", "simultaneous_linear")),
            init_centers=fit_cfg.get("double_init_centers"),
            center_bounds=fit_cfg.get("double_center_bounds"),
            init_sigmas=fit_cfg.get("double_init_sigmas"),
            sigma_bounds=fit_cfg.get("double_sigma_bounds"),
            local_bg_fraction=float(fit_cfg.get("double_local_bg_fraction", 0.12)),
            robust_loss=str(fit_cfg.get("double_loss", "linear")),
            robust_f_scale=float(fit_cfg.get("double_f_scale", 1.0)),
        )
        res = choose_fit_result(single_res, double_res, fit_cfg.get("auto_criteria", {}))
    else:
        raise ValueError(f"Unsupported fit mode: {mode}")

    if (not res.get("success", False)) and fallback:
        fallback_res = fallback_peak_max(spec.x, spec.y_corrected, fit_window, exclude_ranges)
        fallback_res["message"] = (res.get("message", "") + " | fallback_to_max").strip(" |")
        res = fallback_res

    res = copy.deepcopy(res)
    res["label"] = spec.label
    res["pressure"] = spec.pressure
    if "fit_x" in res:
        res["fit_x"] = np.asarray(res["fit_x"], dtype=float)
    if "fit_y" in res:
        res["fit_y"] = np.asarray(res["fit_y"], dtype=float)
    if res.get("fit_yhat") is not None:
        res["fit_yhat"] = np.asarray(res["fit_yhat"], dtype=float)
    return res


def sort_spectra(spectra: List[Spectrum], sort_by_pressure: bool) -> List[Spectrum]:
    return sorted(spectra, key=lambda s: s.pressure) if sort_by_pressure else list(spectra)




def _series_edges(values: np.ndarray) -> np.ndarray:
    vals = np.asarray(values, dtype=float)
    if vals.ndim != 1 or len(vals) == 0:
        raise ValueError("values must be a non-empty 1D array")
    if len(vals) == 1:
        v = float(vals[0])
        return np.array([v - 0.5, v + 0.5], dtype=float)
    diffs = np.diff(vals)
    edges = np.empty(len(vals) + 1, dtype=float)
    edges[1:-1] = vals[:-1] + 0.5 * diffs
    edges[0] = vals[0] - 0.5 * diffs[0]
    edges[-1] = vals[-1] + 0.5 * diffs[-1]
    return edges


def _clean_axes(ax: plt.Axes) -> None:
    ax.grid(False)
    ax.minorticks_on()
    ax.tick_params(which="major", direction="in", length=5, width=1.0)
    ax.tick_params(which="minor", direction="in", length=3, width=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def _unweighted_linear_fit(x: np.ndarray, y: np.ndarray) -> Optional[np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return None
    x_use = x[mask]
    y_use = y[mask]
    if np.unique(x_use).size < 2:
        return None
    try:
        return np.polyfit(x_use, y_use, 1)
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        return None


def _weighted_linear_fit(x: np.ndarray, y: np.ndarray, yerr: np.ndarray) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    yerr = np.asarray(yerr, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    if np.count_nonzero(mask) < 3:
        return None
    x_use = x[mask]
    y_use = y[mask]
    yerr_use = yerr[mask]
    if np.unique(x_use).size < 2:
        return None
    try:
        coeff, cov = np.polyfit(x_use, y_use, 1, w=1.0 / yerr_use, cov=True)
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        return None
    if not np.all(np.isfinite(coeff)) or not np.all(np.isfinite(cov)):
        return None
    return coeff, cov


def _build_trendline_segment_masks(x: np.ndarray, segment_cfgs: Optional[Sequence[Dict[str, Any]]]) -> List[Tuple[np.ndarray, Dict[str, Any]]]:
    x = np.asarray(x, dtype=float)
    if not segment_cfgs:
        return [(np.isfinite(x), {})]

    segments: List[Tuple[np.ndarray, Dict[str, Any]]] = []
    for raw_cfg in segment_cfgs:
        cfg = dict(raw_cfg) if isinstance(raw_cfg, dict) else {}
        mask = np.isfinite(x)
        if cfg.get('pressure_min') is not None:
            mask &= (x >= float(cfg['pressure_min']))
        if cfg.get('pressure_max') is not None:
            mask &= (x <= float(cfg['pressure_max']))
        segments.append((mask, cfg))
    return segments


def _resolve_mask_ranges_for_item(mask_ranges: Any, spec: Spectrum, idx: int) -> List[Tuple[float, float]]:
    if mask_ranges is None:
        return []
    if isinstance(mask_ranges, dict):
        if spec.label in mask_ranges:
            return _normalize_mask_ranges(mask_ranges.get(spec.label))
        if idx in mask_ranges:
            return _normalize_mask_ranges(mask_ranges.get(idx))
        return []
    if isinstance(mask_ranges, (list, tuple)) and len(mask_ranges) > 0:
        first = mask_ranges[0]
        if isinstance(first, (list, tuple)) and len(first) == 2 and np.isscalar(first[0]) and np.isscalar(first[1]):
            return _normalize_mask_ranges(mask_ranges)
    return []


def plot_waterfall(
    spectra: List[Spectrum],
    outpath: Path,
    window: Optional[Tuple[float, float]],
    mask_ranges: Any,
    vertical_offset: float,
    linewidth: float,
    sort_by_pressure: bool,
    dpi: int,
    line_color: str = 'black',
    normalize_each: bool = False,
    figsize: Tuple[float, float] = (8.6, 12.0),
    label_fontsize: float = 8.0,
    right_label_offset_frac: float = 0.025,
    right_label_inside: bool = True,
    right_label_pad_frac: float = 0.012,
    pressure_label_mode: str = 'all',
    pressure_label_every: int = 3,
    pressure_label_include_first_last: bool = True,
    hide_y_ticks: bool = True,
    title: str = 'Waterfall diagram (ALS BG-subtracted)',
) -> None:
    items = sort_spectra(spectra, sort_by_pressure)

    prepared: List[Tuple[Spectrum, np.ndarray, np.ndarray]] = []
    robust_ranges: List[float] = []
    for idx, spec in enumerate(items):
        item_masks = _resolve_mask_ranges_for_item(mask_ranges, spec, idx)
        mask = in_window(spec.x, window)
        x = spec.x[mask]
        y = spec.y_corrected[mask]
        y = apply_nan_masks(x, y, item_masks)
        finite = np.isfinite(y)
        if finite.sum() == 0:
            continue

        yp = y.copy()
        if normalize_each:
            ymin = np.nanmin(yp)
            ymax = np.nanmax(yp)
            if np.isfinite(ymax - ymin) and ymax > ymin:
                yp = (yp - ymin) / (ymax - ymin)
            else:
                yp[:] = 0.0

        prepared.append((spec, x, yp))

        vals = yp[np.isfinite(yp)]
        if vals.size >= 8:
            q10, q90 = np.nanpercentile(vals, [10.0, 90.0])
            rr = float(q90 - q10)
            if np.isfinite(rr) and rr > 0:
                robust_ranges.append(rr)

    if not prepared:
        return

    if normalize_each:
        step = float(vertical_offset)
    else:
        ref_range = float(np.nanmedian(robust_ranges)) if robust_ranges else 1.0
        if not np.isfinite(ref_range) or ref_range <= 0:
            ref_range = 1.0
        step = float(vertical_offset) * ref_range

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    offset = 0.0
    pressure_label_mode = str(pressure_label_mode).lower()
    pressure_label_every = max(int(pressure_label_every), 1)
    n_prepared = len(prepared)

    for idx, (spec, x, yp) in enumerate(prepared):
        ax.plot(x, yp + offset, lw=linewidth, color=line_color)

        xr = float(np.nanmax(x)) if len(x) else 0.0
        xspan = (window[1] - window[0]) if window else max(float(np.nanmax(x) - np.nanmin(x)), 1.0)
        finite = np.isfinite(yp)
        if np.any(finite):
            vals = yp[finite]
            n_tail = min(12, len(vals))
            y_end = float(np.nanmedian(vals[-n_tail:]) + offset)
        else:
            y_end = float(offset)

        should_label = False
        if pressure_label_mode == 'all':
            should_label = True
        elif pressure_label_mode == 'sparse':
            should_label = (idx % pressure_label_every == 0)
            if pressure_label_include_first_last and idx in {0, n_prepared - 1}:
                should_label = True
        elif pressure_label_mode == 'none':
            should_label = False

        if should_label:
            if right_label_inside:
                if window is not None:
                    x_text = float(window[1]) - right_label_pad_frac * xspan
                else:
                    x_text = xr - right_label_pad_frac * xspan
                ha = 'right'
            else:
                x_text = xr + right_label_offset_frac * xspan
                ha = 'left'

            ax.text(x_text, y_end, f'{spec.pressure:.3f} GPa', ha=ha, va='center', fontsize=label_fontsize)
        offset += step

    ax.set_xlabel('Raman shift (cm$^{-1}$)')
    ax.set_ylabel('Offset BG-subtracted intensity')
    ax.set_title(title)
    if hide_y_ticks:
        ax.set_yticks([])
        ax.set_yticklabels([])
    _clean_axes(ax)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def plot_heatmap(
    xgrid: np.ndarray,
    pressures: np.ndarray,
    zmat: np.ndarray,
    outpath: Path,
    dpi: int,
    title: str = "Heatmap (ALS BG-subtracted)",
) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    xedges = _series_edges(xgrid)
    yedges = _series_edges(pressures)
    im = ax.pcolormesh(
        xedges,
        yedges,
        zmat,
        shading="flat",
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Normalized intensity")
    ax.set_xlabel("Raman shift (cm$^{-1}$)")
    ax.set_ylabel("Pressure (GPa)")
    ax.set_title(title)
    _clean_axes(ax)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)




def integrate_spectrum_window(
    x: np.ndarray,
    y: np.ndarray,
    window: Sequence[float],
    mask_ranges: Optional[Sequence[Tuple[float, float]]] = None,
    area_mode: str = 'signed',
) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(window) != 2:
        raise ValueError('Integration window must contain exactly two values.')
    lo, hi = float(window[0]), float(window[1])
    if hi < lo:
        lo, hi = hi, lo
    sel = (x >= lo) & (x <= hi)
    if not np.any(sel):
        return np.nan
    xs = x[sel]
    ys = y[sel]
    ys = apply_nan_masks(xs, ys, mask_ranges)
    mode = str(area_mode).strip().lower()
    if mode in {'positive', 'positive_only', 'clip_negative'}:
        ys = np.where(np.isfinite(ys), np.maximum(ys, 0.0), ys)
    elif mode in {'absolute', 'abs'}:
        ys = np.where(np.isfinite(ys), np.abs(ys), ys)
    elif mode not in {'signed', 'raw'}:
        raise ValueError(f'Unsupported integration area_mode: {area_mode}')
    finite = np.isfinite(xs) & np.isfinite(ys)
    if finite.sum() < 2:
        return np.nan
    segments = _split_contiguous_segments(xs[finite], ys[finite])
    total = 0.0
    used = False
    for xseg, yseg in segments:
        if len(xseg) < 2:
            continue
        total += float(np.trapz(yseg, xseg))
        used = True
    return total if used else np.nan



def compute_spectral_weight_tables(
    spectra: List[Spectrum],
    spectral_cfg: Dict[str, Any],
    cfg: Dict[str, Any],
    spec_order: Dict[str, int],
    sort_by_pressure: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    windows_cfg = spectral_cfg.get('windows', [])
    if not windows_cfg:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    items = sort_spectra(spectra, sort_by_pressure)
    intensity_rows: List[Dict[str, Any]] = []
    fraction_rows: List[Dict[str, Any]] = []
    ratio_rows: List[Dict[str, Any]] = []

    window_names: List[str] = []
    window_ranges: Dict[str, Tuple[float, float]] = {}
    for item in windows_cfg:
        name = str(item.get('name', '')).strip()
        rng = item.get('range') or item.get('window')
        if not name or rng is None or len(rng) != 2:
            continue
        lo, hi = float(rng[0]), float(rng[1])
        if hi < lo:
            lo, hi = hi, lo
        window_names.append(name)
        window_ranges[name] = (lo, hi)

    if not window_names:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    normalize_to_sum = bool(spectral_cfg.get('normalize_to_sum_of_windows', True))
    ratio_pairs = spectral_cfg.get('ratio_pairs', [])
    area_mode = str(spectral_cfg.get('area_mode', 'signed')).strip().lower()

    for spec in items:
        spec_idx = spec_order[spec.label]
        integration_masks = get_target_mask_ranges(cfg, 'integration', spec_idx, spec.pressure, spec.label)
        base_row: Dict[str, Any] = {'label': spec.label, 'pressure': float(spec.pressure)}
        values: Dict[str, float] = {}
        for name in window_names:
            lo, hi = window_ranges[name]
            val = integrate_spectrum_window(
                spec.x,
                spec.y_corrected,
                (lo, hi),
                integration_masks,
                area_mode=area_mode,
            )
            values[name] = float(val) if np.isfinite(val) else np.nan
            base_row[f'S_{name}'] = values[name]
        intensity_rows.append(base_row)

        frac_row: Dict[str, Any] = {'label': spec.label, 'pressure': float(spec.pressure)}
        total = 0.0
        have_any = False
        for name in window_names:
            val = values[name]
            if np.isfinite(val):
                total += float(val)
                have_any = True
        if (not have_any) or (not normalize_to_sum) or (not np.isfinite(total)) or abs(total) <= 0.0:
            for name in window_names:
                frac_row[f'f_{name}'] = np.nan
        else:
            for name in window_names:
                val = values[name]
                frac_row[f'f_{name}'] = (float(val) / total) if np.isfinite(val) else np.nan
        frac_row['S_total_windows'] = total if have_any and np.isfinite(total) else np.nan
        fraction_rows.append(frac_row)

        ratio_row: Dict[str, Any] = {'label': spec.label, 'pressure': float(spec.pressure)}
        for pair in ratio_pairs:
            pair_name = str(pair.get('name', '')).strip()
            numerator = str(pair.get('numerator', '')).strip()
            denominator = str(pair.get('denominator', '')).strip()
            if not pair_name or numerator not in values or denominator not in values:
                continue
            num = values.get(numerator, np.nan)
            den = values.get(denominator, np.nan)
            if (not np.isfinite(num)) or (not np.isfinite(den)) or den == 0:
                ratio_row[pair_name] = np.nan
            else:
                ratio_row[pair_name] = float(num) / float(den)
        ratio_rows.append(ratio_row)

    df_intensity = pd.DataFrame(intensity_rows)
    df_fraction = pd.DataFrame(fraction_rows)
    df_ratio = pd.DataFrame(ratio_rows)
    return df_intensity, df_fraction, df_ratio



def _plot_series_vs_pressure(
    df: pd.DataFrame,
    prefix: Optional[str],
    ylabel: str,
    title: str,
    outpath: Path,
    dpi: int,
) -> None:
    if df.empty or 'pressure' not in df.columns:
        return
    all_numeric_cols = [c for c in df.columns if c not in {'label', 'pressure'} and np.isfinite(pd.to_numeric(df[c], errors='coerce')).any()]
    if prefix is None:
        series_cols = all_numeric_cols
    elif prefix == '':
        series_cols = all_numeric_cols
    else:
        series_cols = [c for c in all_numeric_cols if c.startswith(prefix)]
    if not series_cols:
        return
    work = df.copy().sort_values('pressure')
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for col in series_cols:
        y = pd.to_numeric(work[col], errors='coerce').to_numpy(float)
        x = pd.to_numeric(work['pressure'], errors='coerce').to_numpy(float)
        finite = np.isfinite(x) & np.isfinite(y)
        if not np.any(finite):
            continue
        label = col[len(prefix):] if col.startswith(prefix) else col
        ax.plot(x[finite], y[finite], marker='o', markersize=4.4, linewidth=1.5, label=label)
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    _clean_axes(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)



def plot_integrated_intensities(df: pd.DataFrame, outpath: Path, dpi: int) -> None:
    _plot_series_vs_pressure(
        df=df,
        prefix='S_',
        ylabel='Integrated intensity (a.u.·cm$^{-1}$)',
        title='Integrated spectral weight vs pressure',
        outpath=outpath,
        dpi=dpi,
    )



def plot_spectral_fractions(df: pd.DataFrame, outpath: Path, dpi: int) -> None:
    _plot_series_vs_pressure(
        df=df,
        prefix='f_',
        ylabel='Fraction of total integrated intensity',
        title='Fractional spectral weight vs pressure',
        outpath=outpath,
        dpi=dpi,
    )



def plot_spectral_ratios(df: pd.DataFrame, outpath: Path, dpi: int) -> None:
    _plot_series_vs_pressure(
        df=df,
        prefix='',
        ylabel='Integrated-intensity ratio',
        title='Spectral-weight ratio vs pressure',
        outpath=outpath,
        dpi=dpi,
    )



def plot_pplot_center(df_peak: pd.DataFrame, outpath: Path, dpi: int, plot_cfg: Optional[Dict[str, Any]] = None) -> None:
    good = df_peak[df_peak["success"] == True].copy()  # noqa: E712
    if good.empty:
        return

    plot_cfg = dict(plot_cfg or {})
    segment_cfgs = plot_cfg.get('trendline_segments', [])
    weighted_lw = float(plot_cfg.get('weighted_linewidth', 1.7))
    unweighted_lw = float(plot_cfg.get('unweighted_linewidth', 1.5))
    weighted_ls = str(plot_cfg.get('weighted_linestyle', '-'))
    unweighted_ls = str(plot_cfg.get('unweighted_linestyle', '--'))
    fig_w = float(plot_cfg.get('figsize', (7.2, 5.0))[0]) if isinstance(plot_cfg.get('figsize', (7.2, 5.0)), (list, tuple)) else 7.2
    fig_h = float(plot_cfg.get('figsize', (7.2, 5.0))[1]) if isinstance(plot_cfg.get('figsize', (7.2, 5.0)), (list, tuple)) else 5.0

    good = good.sort_values("pressure")
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    styles = {
        "center": {"color": "#4C72B0", "marker": "D", "label": "Center"},
        "center1": {"color": "#1f77b4", "marker": "o", "label": "Peak 1"},
        "center2": {"color": "#C44E52", "marker": "s", "label": "Peak 2"},
    }

    series_order = []
    if "center1" in good.columns and np.isfinite(good["center1"]).any():
        series_order.append("center1")
    if "center2" in good.columns and np.isfinite(good["center2"]).any():
        series_order.append("center2")
    if not series_order and "center" in good.columns and np.isfinite(good["center"]).any():
        series_order.append("center")

    for key in series_order:
        y = good[key].to_numpy(float)
        yerr_col = f"{key}_se" if f"{key}_se" in good.columns else ("center_se" if key == "center" and "center_se" in good.columns else None)
        yerr = good[yerr_col].to_numpy(float) if yerr_col is not None else np.full_like(y, np.nan, dtype=float)
        finite_err = np.isfinite(yerr) & (yerr > 0)
        x = good["pressure"].to_numpy(float)
        st = styles[key]

        if np.any(finite_err):
            ax.errorbar(
                x[finite_err], y[finite_err], yerr=yerr[finite_err],
                fmt=st["marker"], linestyle="none", markersize=5.6,
                markerfacecolor="white", markeredgecolor=st["color"], markeredgewidth=1.2,
                ecolor=st["color"], elinewidth=1.0, capsize=2.5, label=st["label"],
            )
        if np.any(~finite_err & np.isfinite(y)):
            ax.plot(
                x[~finite_err & np.isfinite(y)], y[~finite_err & np.isfinite(y)],
                linestyle="none", marker=st["marker"], markersize=5.2,
                markerfacecolor=st["color"], markeredgecolor="black", markeredgewidth=0.3, alpha=0.85,
            )

        for segment_mask, segment_cfg in _build_trendline_segment_masks(x, segment_cfgs):
            seg_mask = segment_mask & np.isfinite(x) & np.isfinite(y)
            if np.count_nonzero(seg_mask) < 2:
                continue
            x_seg = x[seg_mask]
            y_seg = y[seg_mask]
            yerr_seg = yerr[seg_mask]
            finite_err_seg = np.isfinite(yerr_seg) & (yerr_seg > 0)

            fit = _weighted_linear_fit(x_seg, y_seg, yerr_seg) if np.count_nonzero(finite_err_seg) >= 3 else None
            if fit is not None:
                coeff, _ = fit
                xx = np.linspace(np.nanmin(x_seg), np.nanmax(x_seg), 200)
                yy = coeff[0] * xx + coeff[1]
                ax.plot(
                    xx, yy,
                    color=st["color"],
                    lw=float(segment_cfg.get('linewidth', weighted_lw)),
                    ls=str(segment_cfg.get('linestyle', weighted_ls)),
                    alpha=float(segment_cfg.get('alpha', 0.95)),
                )
            else:
                coeff = _unweighted_linear_fit(x_seg, y_seg)
                if coeff is not None:
                    xx = np.linspace(np.nanmin(x_seg), np.nanmax(x_seg), 200)
                    yy = coeff[0] * xx + coeff[1]
                    ax.plot(
                        xx, yy,
                        color=st["color"],
                        lw=float(segment_cfg.get('linewidth', unweighted_lw)),
                        ls=str(segment_cfg.get('linestyle', unweighted_ls)),
                        alpha=float(segment_cfg.get('alpha', 0.9)),
                    )

    ax.set_xlabel("Pressure (GPa)")
    ax.set_ylabel("Peak center (cm$^{-1}$)")
    ax.set_title(str(plot_cfg.get('title', "Peak center vs pressure")))
    _clean_axes(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)



def plot_pplot_width(df_peak: pd.DataFrame, outpath: Path, dpi: int, plot_cfg: Optional[Dict[str, Any]] = None) -> None:
    good = df_peak[df_peak["success"] == True].copy()  # noqa: E712
    if good.empty:
        return

    plot_cfg = dict(plot_cfg or {})
    segment_cfgs = plot_cfg.get('trendline_segments', [])
    line_lw = float(plot_cfg.get('linewidth', 1.6))
    line_ls = str(plot_cfg.get('linestyle', '-'))
    fig_w = float(plot_cfg.get('figsize', (7.2, 5.0))[0]) if isinstance(plot_cfg.get('figsize', (7.2, 5.0)), (list, tuple)) else 7.2
    fig_h = float(plot_cfg.get('figsize', (7.2, 5.0))[1]) if isinstance(plot_cfg.get('figsize', (7.2, 5.0)), (list, tuple)) else 5.0

    good = good.sort_values("pressure")
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    styles = {
        "fwhm": {"color": "#4C72B0", "marker": "D", "label": "FWHM"},
        "fwhm1": {"color": "#1f77b4", "marker": "o", "label": "FWHM 1"},
        "fwhm2": {"color": "#C44E52", "marker": "s", "label": "FWHM 2"},
    }

    series_order = []
    if "fwhm1" in good.columns and np.isfinite(good["fwhm1"]).any():
        series_order.append("fwhm1")
    if "fwhm2" in good.columns and np.isfinite(good["fwhm2"]).any():
        series_order.append("fwhm2")
    if not series_order and "fwhm" in good.columns and np.isfinite(good["fwhm"]).any():
        series_order.append("fwhm")

    x = good["pressure"].to_numpy(float)
    for key in series_order:
        y = good[key].to_numpy(float)
        finite = np.isfinite(y)
        if not np.any(finite):
            continue
        st = styles[key]
        ax.plot(
            x[finite], y[finite], linestyle="none", marker=st["marker"], markersize=5.4,
            markerfacecolor="white", markeredgecolor=st["color"], markeredgewidth=1.2, label=st["label"],
        )

        for segment_mask, segment_cfg in _build_trendline_segment_masks(x, segment_cfgs):
            seg_mask = segment_mask & np.isfinite(x) & np.isfinite(y)
            if np.count_nonzero(seg_mask) < 2:
                continue
            x_seg = x[seg_mask]
            y_seg = y[seg_mask]
            coeff = _unweighted_linear_fit(x_seg, y_seg)
            if coeff is not None:
                xx = np.linspace(np.nanmin(x_seg), np.nanmax(x_seg), 200)
                yy = coeff[0] * xx + coeff[1]
                ax.plot(
                    xx, yy,
                    color=st["color"],
                    lw=float(segment_cfg.get('linewidth', line_lw)),
                    ls=str(segment_cfg.get('linestyle', line_ls)),
                    alpha=float(segment_cfg.get('alpha', 0.95)),
                )

    ax.set_xlabel("Pressure (GPa)")
    ax.set_ylabel("FWHM (cm$^{-1}$)")
    ax.set_title(str(plot_cfg.get('title', "FWHM vs pressure")))
    _clean_axes(ax)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)



def plot_representative(
    spec: Spectrum | Sequence[Spectrum],
    peak_result: Optional[Dict[str, Any]] | Sequence[Optional[Dict[str, Any]]],
    outpath: Path,
    analysis_window: Optional[Tuple[float, float]],
    fit_cfg: Dict[str, Any],
    dpi: int,
) -> None:
    specs = list(spec) if isinstance(spec, (list, tuple)) else [spec]
    peak_results = list(peak_result) if isinstance(peak_result, (list, tuple)) else [peak_result]
    if len(peak_results) < len(specs):
        peak_results = peak_results + [None] * (len(specs) - len(peak_results))

    nrows = len(specs)
    fig, axes = plt.subplots(nrows, 2, figsize=(9.4, max(3.0 * nrows, 4.8)), squeeze=False, sharex=False)

    for i, (spec_i, peak_i) in enumerate(zip(specs, peak_results)):
        ax0, ax1 = axes[i, 0], axes[i, 1]
        mask = in_window(spec_i.x, analysis_window)
        xplot = spec_i.x[mask]
        yraw = spec_i.y_rate_norm[mask]
        ybase = spec_i.baseline[mask]

        ax0.plot(xplot, yraw, lw=1.0, color="black", label="Raw spectrum")
        ax0.plot(xplot, ybase, lw=1.2, color="#d62728", label="ALS baseline")
        ax0.set_ylabel("Intensity (a.u.)")
        ax0.text(0.02, 0.93, f"P = {spec_i.pressure:.2f} GPa", transform=ax0.transAxes, ha="left", va="top")
        if i == 0:
            ax0.set_title("Baseline correction")
            ax0.legend(frameon=False, loc="upper right")
        _clean_axes(ax0)

        display_masks = peak_i.get('display_mask_ranges', fit_cfg.get('display_mask_ranges', [])) if peak_i is not None else fit_cfg.get('display_mask_ranges', [])
        xmask = spec_i.x[mask]
        ycorr_mask = spec_i.y_corrected[mask]
        ycorr_plot = apply_nan_masks(xmask, ycorr_mask, display_masks)

        fit_x = peak_i.get("fit_x") if peak_i is not None else None
        fit_bg = peak_i.get("fit_y_bg") if peak_i is not None else None
        fit_yhat = peak_i.get("fit_yhat") if peak_i is not None else None
        peak1 = peak_i.get("peak1_y") if peak_i is not None else None
        peak2 = peak_i.get("peak2_y") if peak_i is not None else None

        plotted_local_bg_removed = False
        if fit_x is not None and fit_bg is not None:
            fit_x = np.asarray(fit_x, dtype=float)
            fit_bg = np.asarray(fit_bg, dtype=float)
            fit_bg_interp = np.full_like(xmask, np.nan, dtype=float)
            if fit_x.size == 1:
                fit_bg_interp[np.isclose(xmask, fit_x[0])] = fit_bg[0]
            elif fit_x.size >= 2:
                domain = (xmask >= np.nanmin(fit_x)) & (xmask <= np.nanmax(fit_x))
                fit_bg_interp[domain] = np.interp(xmask[domain], fit_x, fit_bg)
            y_local_bg_removed = ycorr_mask - fit_bg_interp
            y_local_bg_removed = apply_nan_masks(xmask, y_local_bg_removed, display_masks)
            if np.isfinite(y_local_bg_removed).any():
                ax1.plot(xmask, y_local_bg_removed, lw=1.0, color="black", label="Spectrum - fitted bg")
                plotted_local_bg_removed = True

        if not plotted_local_bg_removed:
            ax1.plot(xmask, ycorr_plot, lw=1.0, color="black", label="Corrected spectrum")

        ax1.axhline(0.0, color="0.75", lw=0.8)

        if peak_i is not None and peak_i.get("success", False):
            if fit_x is not None and fit_yhat is not None and fit_bg is not None:
                total_minus_bg = np.asarray(fit_yhat, dtype=float) - np.asarray(fit_bg, dtype=float)
                ax1.plot(fit_x, total_minus_bg, lw=1.4, color="#d62728", label="Peak 1 + Peak 2")
            if fit_x is not None and peak1 is not None:
                ax1.plot(fit_x, peak1, lw=1.2, ls="--", color="#1f77b4", label="Peak 1")
            if fit_x is not None and peak2 is not None:
                ax1.plot(fit_x, peak2, lw=1.2, ls="--", color="#2ca02c", label="Peak 2")

            metric = peak_i.get("r2")
            metric_text = f"R$^2$ = {metric:.3f}" if metric is not None and np.isfinite(metric) else str(peak_i.get("mode_used", "fit"))
            quality = "reliable" if peak_i.get("success", False) else "failed"
            ax1.text(0.02, 0.93, f"{metric_text} | {quality}", transform=ax1.transAxes, ha="left", va="top")
        if i == 0:
            ax1.set_title("Representative fit (local bg removed in fit window)")
            ax1.legend(frameon=False, loc="upper right")
        _clean_axes(ax1)

    axes[-1, 0].set_xlabel("Raman shift (cm$^{-1}$)")
    axes[-1, 1].set_xlabel("Raman shift (cm$^{-1}$)")
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)



def spectra_to_wide_xy(spectra: List[Spectrum], value_attr: str) -> pd.DataFrame:
    max_len = max(len(s.x) for s in spectra)
    cols: List[str] = []
    data: Dict[str, np.ndarray] = {}
    for s in spectra:
        xcol = np.full(max_len, np.nan)
        ycol = np.full(max_len, np.nan)
        xv = s.x
        yv = getattr(s, value_attr)
        xcol[: len(xv)] = xv
        ycol[: len(yv)] = yv
        data[f"{s.label}_x"] = xcol
        data[f"{s.label}_y"] = ycol
        cols.extend([f"{s.label}_x", f"{s.label}_y"])
    return pd.DataFrame(data)[cols]


def serializable_value(v: Any) -> Any:
    if isinstance(v, (dict, list, tuple)):
        return json.dumps(v, ensure_ascii=False)
    return v


def _auto_representative_indices(nspec: int, count: int) -> List[int]:
    if nspec <= 0 or count <= 0:
        return [0]
    step = int(math.ceil(nspec / float(count)))
    idxs = []
    cur = nspec - 1
    while len(idxs) < count and cur >= 0:
        idxs.append(cur)
        cur -= step
    if not idxs:
        idxs = [nspec - 1]
    return idxs


def _select_indices_by_pressure(processed: Sequence[Spectrum], target_pressures: Sequence[float]) -> List[int]:
    if not processed:
        return []
    available = np.asarray([float(s.pressure) for s in processed], dtype=float)
    chosen: List[int] = []
    used: set[int] = set()
    for pt in target_pressures:
        diffs = np.abs(available - float(pt))
        order = np.argsort(diffs)
        pick = None
        for idx in order:
            ii = int(idx)
            if ii not in used:
                pick = ii
                break
        if pick is None and len(order) > 0:
            pick = int(order[0])
        if pick is not None:
            chosen.append(pick)
            used.add(pick)
    return sorted(chosen, key=lambda i: float(processed[i].pressure))


def plot_selected_normalized_main_band(
    spectra: Sequence[Spectrum],
    outpath: Path,
    panel_cfg: Dict[str, Any],
    cfg: Dict[str, Any],
    spec_order: Dict[str, int],
    dpi: int,
) -> pd.DataFrame:
    if not spectra:
        return pd.DataFrame()

    window_raw = panel_cfg.get('window', cfg.get('analysis_window'))
    if not isinstance(window_raw, (list, tuple)) or len(window_raw) != 2:
        raise ValueError('selected_main_band.window must contain exactly two values.')
    xwin = (float(window_raw[0]), float(window_raw[1]))

    norm_raw = panel_cfg.get('normalize_window', window_raw)
    if not isinstance(norm_raw, (list, tuple)) or len(norm_raw) != 2:
        raise ValueError('selected_main_band.normalize_window must contain exactly two values.')
    nwin = (float(norm_raw[0]), float(norm_raw[1]))

    selected_indices = panel_cfg.get('selected_indices')
    selected_pressures = panel_cfg.get('selected_pressures')
    if selected_indices:
        idxs = sorted({int(v) for v in selected_indices if 0 <= int(v) < len(spectra)}, key=lambda i: float(spectra[i].pressure))
    elif selected_pressures:
        idxs = _select_indices_by_pressure(spectra, [float(v) for v in selected_pressures])
    else:
        max_count = int(panel_cfg.get('max_count', 5))
        idxs = _auto_representative_indices(len(spectra), max_count)
        idxs = sorted({int(v) for v in idxs if 0 <= int(v) < len(spectra)}, key=lambda i: float(spectra[i].pressure))

    if not idxs:
        return pd.DataFrame()

    mode = str(panel_cfg.get('style', 'stacked')).strip().lower()
    offset_step = float(panel_cfg.get('offset_step', 0.85))
    lw = float(panel_cfg.get('linewidth', 1.2))
    color = str(panel_cfg.get('line_color', 'black'))
    figsize = tuple(panel_cfg.get('figsize', (6.4, 4.6)))
    base_label_mode = str(panel_cfg.get('label_mode', 'legend')).strip().lower()
    right_pad_frac = float(panel_cfg.get('right_label_pad_frac', 0.02))
    title = str(panel_cfg.get('title', 'Selected normalized high-resolution spectra in the main-band region'))
    clean_title = str(panel_cfg.get('clean_title', title))
    hide_y_ticks = bool(panel_cfg.get('hide_y_ticks', True))
    ylabel = 'Normalized intensity' if mode == 'overlay' else 'Normalized intensity + offset'

    rows: List[Dict[str, Any]] = []
    xmin, xmax = xwin
    xspan = max(xmax - xmin, 1.0)
    traces: List[Dict[str, Any]] = []

    for j, idx in enumerate(idxs):
        spec = spectra[idx]
        target_masks = get_target_mask_ranges(cfg, 'selected_main_band', spec_order.get(spec.label), spec.pressure, spec.label)
        yplot = apply_nan_masks(spec.x, spec.y_corrected, target_masks)
        mask = in_window(spec.x, xwin) & np.isfinite(yplot)
        if np.count_nonzero(mask) < 3:
            continue
        xsel = np.asarray(spec.x[mask], dtype=float)
        ysel = np.asarray(yplot[mask], dtype=float)

        nmask = in_window(xsel, nwin) & np.isfinite(ysel)
        if np.count_nonzero(nmask) < 3:
            nmask = np.isfinite(ysel)
        if np.count_nonzero(nmask) < 3:
            continue

        yref = ysel[nmask]
        ymin = float(np.nanmin(yref))
        ymax = float(np.nanmax(yref))
        if (not np.isfinite(ymin)) or (not np.isfinite(ymax)) or ymax <= ymin:
            ynorm = np.zeros_like(ysel)
        else:
            ynorm = (ysel - ymin) / (ymax - ymin)

        ydraw = ynorm if mode == 'overlay' else ynorm + j * offset_step
        label = str(panel_cfg.get('label_template', '{pressure:.2f} GPa')).format(label=spec.label, pressure=spec.pressure)
        traces.append({'x': xsel, 'y': ydraw, 'label': label, 'pressure': float(spec.pressure), 'index': j})

        rows.append({
            'label': spec.label,
            'pressure': float(spec.pressure),
            'x_min': xmin,
            'x_max': xmax,
            'norm_x_min': nwin[0],
            'norm_x_max': nwin[1],
        })

    if not traces:
        return pd.DataFrame(rows)

    def _render_selected_main_band(fig_title: str, label_mode: str, fig_out: Path) -> None:
        fig, ax = plt.subplots(figsize=figsize)
        for item in traces:
            ax.plot(item['x'], item['y'], lw=lw, color=color, label=item['label'])
            if label_mode == 'right' and mode != 'overlay':
                xtext = xmax - right_pad_frac * xspan
                finite = np.isfinite(item['y'])
                if np.any(finite):
                    ytext = float(item['y'][finite][-1])
                    ax.text(xtext, ytext, item['label'], ha='right', va='center', fontsize=float(panel_cfg.get('label_fontsize', 8.0)))
        ax.set_xlim(xwin)
        ax.set_xlabel('Raman shift (cm$^{-1}$)')
        ax.set_ylabel(ylabel)
        ax.set_title(fig_title)
        if hide_y_ticks:
            ax.set_yticks([])
        _clean_axes(ax)
        if label_mode == 'legend':
            ax.legend(frameon=False, loc=str(panel_cfg.get('legend_loc', 'upper right')), fontsize=float(panel_cfg.get('legend_fontsize', 8.0)))
        fig.tight_layout()
        fig.savefig(fig_out, dpi=dpi, bbox_inches='tight')
        plt.close(fig)

    _render_selected_main_band(title, base_label_mode, outpath)

    if bool(panel_cfg.get('output_unlabeled_copy', False)):
        clean_name = str(panel_cfg.get('clean_output_name', outpath.stem + '_clean.png'))
        clean_path = outpath.with_name(clean_name)
        _render_selected_main_band(clean_title, 'none', clean_path)

    return pd.DataFrame(rows)


def process_profile(profile: Dict[str, Any]) -> Dict[str, Any]:
    cfg = deep_copy_config(profile)
    progress_enabled = bool(cfg.get('progress', True))
    outdir = resolve_output_dir_from_cfg(cfg)
    cfg['output_dir'] = str(outdir)
    ensure_dir(outdir)

    excel_path = Path(cfg['excel_path'])
    if not excel_path.exists():
        raise FileNotFoundError(f'Excel file not found: {excel_path}')

    stage_bar = make_progress_bar(total=9, desc=f"{cfg.get('name', 'profile')} stages", enabled=progress_enabled, leave=True)
    try:
        stage_bar.set_postfix_str('read excel')
        raw = read_excel_spectra(
            excel_path=str(excel_path),
            sheet_name=cfg.get('sheet_name', 0),
            read_mode=cfg.get('read_mode', 'auto'),
            labels=cfg.get('spectrum_labels'),
        )
        nspec = len(raw)
        pressures = np.asarray(cfg.get('pressure_list', []), dtype=float)
        if len(pressures) != nspec:
            raise ValueError(f'Pressure count mismatch: found {nspec} spectra, but got {len(pressures)} pressures.')
        stage_bar.update(1)

        stage_bar.set_postfix_str('baseline / correction')
        rate_factors = build_rate_factors(cfg.get('rate_normalization', {}), nspec)
        als_cfg = cfg.get('als', {})
        zero_ref_cfg = cfg.get('zero_reference', {})

        processed: List[Spectrum] = []
        spec_order: Dict[str, int] = {}
        raw_iter = progress_iter(raw, desc=f"{cfg.get('name', 'profile')} preprocess", total=nspec, enabled=progress_enabled, leave=False)
        for i, (label, x, yraw) in enumerate(raw_iter):
            y_rate = yraw.copy() if rate_factors is None else yraw / rate_factors[i]
            baseline_masks = get_target_mask_ranges(cfg, 'baseline', i, float(pressures[i]), label)
            y_for_baseline = fill_masked_regions_linear(x, y_rate, baseline_masks)
            baseline = als_baseline(
                y_for_baseline,
                lam=float(als_cfg.get('lambda', 1e7)),
                p=float(als_cfg.get('p', 0.01)),
                niter=int(als_cfg.get('niter', 30)),
            )
            ycorr = y_rate - baseline
            item_zero_cfg = copy.deepcopy(zero_ref_cfg) if isinstance(zero_ref_cfg, dict) else {}
            extra_excludes = get_target_mask_ranges(cfg, 'zero_reference', i, float(pressures[i]), label)
            item_zero_cfg['exclude_ranges'] = _normalize_mask_ranges(list(item_zero_cfg.get('exclude_ranges', [])) + extra_excludes)
            ycorr, zero_offset = maybe_apply_zero_reference(x, ycorr, item_zero_cfg)
            spec = Spectrum(
                label=label,
                pressure=float(pressures[i]),
                x=np.asarray(x, dtype=float),
                y_raw=np.asarray(yraw, dtype=float),
                y_rate_norm=np.asarray(y_rate, dtype=float),
                baseline=np.asarray(baseline, dtype=float),
                y_corrected=np.asarray(ycorr, dtype=float),
                zero_reference_offset=float(zero_offset),
            )
            processed.append(spec)
            spec_order[label] = i
        stage_bar.update(1)

        plot_cfg = cfg.get('plots', {})
        common_plot_window = tuple(cfg['analysis_window']) if cfg.get('analysis_window') else None
        sort_by_pressure = bool(cfg.get('sort_by_pressure_for_plots', True))
        dpi = int(cfg.get('dpi', 300))

        stage_bar.set_postfix_str('waterfall')
        if plot_cfg.get('waterfall', True):
            waterfall_cfg = cfg.get('waterfall', {})
            waterfall_mask_map = {
                spec.label: get_target_mask_ranges(cfg, 'waterfall', spec_order[spec.label], spec.pressure, spec.label)
                for spec in processed
            }
            plot_waterfall(
                spectra=processed,
                outpath=outdir / '01_waterfall_bg_subtracted.png',
                window=common_plot_window,
                mask_ranges=waterfall_mask_map,
                vertical_offset=float(waterfall_cfg.get('vertical_offset', 1.00)),
                linewidth=float(waterfall_cfg.get('linewidth', 1.0)),
                sort_by_pressure=sort_by_pressure,
                dpi=dpi,
                line_color=str(waterfall_cfg.get('line_color', 'black')),
                normalize_each=bool(waterfall_cfg.get('normalize_each', False)),
                figsize=tuple(waterfall_cfg.get('figsize', (8.6, 12.0))),
                label_fontsize=float(waterfall_cfg.get('label_fontsize', 8.0)),
                right_label_offset_frac=float(waterfall_cfg.get('right_label_offset_frac', 0.025)),
                right_label_inside=bool(waterfall_cfg.get('right_label_inside', True)),
                right_label_pad_frac=float(waterfall_cfg.get('right_label_pad_frac', 0.012)),
                pressure_label_mode=str(waterfall_cfg.get('pressure_label_mode', 'all')),
                pressure_label_every=int(waterfall_cfg.get('pressure_label_every', 3)),
                pressure_label_include_first_last=bool(waterfall_cfg.get('pressure_label_include_first_last', True)),
                hide_y_ticks=bool(waterfall_cfg.get('hide_y_ticks', True)),
                title=str(waterfall_cfg.get('title', 'Waterfall diagram (ALS BG-subtracted)')),
            )
            if bool(waterfall_cfg.get('output_unlabeled_copy', True)):
                plot_waterfall(
                    spectra=processed,
                    outpath=outdir / '01_waterfall_bg_subtracted_clean.png',
                    window=common_plot_window,
                    mask_ranges=waterfall_mask_map,
                    vertical_offset=float(waterfall_cfg.get('vertical_offset', 1.00)),
                    linewidth=float(waterfall_cfg.get('linewidth', 1.0)),
                    sort_by_pressure=sort_by_pressure,
                    dpi=dpi,
                    line_color=str(waterfall_cfg.get('line_color', 'black')),
                    normalize_each=bool(waterfall_cfg.get('normalize_each', False)),
                    figsize=tuple(waterfall_cfg.get('figsize', (8.6, 12.0))),
                    label_fontsize=float(waterfall_cfg.get('label_fontsize', 8.0)),
                    right_label_offset_frac=float(waterfall_cfg.get('right_label_offset_frac', 0.025)),
                    right_label_inside=bool(waterfall_cfg.get('right_label_inside', True)),
                    right_label_pad_frac=float(waterfall_cfg.get('right_label_pad_frac', 0.012)),
                    pressure_label_mode='none',
                    pressure_label_every=int(waterfall_cfg.get('pressure_label_every', 3)),
                    pressure_label_include_first_last=False,
                    hide_y_ticks=bool(waterfall_cfg.get('hide_y_ticks', True)),
                    title=str(waterfall_cfg.get('clean_title', 'Waterfall diagram (clean)')),
                )
        stage_bar.update(1)

        heatmap_items = sort_spectra(processed, sort_by_pressure)
        xgrid = None
        pgrid = None
        zrows: List[np.ndarray] = []
        xyz_rows: List[Dict[str, float]] = []
        heatmap_cfg = cfg.get('heatmap', {})
        stage_bar.set_postfix_str('heatmap')
        if plot_cfg.get('heatmap', True):
            analysis_window = tuple(heatmap_cfg.get('window', cfg['analysis_window']))
            xgrid = make_common_grid(analysis_window, float(heatmap_cfg.get('grid_step', 1.0)))
            pgrid = np.array([s.pressure for s in heatmap_items], dtype=float)
            norm_mode = str(heatmap_cfg.get('normalize_mode', 'per_spectrum_minmax'))
            norm_window_rules = heatmap_cfg.get('normalization_windows_by_pressure', [])
            fill_rules = heatmap_cfg.get('conditional_fill_regions', [])
            heatmap_iter = progress_iter(heatmap_items, desc=f"{cfg.get('name', 'profile')} heatmap", total=len(heatmap_items), enabled=progress_enabled, leave=False)
            for spec in heatmap_iter:
                spec_idx = spec_order[spec.label]
                heatmap_masks = get_target_mask_ranges(cfg, 'heatmap', spec_idx, spec.pressure, spec.label)
                y_interp = interpolate_spectrum(spec.x, spec.y_corrected, xgrid)
                y_interp = apply_nan_masks(xgrid, y_interp, heatmap_masks)

                norm_window = _select_heatmap_norm_window(spec.pressure, analysis_window, norm_window_rules)
                if norm_mode.lower() in {'minmax', 'per_spectrum_minmax', 'max', 'peakmax', 'per_spectrum_max'}:
                    norm_mask = in_window(xgrid, norm_window) & np.isfinite(y_interp)
                    norm_mask = _exclude_heatmap_rules_from_norm(norm_mask, xgrid, spec.pressure, fill_rules)
                    y_norm = np.full_like(y_interp, np.nan, dtype=float)
                    if np.any(norm_mask):
                        ref = y_interp[norm_mask]
                        if norm_mode.lower() in {'max', 'peakmax', 'per_spectrum_max'}:
                            ymax = np.nanmax(ref)
                            y_norm[np.isfinite(y_interp)] = 0.0 if (not np.isfinite(ymax) or ymax == 0) else y_interp[np.isfinite(y_interp)] / ymax
                        else:
                            ymin = np.nanmin(ref)
                            ymax = np.nanmax(ref)
                            if (not np.isfinite(ymin)) or (not np.isfinite(ymax)) or ymax <= ymin:
                                y_norm[np.isfinite(y_interp)] = 0.0
                            else:
                                y_norm[np.isfinite(y_interp)] = (y_interp[np.isfinite(y_interp)] - ymin) / (ymax - ymin)
                        y_norm[np.isfinite(y_norm)] = np.clip(y_norm[np.isfinite(y_norm)], 0.0, 1.0)
                    else:
                        y_norm[np.isfinite(y_interp)] = 0.0
                else:
                    y_norm = normalize_vector(y_interp, norm_mode)

                y_norm = _apply_heatmap_fill_rules(y_norm, xgrid, spec.pressure, fill_rules)
                zrows.append(y_norm)
                for xv, yv in zip(xgrid, y_norm):
                    xyz_rows.append({'pressure': float(spec.pressure), 'x': float(xv), 'intensity': float(yv) if np.isfinite(yv) else np.nan})
            zmat = np.vstack(zrows)
            plot_heatmap(xgrid, pgrid, zmat, outdir / '02_heatmap.png', dpi=dpi)
        else:
            zmat = None
        stage_bar.update(1)

        df_peak = pd.DataFrame()
        peak_cfg = cfg.get('peak_tracking', {})
        peak_results: List[Dict[str, Any]] = []
        peak_items = sort_spectra(processed, sort_by_pressure)
        stage_bar.set_postfix_str('peak tracking')
        if plot_cfg.get('pplot', True) and peak_cfg.get('enabled', True):
            peak_iter = progress_iter(peak_items, desc=f"{cfg.get('name', 'profile')} peak fit", total=len(peak_items), enabled=progress_enabled, leave=False)
            for spec in peak_iter:
                spec_idx = spec_order[spec.label]
                local_peak_cfg = resolve_peak_fit_cfg(peak_cfg, spec_idx=spec_idx, pressure=spec.pressure, label=spec.label)
                local_peak_cfg['exclude_ranges'] = _normalize_mask_ranges(list(local_peak_cfg.get('exclude_ranges', [])) + get_target_mask_ranges(cfg, 'fit', spec_idx, spec.pressure, spec.label))
                local_peak_cfg['display_mask_ranges'] = _normalize_mask_ranges(list(local_peak_cfg.get('display_mask_ranges', [])) + get_target_mask_ranges(cfg, 'display', spec_idx, spec.pressure, spec.label))
                res = run_peak_fit(spec, local_peak_cfg)
                res['display_mask_ranges'] = local_peak_cfg.get('display_mask_ranges', [])
                peak_results.append(res)
            df_peak = pd.DataFrame([{k: v for k, v in r.items() if k not in {'fit_x', 'fit_y', 'fit_yhat', 'fit_y_bg', 'peak1_y', 'peak2_y', 'display_mask_ranges'}} for r in peak_results])
            plot_pplot_center(df_peak, outdir / '03_pplot_peak_center.png', dpi=dpi, plot_cfg=cfg.get('pplot_center', {}))
            plot_pplot_width(df_peak, outdir / '04_pplot_fwhm.png', dpi=dpi, plot_cfg=cfg.get('pplot_fwhm', {}))
        stage_bar.update(1)

        spectral_cfg = cfg.get('spectral_weight', {})
        df_spectral_intensity = pd.DataFrame()
        df_spectral_fraction = pd.DataFrame()
        df_spectral_ratio = pd.DataFrame()
        stage_bar.set_postfix_str('spectral weight')
        if bool(spectral_cfg.get('enabled', False)):
            df_spectral_intensity, df_spectral_fraction, df_spectral_ratio = compute_spectral_weight_tables(
                spectra=processed,
                spectral_cfg=spectral_cfg,
                cfg=cfg,
                spec_order=spec_order,
                sort_by_pressure=sort_by_pressure,
            )
            if plot_cfg.get('spectral_weight', True):
                if not df_spectral_intensity.empty:
                    plot_integrated_intensities(df_spectral_intensity, outdir / '06_integrated_intensity.png', dpi=dpi)
                if not df_spectral_fraction.empty:
                    plot_spectral_fractions(df_spectral_fraction, outdir / '07_spectral_fraction.png', dpi=dpi)
                if not df_spectral_ratio.empty and len([c for c in df_spectral_ratio.columns if c not in {'label', 'pressure'}]) > 0:
                    plot_spectral_ratios(df_spectral_ratio, outdir / '08_spectral_ratio.png', dpi=dpi)
        stage_bar.update(1)

        df_selected_main_band = pd.DataFrame()
        selected_main_band_cfg = cfg.get('selected_main_band', {})
        stage_bar.set_postfix_str('selected main band')
        if bool(selected_main_band_cfg.get('enabled', False)):
            df_selected_main_band = plot_selected_normalized_main_band(
                spectra=sort_spectra(processed, sort_by_pressure),
                outpath=outdir / str(selected_main_band_cfg.get('output_name', '09_selected_main_band.png')),
                panel_cfg=selected_main_band_cfg,
                cfg=cfg,
                spec_order=spec_order,
                dpi=dpi,
            )
        stage_bar.update(1)

        rep_indices_cfg = cfg.get('representative_spectrum_indices')
        if rep_indices_cfg is None and int(cfg.get('representative_count_from_end', 0)) > 0:
            rep_indices = _auto_representative_indices(len(processed), int(cfg.get('representative_count_from_end', 4)))
        elif rep_indices_cfg is None:
            rep_indices = [int(cfg.get('representative_spectrum_index', 0))]
        else:
            rep_indices = [int(v) for v in rep_indices_cfg]
        if not rep_indices:
            raise ValueError('representative_spectrum_indices must contain at least one index.')
        rep_specs = []
        rep_peaks = []
        for rep_idx in rep_indices:
            if not (0 <= rep_idx < len(processed)):
                raise IndexError('representative spectrum index is out of range.')
            rep_spec = processed[rep_idx]
            rep_specs.append(rep_spec)
            rep_peak = None
            if peak_results:
                for r in peak_results:
                    if r.get('label') == rep_spec.label:
                        rep_peak = r
                        break
            rep_peaks.append(rep_peak)
        stage_bar.set_postfix_str('representative plot')
        if plot_cfg.get('representative', True):
            plot_representative(rep_specs, rep_peaks, outdir / '05_representative_baseline_and_fit.png', common_plot_window, peak_cfg, dpi)
        stage_bar.update(1)

        stage_bar.set_postfix_str('excel export')
        if cfg.get('save_excel_summary', True):
            summary_path = outdir / cfg.get('excel_summary_name', 'raman_processed_summary.xlsx')
            meta_rows = [{'key': k, 'value': serializable_value(v)} for k, v in cfg.items()]
            meta_rows.append({'key': 'n_spectra', 'value': nspec})
            meta_rows.append({'key': 'read_excel_path', 'value': str(excel_path.resolve())})

            df_meta = pd.DataFrame(meta_rows)
            df_raw = spectra_to_wide_xy(processed, 'y_raw')
            df_rate = spectra_to_wide_xy(processed, 'y_rate_norm')
            df_base = spectra_to_wide_xy(processed, 'baseline')
            df_corr = spectra_to_wide_xy(processed, 'y_corrected')
            df_xyz = pd.DataFrame(xyz_rows)
            df_zero = pd.DataFrame([{
                'label': spec.label,
                'pressure': spec.pressure,
                'zero_reference_offset': spec.zero_reference_offset,
            } for spec in processed])

            with pd.ExcelWriter(summary_path, engine='openpyxl') as writer:
                df_meta.to_excel(writer, sheet_name='Metadata', index=False)
                df_raw.to_excel(writer, sheet_name='Raw', index=False)
                df_rate.to_excel(writer, sheet_name='RateNorm', index=False)
                df_base.to_excel(writer, sheet_name='Baseline', index=False)
                df_corr.to_excel(writer, sheet_name='Corrected', index=False)
                if not df_zero.empty:
                    df_zero.to_excel(writer, sheet_name='ZeroReference', index=False)
                if xgrid is not None and zrows:
                    norm_data = {'x': xgrid}
                    for spec, row in zip(heatmap_items, zrows):
                        norm_data[f'{spec.label}_P_{spec.pressure:.6f}'] = row
                    pd.DataFrame(norm_data).to_excel(writer, sheet_name='Normalized', index=False)
                    df_xyz.to_excel(writer, sheet_name='Heatmap_XYZ', index=False)
                if not df_peak.empty:
                    df_peak.to_excel(writer, sheet_name='PeakTracking', index=False)
                if not df_spectral_intensity.empty:
                    df_spectral_intensity.to_excel(writer, sheet_name='IntegratedIntensity', index=False)
                if not df_spectral_fraction.empty:
                    df_spectral_fraction.to_excel(writer, sheet_name='SpectralFractions', index=False)
                if not df_spectral_ratio.empty and len(df_spectral_ratio.columns) > 2:
                    df_spectral_ratio.to_excel(writer, sheet_name='SpectralRatios', index=False)
                if 'df_selected_main_band' in locals() and not df_selected_main_band.empty:
                    df_selected_main_band.to_excel(writer, sheet_name='SelectedMainBand', index=False)
        stage_bar.update(1)
    finally:
        stage_bar.close()

    return {
        'config': cfg,
        'n_spectra': nspec,
        'output_dir': str(outdir.resolve()),
        'peak_table': df_peak,
        'spectral_intensity_table': df_spectral_intensity if 'df_spectral_intensity' in locals() else pd.DataFrame(),
        'spectral_fraction_table': df_spectral_fraction if 'df_spectral_fraction' in locals() else pd.DataFrame(),
        'spectral_ratio_table': df_spectral_ratio if 'df_spectral_ratio' in locals() else pd.DataFrame(),
        'selected_main_band_table': df_selected_main_band if 'df_selected_main_band' in locals() else pd.DataFrame(),
    }



def _safe_numeric_df_for_compare(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    numeric_cols = [
        'pressure', 'center', 'center_se', 'center1', 'center1_se', 'center2', 'center2_se',
        'effective_center', 'effective_center_se', 'midpoint_center', 'fwhm', 'fwhm1', 'fwhm2',
        'sigma', 'sigma1', 'sigma2', 'amplitude', 'amplitude1', 'amplitude2', 'area', 'area1', 'area2',
        'bg_slope', 'bg_intercept', 'peak_separation', 'rss', 'aic', 'bic', 'n_points'
    ]
    for col in numeric_cols:
        if col in out.columns:
            try:
                out[col] = pd.to_numeric(out[col])
            except (ValueError, TypeError):
                pass
    return out


def _prepare_peak_table_for_compare(df: pd.DataFrame) -> pd.DataFrame:
    out = _safe_numeric_df_for_compare(df)
    if 'label' not in out.columns:
        out['label'] = [f"S{i+1:02d}" for i in range(len(out))]
    if 'mode_used' not in out.columns:
        out['mode_used'] = ''
    if 'success' not in out.columns:
        out['success'] = True
    if 'pressure' not in out.columns:
        raise ValueError('Peak table must contain a pressure column for comparison.')
    out = out.sort_values(['pressure', 'label']).reset_index(drop=True)
    return out


def _metric_se_column(metric: str) -> str:
    if metric.endswith('_se'):
        return metric
    return f'{metric}_se'


def _get_metric_frame(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    se_col = _metric_se_column(metric)
    cols = ['label', 'pressure', metric]
    if se_col in df.columns:
        cols.append(se_col)
    if 'mode_used' in df.columns:
        cols.append('mode_used')
    out = df.loc[:, [c for c in cols if c in df.columns]].copy()
    if metric not in out.columns:
        return pd.DataFrame(columns=['label', 'pressure', 'value', 'value_se', 'mode_used'])
    out = out.rename(columns={metric: 'value'})
    if se_col in out.columns:
        out = out.rename(columns={se_col: 'value_se'})
    else:
        out['value_se'] = np.nan
    if 'mode_used' not in out.columns:
        out['mode_used'] = ''
    out = out[np.isfinite(out['pressure']) & np.isfinite(out['value'])].copy()
    out = out.sort_values(['pressure', 'label']).reset_index(drop=True)
    return out[['label', 'pressure', 'value', 'value_se', 'mode_used']]


def _interp_linear_safe(x_src: np.ndarray, y_src: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    x_src = np.asarray(x_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    x_new = np.asarray(x_new, dtype=float)
    mask = np.isfinite(x_src) & np.isfinite(y_src)
    if np.count_nonzero(mask) == 0:
        return np.full_like(x_new, np.nan, dtype=float)
    xs = x_src[mask]
    ys = y_src[mask]
    order = np.argsort(xs)
    xs = xs[order]
    ys = ys[order]
    xs, idx = np.unique(xs, return_index=True)
    ys = ys[idx]
    if len(xs) == 1:
        out = np.full_like(x_new, np.nan, dtype=float)
        only = np.isclose(x_new, xs[0])
        out[only] = ys[0]
        return out
    out = np.interp(x_new, xs, ys)
    outside = (x_new < xs[0]) | (x_new > xs[-1])
    out[outside] = np.nan
    return out


def _weighted_linear_fit_general(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> Optional[Dict[str, float]]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0)
    if np.count_nonzero(mask) < 2:
        return None
    xv = x[mask]
    yv = y[mask]
    wv = weights[mask]
    X = np.column_stack([xv, np.ones_like(xv)])
    sw = np.sqrt(wv)
    Xw = X * sw[:, None]
    yw = yv * sw
    try:
        beta, _, _, _ = np.linalg.lstsq(Xw, yw, rcond=None)
    except np.linalg.LinAlgError:
        return None
    yhat = X @ beta
    resid = yv - yhat
    dof = max(len(yv) - 2, 1)
    s2 = float(np.sum(wv * resid * resid) / dof)
    try:
        cov = np.linalg.inv(Xw.T @ Xw) * s2
    except np.linalg.LinAlgError:
        cov = np.full((2, 2), np.nan)
    slope = float(beta[0])
    intercept = float(beta[1])
    slope_se = float(np.sqrt(cov[0, 0])) if np.isfinite(cov[0, 0]) and cov[0, 0] >= 0 else np.nan
    intercept_se = float(np.sqrt(cov[1, 1])) if np.isfinite(cov[1, 1]) and cov[1, 1] >= 0 else np.nan
    ss_res = float(np.sum((yv - yhat) ** 2))
    ss_tot = float(np.sum((yv - np.mean(yv)) ** 2))
    r2 = np.nan if ss_tot <= 0 else 1.0 - ss_res / ss_tot
    return {
        'slope': slope,
        'intercept': intercept,
        'slope_se': slope_se,
        'intercept_se': intercept_se,
        'n': int(len(yv)),
        'r2': float(r2) if np.isfinite(r2) else np.nan,
    }


def robust_linear_fit(
    x: np.ndarray,
    y: np.ndarray,
    yerr: Optional[np.ndarray] = None,
    max_iter: int = 50,
    huber_c: float = 1.345,
) -> Optional[Dict[str, float]]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if yerr is None:
        base_weights = np.ones_like(x, dtype=float)
    else:
        yerr = np.asarray(yerr, dtype=float)
        base_weights = np.ones_like(x, dtype=float)
        mask_err = np.isfinite(yerr) & (yerr > 0)
        base_weights[mask_err] = 1.0 / (yerr[mask_err] ** 2)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(base_weights) & (base_weights > 0)
    if np.count_nonzero(mask) < 2:
        return None
    xv = x[mask]
    yv = y[mask]
    bw = base_weights[mask]
    fit = _weighted_linear_fit_general(xv, yv, bw)
    if fit is None:
        return None
    slope = fit['slope']
    intercept = fit['intercept']
    weights = bw.copy()
    for _ in range(max_iter):
        resid = yv - (slope * xv + intercept)
        med = np.nanmedian(resid)
        mad = np.nanmedian(np.abs(resid - med))
        scale = 1.4826 * mad
        if not np.isfinite(scale) or scale <= 1e-12:
            break
        u = resid / (huber_c * scale)
        huber_w = np.ones_like(u)
        big = np.abs(u) > 1.0
        huber_w[big] = 1.0 / np.abs(u[big])
        new_weights = bw * huber_w
        new_fit = _weighted_linear_fit_general(xv, yv, new_weights)
        if new_fit is None:
            break
        ds = abs(new_fit['slope'] - slope)
        di = abs(new_fit['intercept'] - intercept)
        slope = new_fit['slope']
        intercept = new_fit['intercept']
        fit = new_fit
        weights = new_weights
        if ds < 1e-8 and di < 1e-6:
            break
    fit = dict(fit)
    fit['method'] = 'robust_huber'
    fit['weight_sum'] = float(np.sum(weights))
    return fit


def _pair_union_grid(
    left: pd.DataFrame,
    right: pd.DataFrame,
) -> Tuple[np.ndarray, float, float]:
    lmin = float(np.nanmin(left['pressure'].to_numpy(dtype=float)))
    lmax = float(np.nanmax(left['pressure'].to_numpy(dtype=float)))
    rmin = float(np.nanmin(right['pressure'].to_numpy(dtype=float)))
    rmax = float(np.nanmax(right['pressure'].to_numpy(dtype=float)))
    lo = max(lmin, rmin)
    hi = min(lmax, rmax)
    if hi < lo:
        return np.array([], dtype=float), lo, hi
    xl = left['pressure'].to_numpy(dtype=float)
    xr = right['pressure'].to_numpy(dtype=float)
    grid = np.unique(np.concatenate([xl[(xl >= lo) & (xl <= hi)], xr[(xr >= lo) & (xr <= hi)]]))
    grid = grid[np.isfinite(grid)]
    return np.sort(grid), lo, hi


def _pair_series_on_grid(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    xgrid, lo, hi = _pair_union_grid(left, right)
    if len(xgrid) == 0:
        return pd.DataFrame(columns=[
            'pressure', 'left_value', 'left_value_se', 'right_value', 'right_value_se',
            'delta_left_minus_right', 'left_observed', 'right_observed'
        ])
    xl = left['pressure'].to_numpy(dtype=float)
    yl = left['value'].to_numpy(dtype=float)
    yls = left['value_se'].to_numpy(dtype=float) if 'value_se' in left.columns else np.full(len(left), np.nan)
    xr = right['pressure'].to_numpy(dtype=float)
    yr = right['value'].to_numpy(dtype=float)
    yrs = right['value_se'].to_numpy(dtype=float) if 'value_se' in right.columns else np.full(len(right), np.nan)
    left_interp = _interp_linear_safe(xl, yl, xgrid)
    left_se_interp = _interp_linear_safe(xl, yls, xgrid)
    right_interp = _interp_linear_safe(xr, yr, xgrid)
    right_se_interp = _interp_linear_safe(xr, yrs, xgrid)
    out = pd.DataFrame({
        'pressure': xgrid,
        'left_value': left_interp,
        'left_value_se': left_se_interp,
        'right_value': right_interp,
        'right_value_se': right_se_interp,
    })
    out['left_observed'] = np.isin(np.round(out['pressure'].to_numpy(dtype=float), 10), np.round(xl, 10))
    out['right_observed'] = np.isin(np.round(out['pressure'].to_numpy(dtype=float), 10), np.round(xr, 10))
    out['delta_left_minus_right'] = out['left_value'] - out['right_value']
    out['delta_se'] = np.sqrt(np.square(out['left_value_se']) + np.square(out['right_value_se']))
    out = out[np.isfinite(out['left_value']) & np.isfinite(out['right_value'])].reset_index(drop=True)
    return out


def _pair_correlation(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    if np.count_nonzero(mask) < 2:
        return np.nan
    av = a[mask]
    bv = b[mask]
    if np.nanstd(av) <= 0 or np.nanstd(bv) <= 0:
        return np.nan
    return float(np.corrcoef(av, bv)[0, 1])


def _trend_sign_changes(y: np.ndarray) -> int:
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if len(y) < 4:
        return 0
    dy = np.diff(y)
    dy = dy[np.isfinite(dy)]
    if len(dy) < 3:
        return 0
    thr = max(np.nanmedian(np.abs(dy)) * 0.25, 1e-9)
    signs = []
    for v in dy:
        if abs(v) < thr:
            continue
        signs.append(np.sign(v))
    if len(signs) < 3:
        return 0
    changes = 0
    prev = signs[0]
    for s in signs[1:]:
        if s != prev:
            changes += 1
        prev = s
    return int(changes)


def _leave_one_out_slope_deltas(
    x: np.ndarray,
    yl: np.ndarray,
    yls: np.ndarray,
    yr: np.ndarray,
    yrs: np.ndarray,
) -> Dict[str, Any]:
    x = np.asarray(x, dtype=float)
    yl = np.asarray(yl, dtype=float)
    yr = np.asarray(yr, dtype=float)
    yls = np.asarray(yls, dtype=float)
    yrs = np.asarray(yrs, dtype=float)
    mask = np.isfinite(x) & np.isfinite(yl) & np.isfinite(yr)
    x = x[mask]
    yl = yl[mask]
    yr = yr[mask]
    yls = yls[mask] if len(yls) == len(mask) else np.full(len(x), np.nan)
    yrs = yrs[mask] if len(yrs) == len(mask) else np.full(len(x), np.nan)
    if len(x) < 4:
        return {
            'n_loo': 0,
            'loo_slope_delta_min': np.nan,
            'loo_slope_delta_max': np.nan,
            'loo_slope_delta_median': np.nan,
            'loo_sign_agreement': np.nan,
            'loo_sig_agreement': np.nan,
        }
    full_l = robust_linear_fit(x, yl, yls)
    full_r = robust_linear_fit(x, yr, yrs)
    if full_l is None or full_r is None:
        return {
            'n_loo': 0,
            'loo_slope_delta_min': np.nan,
            'loo_slope_delta_max': np.nan,
            'loo_slope_delta_median': np.nan,
            'loo_sign_agreement': np.nan,
            'loo_sig_agreement': np.nan,
        }
    full_delta = full_l['slope'] - full_r['slope']
    full_scale = math.sqrt(
        (full_l.get('slope_se', np.nan) ** 2 if np.isfinite(full_l.get('slope_se', np.nan)) else 0.0) +
        (full_r.get('slope_se', np.nan) ** 2 if np.isfinite(full_r.get('slope_se', np.nan)) else 0.0)
    )
    full_sig = abs(full_delta / full_scale) >= 2.0 if np.isfinite(full_scale) and full_scale > 0 else False
    vals = []
    sign_match = 0
    sig_match = 0
    n_ok = 0
    full_sign = 0.0 if abs(full_delta) < 1e-12 else float(np.sign(full_delta))
    for i in range(len(x)):
        keep = np.ones(len(x), dtype=bool)
        keep[i] = False
        fit_l = robust_linear_fit(x[keep], yl[keep], yls[keep])
        fit_r = robust_linear_fit(x[keep], yr[keep], yrs[keep])
        if fit_l is None or fit_r is None:
            continue
        delta = fit_l['slope'] - fit_r['slope']
        vals.append(delta)
        n_ok += 1
        sign = 0.0 if abs(delta) < 1e-12 else float(np.sign(delta))
        if sign == full_sign:
            sign_match += 1
        scale = math.sqrt(
            (fit_l.get('slope_se', np.nan) ** 2 if np.isfinite(fit_l.get('slope_se', np.nan)) else 0.0) +
            (fit_r.get('slope_se', np.nan) ** 2 if np.isfinite(fit_r.get('slope_se', np.nan)) else 0.0)
        )
        sig = abs(delta / scale) >= 2.0 if np.isfinite(scale) and scale > 0 else False
        if sig == full_sig:
            sig_match += 1
    if not vals:
        return {
            'n_loo': 0,
            'loo_slope_delta_min': np.nan,
            'loo_slope_delta_max': np.nan,
            'loo_slope_delta_median': np.nan,
            'loo_sign_agreement': np.nan,
            'loo_sig_agreement': np.nan,
        }
    arr = np.asarray(vals, dtype=float)
    return {
        'n_loo': int(n_ok),
        'loo_slope_delta_min': float(np.nanmin(arr)),
        'loo_slope_delta_max': float(np.nanmax(arr)),
        'loo_slope_delta_median': float(np.nanmedian(arr)),
        'loo_sign_agreement': float(sign_match / n_ok) if n_ok else np.nan,
        'loo_sig_agreement': float(sig_match / n_ok) if n_ok else np.nan,
    }


def _segment_quality(x: np.ndarray, y: np.ndarray, full_span: Optional[float] = None) -> Dict[str, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 2:
        return {'corr_abs': np.nan, 'sign_changes': np.nan, 'span': np.nan, 'score': -np.inf}
    corr = _pair_correlation(x, y)
    corr_abs = abs(corr) if np.isfinite(corr) else 0.0
    changes = float(_trend_sign_changes(y))
    span = float(x[-1] - x[0]) if len(x) >= 2 else 0.0
    if not np.isfinite(full_span) or full_span <= 0:
        span_frac = 1.0
    else:
        span_frac = float(np.clip(span / full_span, 0.0, 1.0))
    len_term = min(len(x), 10) / 10.0
    score = 1.75 * corr_abs + 0.45 * span_frac + 0.35 * len_term - 0.60 * changes
    return {'corr_abs': float(corr_abs), 'sign_changes': float(changes), 'span': float(span), 'score': float(score)}


def _select_stable_segment(df: pd.DataFrame, min_len: int = 4) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    base = df.copy().sort_values('pressure').reset_index(drop=True)
    base = base[np.isfinite(base['pressure']) & np.isfinite(base['value'])].reset_index(drop=True)
    meta = {
        'segment_used': False,
        'segment_reason': 'full_series',
        'segment_n': int(len(base)),
        'segment_pressure_min': float(np.nanmin(base['pressure'])) if len(base) else np.nan,
        'segment_pressure_max': float(np.nanmax(base['pressure'])) if len(base) else np.nan,
        'segment_score': np.nan,
        'segment_corr_abs': np.nan,
        'segment_sign_changes': np.nan,
        'full_n': int(len(base)),
        'full_pressure_min': float(np.nanmin(base['pressure'])) if len(base) else np.nan,
        'full_pressure_max': float(np.nanmax(base['pressure'])) if len(base) else np.nan,
        'full_score': np.nan,
        'full_corr_abs': np.nan,
        'full_sign_changes': np.nan,
    }
    if len(base) < max(min_len, 3):
        return base, meta
    x = base['pressure'].to_numpy(dtype=float)
    y = base['value'].to_numpy(dtype=float)
    full_span = float(x[-1] - x[0]) if len(x) >= 2 else np.nan
    full_q = _segment_quality(x, y, full_span=full_span)
    meta.update({
        'segment_score': full_q['score'],
        'segment_corr_abs': full_q['corr_abs'],
        'segment_sign_changes': full_q['sign_changes'],
        'full_score': full_q['score'],
        'full_corr_abs': full_q['corr_abs'],
        'full_sign_changes': full_q['sign_changes'],
    })
    best = {'i': 0, 'j': len(base), **full_q}
    for i in range(0, len(base) - min_len + 1):
        for j in range(i + min_len, len(base) + 1):
            cand = base.iloc[i:j].reset_index(drop=True)
            q = _segment_quality(cand['pressure'].to_numpy(dtype=float), cand['value'].to_numpy(dtype=float), full_span=full_span)
            better = q['score'] > best['score'] + 1e-12
            tied = abs(q['score'] - best['score']) <= 1e-12 and len(cand) > (best['j'] - best['i'])
            if better or tied:
                best = {'i': i, 'j': j, **q}
    should_consider = (
        (np.isfinite(full_q['sign_changes']) and full_q['sign_changes'] >= 2) or
        (np.isfinite(full_q['corr_abs']) and full_q['corr_abs'] < 0.70)
    )
    seg_len = int(best['j'] - best['i'])
    seg_span_ok = np.isfinite(full_span) and full_span > 0 and best['span'] >= 0.35 * full_span
    seg_len_ok = seg_len >= max(min_len, int(math.ceil(0.45 * len(base))))
    seg_better = best['score'] >= full_q['score'] + 0.20
    if should_consider and seg_len_ok and seg_span_ok and seg_better:
        out = base.iloc[best['i']:best['j']].reset_index(drop=True)
        meta.update({
            'segment_used': True,
            'segment_reason': 'stable_monotonic_subsegment',
            'segment_n': int(len(out)),
            'segment_pressure_min': float(np.nanmin(out['pressure'])) if len(out) else np.nan,
            'segment_pressure_max': float(np.nanmax(out['pressure'])) if len(out) else np.nan,
            'segment_score': best['score'],
            'segment_corr_abs': best['corr_abs'],
            'segment_sign_changes': best['sign_changes'],
        })
        return out, meta
    return base, meta


def _delta_stability_stats(pair_df: pd.DataFrame) -> Dict[str, Any]:
    if pair_df is None or len(pair_df) == 0 or 'delta_left_minus_right' not in pair_df.columns:
        return {
            'delta_sign_consistency': np.nan,
            'delta_majority_sign': 'mixed',
            'delta_sign_changes': np.nan,
            'delta_mean_z': np.nan,
            'delta_median': np.nan,
            'delta_iqr': np.nan,
            'delta_stable_offset': False,
        }
    delta = pair_df['delta_left_minus_right'].to_numpy(dtype=float)
    delta = delta[np.isfinite(delta)]
    if len(delta) == 0:
        return {
            'delta_sign_consistency': np.nan,
            'delta_majority_sign': 'mixed',
            'delta_sign_changes': np.nan,
            'delta_mean_z': np.nan,
            'delta_median': np.nan,
            'delta_iqr': np.nan,
            'delta_stable_offset': False,
        }
    absd = np.abs(delta)
    thr = max(2.0, float(np.nanmedian(absd) * 0.35)) if len(absd) else 2.0
    informative = absd >= thr
    if np.count_nonzero(informative) < min(3, len(delta)):
        informative = np.isfinite(delta)
    dv = delta[informative]
    pos = float(np.mean(dv > 0)) if len(dv) else np.nan
    neg = float(np.mean(dv < 0)) if len(dv) else np.nan
    sign_consistency = max(pos, neg) if (np.isfinite(pos) and np.isfinite(neg)) else np.nan
    if np.isfinite(pos) and np.isfinite(neg):
        if pos > neg:
            majority = 'positive'
        elif neg > pos:
            majority = 'negative'
        else:
            majority = 'mixed'
    else:
        majority = 'mixed'
    delta_changes = _trend_sign_changes(delta)
    if len(delta) >= 3:
        sd = float(np.nanstd(delta, ddof=1))
        sem = sd / math.sqrt(len(delta)) if np.isfinite(sd) and sd > 0 else np.nan
        mean_z = float(np.nanmean(delta) / sem) if np.isfinite(sem) and sem > 0 else np.nan
    else:
        mean_z = np.nan
    q25 = float(np.nanpercentile(delta, 25)) if len(delta) else np.nan
    q75 = float(np.nanpercentile(delta, 75)) if len(delta) else np.nan
    iqr = q75 - q25 if np.isfinite(q25) and np.isfinite(q75) else np.nan
    stable = bool(
        np.isfinite(sign_consistency) and sign_consistency >= 0.75 and
        delta_changes <= 1 and
        np.isfinite(mean_z) and abs(mean_z) >= 2.0
    )
    return {
        'delta_sign_consistency': float(sign_consistency) if np.isfinite(sign_consistency) else np.nan,
        'delta_majority_sign': majority,
        'delta_sign_changes': int(delta_changes),
        'delta_mean_z': mean_z,
        'delta_median': float(np.nanmedian(delta)) if len(delta) else np.nan,
        'delta_iqr': float(iqr) if np.isfinite(iqr) else np.nan,
        'delta_stable_offset': stable,
    }


def _degrade_grade(grade: str, steps: int = 1) -> str:
    order = ['high', 'medium', 'low']
    g = str(grade).lower()
    if g not in order:
        g = 'medium'
    idx = order.index(g)
    idx = min(idx + max(int(steps), 0), len(order) - 1)
    return order[idx]


def _classify_comparison(
    pair_df: pd.DataFrame,
    left_fit: Optional[Dict[str, float]],
    right_fit: Optional[Dict[str, float]],
    fallback_left: float,
    fallback_right: float,
    region_kind: str = 'main_band',
    delta_priority: bool = False,
) -> Dict[str, Any]:
    n_pairs = int(len(pair_df))
    if n_pairs == 0 or left_fit is None or right_fit is None:
        return {
            'evidence_grade': 'low',
            'diagnostic_label': 'insufficient_overlap',
            'recommended_use': 'inspect_manually',
            **_delta_stability_stats(pair_df),
            **{
                'pair_correlation': np.nan,
                'slope_z_score': np.nan,
                'trend_sign_changes_left': np.nan,
                'trend_sign_changes_right': np.nan,
                'n_loo': 0,
                'loo_slope_delta_min': np.nan,
                'loo_slope_delta_max': np.nan,
                'loo_slope_delta_median': np.nan,
                'loo_sign_agreement': np.nan,
                'loo_sig_agreement': np.nan,
            },
        }
    x = pair_df['pressure'].to_numpy(dtype=float)
    yl = pair_df['left_value'].to_numpy(dtype=float)
    yr = pair_df['right_value'].to_numpy(dtype=float)
    yls = pair_df['left_value_se'].to_numpy(dtype=float)
    yrs = pair_df['right_value_se'].to_numpy(dtype=float)
    delta = yl - yr
    mean_delta = float(np.nanmean(delta)) if n_pairs else np.nan
    pair_corr = _pair_correlation(yl, yr)
    slope_delta = float(left_fit['slope'] - right_fit['slope'])
    slope_scale = math.sqrt(
        (left_fit.get('slope_se', np.nan) ** 2 if np.isfinite(left_fit.get('slope_se', np.nan)) else 0.0) +
        (right_fit.get('slope_se', np.nan) ** 2 if np.isfinite(right_fit.get('slope_se', np.nan)) else 0.0)
    )
    slope_z = float(slope_delta / slope_scale) if np.isfinite(slope_scale) and slope_scale > 0 else np.nan
    changes_left = _trend_sign_changes(yl)
    changes_right = _trend_sign_changes(yr)
    loo = _leave_one_out_slope_deltas(x, yl, yls, yr, yrs)
    delta_stats = _delta_stability_stats(pair_df)

    grade = 'high'
    if n_pairs < 5:
        grade = 'medium'
    if n_pairs < 4:
        grade = 'low'
    max_fb = max(float(fallback_left), float(fallback_right))
    if max_fb > 0.10:
        grade = _degrade_grade(grade, 1)
    if max_fb > 0.25:
        grade = _degrade_grade(grade, 1)
    if np.isfinite(pair_corr):
        if pair_corr < 0.70:
            grade = _degrade_grade(grade, 1)
        if pair_corr < 0.40:
            grade = _degrade_grade(grade, 1)
    else:
        grade = _degrade_grade(grade, 1)
    if np.isfinite(loo.get('loo_sign_agreement', np.nan)) and loo['loo_sign_agreement'] < 0.75:
        grade = _degrade_grade(grade, 1)
    if np.isfinite(loo.get('loo_sig_agreement', np.nan)) and loo['loo_sig_agreement'] < 0.60:
        grade = _degrade_grade(grade, 1)

    nonlinear = (changes_left >= 2) or (changes_right >= 2)
    if region_kind == 'branch_sensitive' and nonlinear:
        grade = _degrade_grade(grade, 1)

    if delta_priority:
        if not bool(delta_stats.get('delta_stable_offset', False)):
            grade = _degrade_grade(grade, 1)
        if np.isfinite(delta_stats.get('delta_sign_consistency', np.nan)) and delta_stats['delta_sign_consistency'] < 0.70:
            grade = _degrade_grade(grade, 1)

    if n_pairs < 4:
        label = 'insufficient_overlap'
    elif region_kind == 'branch_sensitive' and (nonlinear or (np.isfinite(pair_corr) and pair_corr < 0.55)):
        label = 'nonlinear_or_branch_switching'
    elif delta_priority and bool(delta_stats.get('delta_stable_offset', False)) and np.isfinite(delta_stats.get('delta_sign_consistency', np.nan)) and delta_stats['delta_sign_consistency'] >= 0.75:
        label = 'hysteresis_like_offset'
    elif np.isfinite(slope_z) and abs(slope_z) >= 2.0 and np.isfinite(loo.get('loo_sign_agreement', np.nan)) and loo['loo_sign_agreement'] >= 0.75 and (not delta_priority or np.isfinite(delta_stats.get('delta_sign_consistency', np.nan)) and delta_stats['delta_sign_consistency'] >= 0.60):
        label = 'slope_difference'
    elif np.isfinite(mean_delta) and abs(mean_delta) >= 5.0 and (not np.isfinite(slope_z) or abs(slope_z) < 2.0) and (not np.isfinite(pair_corr) or pair_corr >= 0.60):
        label = 'offset_with_similar_slope'
    elif np.isfinite(pair_corr) and pair_corr >= 0.70 and (not np.isfinite(slope_z) or abs(slope_z) < 2.0):
        label = 'trend_consistent'
    else:
        label = 'ambiguous'

    recommended = 'main_text'
    if grade == 'medium' or label in {'offset_with_similar_slope', 'ambiguous', 'hysteresis_like_offset'}:
        recommended = 'supporting'
    if grade == 'low' or label in {'nonlinear_or_branch_switching', 'insufficient_overlap'}:
        recommended = 'inspect_manually'

    return {
        'evidence_grade': grade,
        'diagnostic_label': label,
        'pair_correlation': pair_corr,
        'slope_z_score': slope_z,
        'trend_sign_changes_left': changes_left,
        'trend_sign_changes_right': changes_right,
        'recommended_use': recommended,
        **delta_stats,
        **loo,
    }


def _comparison_metric_ylabel(metric: str) -> str:
    metric_l = str(metric).strip().lower()
    if 'center' in metric_l or 'position' in metric_l:
        return 'Peak position (cm$^{-1}$)'
    if 'fwhm' in metric_l or 'width' in metric_l:
        return 'FWHM (cm$^{-1}$)'
    if 'sigma' in metric_l:
        return r'Peak width $\sigma$ (cm$^{-1}$)'
    if 'area' in metric_l:
        return 'Peak area (a.u.)'
    if 'amplitude' in metric_l or 'height' in metric_l:
        return 'Peak amplitude (a.u.)'
    return 'Metric value (a.u.)'


def _comparison_no_tracking_outpath(outpath: Path) -> Path:
    return outpath.with_name(f'{outpath.stem}_no_tracking_axes{outpath.suffix}')


def plot_comparison_metric(

    pair_df: pd.DataFrame,
    left_raw: pd.DataFrame,
    right_raw: pd.DataFrame,
    left_fit: Optional[Dict[str, float]],
    right_fit: Optional[Dict[str, float]],
    outpath: Path,
    analysis_name: str,
    metric: str,
    left_label: str,
    right_label: str,
    extra_text: Optional[str] = None,
    dpi: int = 300,
    segmented_fit_breaks: Optional[Dict[str, Sequence[float]]] = None,
) -> None:
    ensure_dir(outpath.parent)

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.2), sharex=True, gridspec_kw={'height_ratios': [3.1, 1.2]})
    ax = axes[0]
    ax2 = axes[1]

    styles = {
        'left': {'color': '#1f77b4', 'marker': 'o', 'label': left_label},
        'right': {'color': '#C44E52', 'marker': 's', 'label': right_label},
    }

    def _plot_series(
        axh: plt.Axes,
        frame: pd.DataFrame,
        fit: Optional[Dict[str, float]],
        style: Dict[str, Any],
        series_key: str,
    ) -> None:
        x = frame['pressure'].to_numpy(dtype=float)
        y = frame['value'].to_numpy(dtype=float)
        yerr = frame['value_se'].to_numpy(dtype=float) if 'value_se' in frame.columns else np.full_like(y, np.nan, dtype=float)
        finite = np.isfinite(x) & np.isfinite(y)
        finite_err = finite & np.isfinite(yerr) & (yerr > 0)

        if np.any(finite_err):
            axh.errorbar(
                x[finite_err], y[finite_err], yerr=yerr[finite_err],
                fmt=style['marker'], linestyle='none', markersize=5.6,
                markerfacecolor='white', markeredgecolor=style['color'], markeredgewidth=1.2,
                ecolor=style['color'], elinewidth=1.0, capsize=2.5, label=style['label'],
            )
        if np.any(finite & ~finite_err):
            axh.plot(
                x[finite & ~finite_err], y[finite & ~finite_err],
                linestyle='none', marker=style['marker'], markersize=5.2,
                markerfacecolor='white', markeredgecolor=style['color'], markeredgewidth=1.2,
                label=None,
            )

        segment_breaks = []
        if segmented_fit_breaks is not None:
            segment_breaks = [float(v) for v in segmented_fit_breaks.get(series_key, []) if np.isfinite(v)]

        if segment_breaks and np.isfinite(x).sum() >= 2:
            edges = [float(np.nanmin(x[finite])), *sorted(segment_breaks), float(np.nanmax(x[finite]))]
            for idx in range(len(edges) - 1):
                lo = edges[idx]
                hi = edges[idx + 1]
                if idx == 0:
                    seg_mask = finite & (x <= hi)
                elif idx == len(edges) - 2:
                    seg_mask = finite & (x > lo)
                else:
                    seg_mask = finite & (x > lo) & (x <= hi)

                xseg = x[seg_mask]
                yseg = y[seg_mask]
                yerr_seg = yerr[seg_mask]
                if np.count_nonzero(np.isfinite(xseg) & np.isfinite(yseg)) < 2:
                    continue
                seg_fit = robust_linear_fit(xseg, yseg, yerr_seg)
                if seg_fit is None:
                    continue
                xseg_fit = xseg[np.isfinite(xseg)]
                if np.unique(xseg_fit).size < 2:
                    continue
                xx = np.linspace(np.nanmin(xseg_fit), np.nanmax(xseg_fit), 100)
                yy = seg_fit['slope'] * xx + seg_fit['intercept']
                axh.plot(xx, yy, color=style['color'], lw=1.6, ls='-', alpha=0.95)
        elif fit is not None and np.isfinite(x).sum() >= 2:
            xfit = x[np.isfinite(x)]
            if np.unique(xfit).size >= 2:
                xx = np.linspace(np.nanmin(xfit), np.nanmax(xfit), 200)
                yy = fit['slope'] * xx + fit['intercept']
                axh.plot(xx, yy, color=style['color'], lw=1.6, ls='-', alpha=0.95)

    _plot_series(ax, left_raw, left_fit, styles['left'], 'left')
    _plot_series(ax, right_raw, right_fit, styles['right'], 'right')
    ax.legend(frameon=False, fontsize=8, loc='best')
    ax.set_ylabel(_comparison_metric_ylabel(metric))
    _clean_axes(ax)

    delta_x = pair_df['pressure'].to_numpy(dtype=float)
    delta_y = pair_df['delta_left_minus_right'].to_numpy(dtype=float)
    ax2.axhline(0.0, color='0.55', lw=0.9)
    finite_delta = np.isfinite(delta_x) & np.isfinite(delta_y)
    if np.any(finite_delta):
        ax2.plot(
            delta_x[finite_delta], delta_y[finite_delta],
            linestyle='-', lw=1.0, color='black',
            marker='o', markersize=4.2, markerfacecolor='white', markeredgecolor='black', markeredgewidth=0.9,
        )
    ax2.set_xlabel('Pressure (GPa)')
    ax2.set_ylabel(r'$\Delta$')
    _clean_axes(ax2)

    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    plain_outpath = _comparison_no_tracking_outpath(outpath)
    fig_plain, ax_plain = plt.subplots(1, 1, figsize=(7.2, 5.2))
    _plot_series(ax_plain, left_raw, left_fit, styles['left'], 'left')
    _plot_series(ax_plain, right_raw, right_fit, styles['right'], 'right')
    ax_plain.legend(frameon=False, fontsize=8, loc='best')
    ax_plain.set_xlabel('Pressure (GPa)')
    ax_plain.set_ylabel(_comparison_metric_ylabel(metric))
    _clean_axes(ax_plain)
    fig_plain.tight_layout()
    fig_plain.savefig(plain_outpath, dpi=dpi, bbox_inches='tight')
    plt.close(fig_plain)


def compare_peak_tables(
    left_peak: pd.DataFrame,
    right_peak: pd.DataFrame,
    analysis_name: str,
    left_profile: str,
    right_profile: str,
    outdir: Path,
    metrics: Sequence[str],
    region_kind: str = 'main_band',
    delta_priority: bool = False,
    segment_branch_sensitive: bool = False,
    dpi: int = 300,
    progress: bool = True,
    segmented_fit_breaks: Optional[Dict[str, Sequence[float]]] = None,
) -> Dict[str, Any]:
    ensure_dir(outdir)
    left_df = _prepare_peak_table_for_compare(left_peak)
    right_df = _prepare_peak_table_for_compare(right_peak)

    overall_rows: List[Dict[str, Any]] = []
    pair_sheets: Dict[str, pd.DataFrame] = {}
    available_metrics = [m for m in metrics if m in left_df.columns and m in right_df.columns]
    fallback_left = float(np.mean(left_df.get('mode_used', '').astype(str).eq('max'))) if len(left_df) else np.nan
    fallback_right = float(np.mean(right_df.get('mode_used', '').astype(str).eq('max'))) if len(right_df) else np.nan

    metric_iter = progress_iter(available_metrics, desc=f'{analysis_name} metrics', total=len(available_metrics), enabled=progress, leave=False)
    for metric in metric_iter:
        lmet_raw = _get_metric_frame(left_df, metric)
        rmet_raw = _get_metric_frame(right_df, metric)
        left_segment_meta = {
            'segment_used': False, 'segment_reason': 'full_series', 'segment_n': int(len(lmet_raw)),
            'segment_pressure_min': float(np.nanmin(lmet_raw['pressure'])) if len(lmet_raw) else np.nan,
            'segment_pressure_max': float(np.nanmax(lmet_raw['pressure'])) if len(lmet_raw) else np.nan,
            'segment_score': np.nan, 'segment_corr_abs': np.nan, 'segment_sign_changes': np.nan,
            'full_n': int(len(lmet_raw)), 'full_pressure_min': float(np.nanmin(lmet_raw['pressure'])) if len(lmet_raw) else np.nan,
            'full_pressure_max': float(np.nanmax(lmet_raw['pressure'])) if len(lmet_raw) else np.nan,
            'full_score': np.nan, 'full_corr_abs': np.nan, 'full_sign_changes': np.nan,
        }
        right_segment_meta = {
            'segment_used': False, 'segment_reason': 'full_series', 'segment_n': int(len(rmet_raw)),
            'segment_pressure_min': float(np.nanmin(rmet_raw['pressure'])) if len(rmet_raw) else np.nan,
            'segment_pressure_max': float(np.nanmax(rmet_raw['pressure'])) if len(rmet_raw) else np.nan,
            'segment_score': np.nan, 'segment_corr_abs': np.nan, 'segment_sign_changes': np.nan,
            'full_n': int(len(rmet_raw)), 'full_pressure_min': float(np.nanmin(rmet_raw['pressure'])) if len(rmet_raw) else np.nan,
            'full_pressure_max': float(np.nanmax(rmet_raw['pressure'])) if len(rmet_raw) else np.nan,
            'full_score': np.nan, 'full_corr_abs': np.nan, 'full_sign_changes': np.nan,
        }

        lmet = lmet_raw.copy()
        rmet = rmet_raw.copy()
        if segment_branch_sensitive and region_kind == 'branch_sensitive':
            lmet, left_segment_meta = _select_stable_segment(lmet_raw)
            rmet, right_segment_meta = _select_stable_segment(rmet_raw)

        pair_df = _pair_series_on_grid(lmet, rmet)
        if len(pair_df) < 4 and (left_segment_meta.get('segment_used') or right_segment_meta.get('segment_used')):
            lmet = lmet_raw.copy()
            rmet = rmet_raw.copy()
            pair_df = _pair_series_on_grid(lmet, rmet)
            left_segment_meta['segment_used'] = False
            left_segment_meta['segment_reason'] = 'segment_reverted_low_overlap'
            right_segment_meta['segment_used'] = False
            right_segment_meta['segment_reason'] = 'segment_reverted_low_overlap'

        left_fit = robust_linear_fit(lmet['pressure'].to_numpy(dtype=float), lmet['value'].to_numpy(dtype=float), lmet['value_se'].to_numpy(dtype=float))
        right_fit = robust_linear_fit(rmet['pressure'].to_numpy(dtype=float), rmet['value'].to_numpy(dtype=float), rmet['value_se'].to_numpy(dtype=float))
        cls = _classify_comparison(
            pair_df,
            left_fit,
            right_fit,
            fallback_left,
            fallback_right,
            region_kind=region_kind,
            delta_priority=delta_priority,
        )

        overlap_min = float(np.nanmin(pair_df['pressure'])) if len(pair_df) else np.nan
        overlap_max = float(np.nanmax(pair_df['pressure'])) if len(pair_df) else np.nan
        overlap_span = float(overlap_max - overlap_min) if np.isfinite(overlap_min) and np.isfinite(overlap_max) else np.nan
        delta = pair_df['delta_left_minus_right'].to_numpy(dtype=float) if len(pair_df) else np.array([], dtype=float)
        row = {
            'analysis': analysis_name,
            'metric': metric,
            'left_profile': left_profile,
            'right_profile': right_profile,
            'left_label': left_profile,
            'right_label': right_profile,
            'overlap_pressure_min': overlap_min,
            'overlap_pressure_max': overlap_max,
            'overlap_pressure_span': overlap_span,
            'n_left_points': int(len(lmet)),
            'n_right_points': int(len(rmet)),
            'n_pairs': int(len(pair_df)),
            'left_fit_available': bool(left_fit is not None),
            'right_fit_available': bool(right_fit is not None),
            'mean_delta_left_minus_right': float(np.nanmean(delta)) if len(delta) else np.nan,
            'mean_abs_delta': float(np.nanmean(np.abs(delta))) if len(delta) else np.nan,
            'std_delta': float(np.nanstd(delta, ddof=1)) if len(delta) >= 2 else np.nan,
            'rmse_delta': float(np.sqrt(np.nanmean(delta ** 2))) if len(delta) else np.nan,
            'pair_correlation': cls.get('pair_correlation', np.nan),
            'fallback_fraction_left': fallback_left,
            'fallback_fraction_right': fallback_right,
            'trend_sign_changes_left': cls.get('trend_sign_changes_left', np.nan),
            'trend_sign_changes_right': cls.get('trend_sign_changes_right', np.nan),
            'delta_sign_consistency': cls.get('delta_sign_consistency', np.nan),
            'delta_majority_sign': cls.get('delta_majority_sign', 'mixed'),
            'delta_sign_changes': cls.get('delta_sign_changes', np.nan),
            'delta_mean_z': cls.get('delta_mean_z', np.nan),
            'delta_median': cls.get('delta_median', np.nan),
            'delta_iqr': cls.get('delta_iqr', np.nan),
            'delta_stable_offset': cls.get('delta_stable_offset', False),
            'n_loo': cls.get('n_loo', np.nan),
            'loo_sign_agreement': cls.get('loo_sign_agreement', np.nan),
            'loo_sig_agreement': cls.get('loo_sig_agreement', np.nan),
            'evidence_grade': cls.get('evidence_grade', 'low'),
            'diagnostic_label': cls.get('diagnostic_label', 'ambiguous'),
            'recommended_use': cls.get('recommended_use', 'inspect_manually'),
            'region_kind': region_kind,
            'delta_priority': bool(delta_priority),
            'segment_branch_sensitive': bool(segment_branch_sensitive),
            'left_segment_used': bool(left_segment_meta.get('segment_used', False)),
            'left_segment_reason': left_segment_meta.get('segment_reason', 'full_series'),
            'left_segment_n': left_segment_meta.get('segment_n', np.nan),
            'left_segment_pressure_min': left_segment_meta.get('segment_pressure_min', np.nan),
            'left_segment_pressure_max': left_segment_meta.get('segment_pressure_max', np.nan),
            'left_segment_score': left_segment_meta.get('segment_score', np.nan),
            'left_full_n': left_segment_meta.get('full_n', np.nan),
            'left_full_score': left_segment_meta.get('full_score', np.nan),
            'left_full_corr_abs': left_segment_meta.get('full_corr_abs', np.nan),
            'left_full_sign_changes': left_segment_meta.get('full_sign_changes', np.nan),
            'right_segment_used': bool(right_segment_meta.get('segment_used', False)),
            'right_segment_reason': right_segment_meta.get('segment_reason', 'full_series'),
            'right_segment_n': right_segment_meta.get('segment_n', np.nan),
            'right_segment_pressure_min': right_segment_meta.get('segment_pressure_min', np.nan),
            'right_segment_pressure_max': right_segment_meta.get('segment_pressure_max', np.nan),
            'right_segment_score': right_segment_meta.get('segment_score', np.nan),
            'right_full_n': right_segment_meta.get('full_n', np.nan),
            'right_full_score': right_segment_meta.get('full_score', np.nan),
            'right_full_corr_abs': right_segment_meta.get('full_corr_abs', np.nan),
            'right_full_sign_changes': right_segment_meta.get('full_sign_changes', np.nan),
        }
        overall_rows.append(row)

        plot_comparison_metric(
            pair_df=pair_df,
            left_raw=lmet,
            right_raw=rmet,
            left_fit=left_fit,
            right_fit=right_fit,
            outpath=outdir / f'{analysis_name}_{metric}.png',
            analysis_name=analysis_name,
            metric=metric,
            left_label=left_profile,
            right_label=right_profile,
            extra_text=None,
            dpi=dpi,
            segmented_fit_breaks=segmented_fit_breaks,
        )
        pair_df = pair_df.copy()
        pair_df['left_profile'] = left_profile
        pair_df['right_profile'] = right_profile
        pair_df['metric'] = metric
        pair_df['analysis'] = analysis_name
        pair_df['delta_priority'] = bool(delta_priority)
        pair_df['left_segment_used'] = bool(left_segment_meta.get('segment_used', False))
        pair_df['right_segment_used'] = bool(right_segment_meta.get('segment_used', False))
        pair_sheets[f'Pairs_{analysis_name}_{metric}'] = pair_df

    return {
        'summary': pd.DataFrame(overall_rows),
        'pair_sheets': pair_sheets,
        'manifest_row': {
            'analysis': analysis_name,
            'left_profile': left_profile,
            'right_profile': right_profile,
            'output_dir': str(outdir.resolve()),
        },
    }


def run_comparison_group(

    group_cfg: Dict[str, Any],
    profiles_dict: Dict[str, Dict[str, Any]],
    excel_map: Dict[str, str],
    output_dir: Path,
) -> Dict[str, Any]:
    ensure_dir(output_dir)
    analyses = group_cfg.get('analyses', [])
    dpi = int(group_cfg.get('dpi', 300))
    progress_enabled = bool(group_cfg.get('progress', True))
    processed_cache: Dict[str, Dict[str, Any]] = {}
    manifest_rows: List[Dict[str, Any]] = []
    overall_tables: List[pd.DataFrame] = []
    analysis_pair_sheets: Dict[str, pd.DataFrame] = {}

    def _run_profile(profile_name: str) -> Dict[str, Any]:
        if profile_name in processed_cache:
            return processed_cache[profile_name]
        if profile_name not in profiles_dict:
            raise KeyError(f'Unknown profile in comparison group: {profile_name}')
        cfg = deep_copy_config(profiles_dict[profile_name])
        cfg['progress'] = progress_enabled
        if profile_name in excel_map:
            cfg['excel_path'] = excel_map[profile_name]
        if not cfg.get('excel_path'):
            raise ValueError(f'Missing excel_path for profile {profile_name}. Provide it in the profile or via --excel-map.')
        result = process_profile(cfg)
        processed_cache[profile_name] = result
        return result

    analysis_iter = progress_iter(analyses, desc=f"{group_cfg.get('name', 'comparison')} analyses", total=len(analyses), enabled=progress_enabled, leave=True)
    for item in analysis_iter:
        analysis_name = str(item['name'])
        left_profile = str(item['left_profile'])
        right_profile = str(item['right_profile'])
        metrics = list(item.get('metrics', ['effective_center']))
        region_kind = str(item.get('region_kind', 'main_band'))
        delta_priority = bool(item.get('delta_priority', False))
        segment_branch_sensitive = bool(item.get('segment_branch_sensitive', False))
        left_res = _run_profile(left_profile)
        right_res = _run_profile(right_profile)
        left_peak = left_res.get('peak_table', pd.DataFrame())
        right_peak = right_res.get('peak_table', pd.DataFrame())
        analysis_outdir = output_dir / analysis_name
        cmp_res = compare_peak_tables(
            left_peak=left_peak,
            right_peak=right_peak,
            analysis_name=analysis_name,
            left_profile=left_profile,
            right_profile=right_profile,
            outdir=analysis_outdir,
            metrics=metrics,
            region_kind=region_kind,
            delta_priority=delta_priority,
            segment_branch_sensitive=segment_branch_sensitive,
            dpi=dpi,
            progress=progress_enabled,
            segmented_fit_breaks=item.get('segmented_fit_breaks'),
        )
        manifest_rows.append(cmp_res['manifest_row'])
        overall_tables.append(cmp_res['summary'])
        for key, value in cmp_res['pair_sheets'].items():
            analysis_pair_sheets[key] = value

    manifest_df = pd.DataFrame(manifest_rows)
    overall_df = pd.concat(overall_tables, ignore_index=True) if overall_tables else pd.DataFrame()
    summary_path = output_dir / 'comparison_summary.xlsx'
    with pd.ExcelWriter(summary_path, engine='openpyxl') as writer:
        manifest_df.to_excel(writer, sheet_name='Manifest', index=False)
        overall_df.to_excel(writer, sheet_name='OverallSummary', index=False)
        if not overall_df.empty:
            for analysis_name in manifest_df.get('analysis', pd.Series(dtype=str)).tolist():
                sub = overall_df[overall_df['analysis'] == analysis_name].copy()
                sub.to_excel(writer, sheet_name=f'Summary_{analysis_name}'[:31], index=False)
        for sheet_name, df_sheet in analysis_pair_sheets.items():
            df_sheet.to_excel(writer, sheet_name=sheet_name[:31], index=False)

    log_payload = {
        'group_name': group_cfg.get('name'),
        'description': group_cfg.get('description', ''),
        'output_dir': str(output_dir.resolve()),
        'analyses': manifest_rows,
        'timestamp': pd.Timestamp.now().isoformat(),
    }
    with open(output_dir / 'comparison_run_log.json', 'w', encoding='utf-8') as f:
        json.dump(log_payload, f, ensure_ascii=False, indent=2)

    return {
        'output_dir': str(output_dir.resolve()),
        'summary_path': str(summary_path.resolve()),
        'manifest': manifest_df,
        'overall_summary': overall_df,
    }
