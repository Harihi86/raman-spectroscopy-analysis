#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class PhotoactivityResult:
    config: Dict[str, Any]
    wavelength_nm: np.ndarray
    wavenumber_cm1: np.ndarray
    counts: np.ndarray
    counts_per_s: np.ndarray
    time_s: np.ndarray
    time_norm: np.ndarray
    frame_index: np.ndarray
    peak_metric: np.ndarray
    background_metric: np.ndarray
    mean_spectrum_counts_per_s: np.ndarray
    early_mean_counts_per_s: np.ndarray
    late_mean_counts_per_s: np.ndarray
    dominant_peak_center_nm: float
    dominant_peak_center_cm1: float
    peak_change_pct: float
    background_change_pct: float


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_time_series_excel(excel_path: str) -> Tuple[np.ndarray, np.ndarray]:
    df = pd.read_excel(excel_path, header=None)
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if df.shape[1] < 2:
        raise ValueError(f"Expected first column = spectral axis and subsequent frame columns; got shape {df.shape}")

    wavelength = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
    counts = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    finite_row = np.isfinite(wavelength) & np.isfinite(counts).any(axis=1)
    wavelength = wavelength[finite_row]
    counts = counts[finite_row, :]

    order = np.argsort(wavelength)
    wavelength = wavelength[order]
    counts = counts[order, :]
    return wavelength, counts


def wavelength_to_wavenumber_cm1(wavelength_nm: np.ndarray, excitation_nm: float) -> np.ndarray:
    wavelength_nm = np.asarray(wavelength_nm, dtype=float)
    excitation_nm = float(excitation_nm)
    if excitation_nm <= 0:
        raise ValueError("excitation_wavelength_nm must be > 0")
    if np.any(wavelength_nm <= 0):
        raise ValueError("All spectral-axis wavelengths must be > 0 nm")
    return 1.0e7 * (1.0 / excitation_nm - 1.0 / wavelength_nm)


def nm_window_to_cm1(window_nm: Sequence[float], excitation_nm: float) -> List[float]:
    w = wavelength_to_wavenumber_cm1(np.asarray(window_nm, dtype=float), excitation_nm)
    return [float(np.min(w)), float(np.max(w))]


def frame_slice(frame_window_1_based: Sequence[int]) -> slice:
    start, end = int(frame_window_1_based[0]), int(frame_window_1_based[1])
    if start < 1 or end < start:
        raise ValueError(f"Invalid frame window: {frame_window_1_based}")
    return slice(start - 1, end)


def format_window(window: Sequence[float]) -> str:
    return f"{float(window[0]):.1f}-{float(window[1]):.1f}"


def integrate_window(x: np.ndarray, y2d: np.ndarray, window: Sequence[float]) -> np.ndarray:
    lo, hi = float(window[0]), float(window[1])
    mask = (x >= lo) & (x <= hi)
    if np.count_nonzero(mask) < 2:
        raise ValueError(f"Integration window {window} contains fewer than 2 spectral points.")
    return np.trapz(y2d[mask, :], x=x[mask], axis=0)


def mean_window(x: np.ndarray, y2d: np.ndarray, window: Sequence[float]) -> np.ndarray:
    lo, hi = float(window[0]), float(window[1])
    mask = (x >= lo) & (x <= hi)
    if np.count_nonzero(mask) < 1:
        raise ValueError(f"Background window {window} contains no spectral points.")
    return np.nanmean(y2d[mask, :], axis=0)


def _norm_to_early(metric: np.ndarray, early_sl: slice) -> np.ndarray:
    early_mean = float(np.nanmean(metric[early_sl]))
    if not np.isfinite(early_mean) or early_mean == 0:
        return np.full_like(metric, np.nan, dtype=float)
    return metric / early_mean


def analyze_single_dataset(cfg: Dict[str, Any]) -> PhotoactivityResult:
    wavelength_nm, counts = read_time_series_excel(cfg["excel_path"])
    _n_spectral, n_frames = counts.shape
    expected = int(cfg.get("n_frames_expected", n_frames))
    if n_frames != expected:
        raise ValueError(f"{cfg['dataset_name']}: expected {expected} frames, found {n_frames}")

    exposure_s = float(cfg["exposure_s_per_frame"])
    if exposure_s <= 0:
        raise ValueError("exposure_s_per_frame must be > 0")

    excitation_nm = float(cfg["excitation_wavelength_nm"])
    wavenumber_cm1 = wavelength_to_wavenumber_cm1(wavelength_nm, excitation_nm)
    order = np.argsort(wavenumber_cm1)
    wavenumber_cm1 = wavenumber_cm1[order]
    wavelength_nm = wavelength_nm[order]
    counts = counts[order, :]

    counts_per_s = counts / exposure_s
    frame_index = np.arange(1, n_frames + 1, dtype=int)
    time_s = (frame_index.astype(float) - 1.0) * exposure_s
    time_norm = time_s / float(time_s[-1]) if float(time_s[-1]) > 0 else np.zeros_like(time_s)

    peak_window_cm1 = cfg.get("peak_window_cm1")
    if peak_window_cm1 is None:
        peak_window_cm1 = nm_window_to_cm1(cfg["peak_window_nm"], excitation_nm)
    tail_window_cm1 = cfg.get("tail_window_cm1")
    if tail_window_cm1 is None:
        tail_window_cm1 = nm_window_to_cm1(cfg["tail_window_nm"], excitation_nm)

    peak_metric = integrate_window(wavenumber_cm1, counts_per_s, peak_window_cm1)
    background_metric = mean_window(wavenumber_cm1, counts_per_s, tail_window_cm1)

    mean_spectrum = np.nanmean(counts_per_s, axis=1)
    peak_lo, peak_hi = [float(v) for v in peak_window_cm1]
    peak_mask = (wavenumber_cm1 >= peak_lo) & (wavenumber_cm1 <= peak_hi)
    dominant_peak_center_cm1 = float(wavenumber_cm1[peak_mask][np.nanargmax(mean_spectrum[peak_mask])])
    dominant_peak_center_nm = float(wavelength_nm[peak_mask][np.nanargmax(mean_spectrum[peak_mask])])

    early_sl = frame_slice(cfg["frame_windows"]["early"])
    late_sl = frame_slice(cfg["frame_windows"]["late"])
    early_mean_spectrum = np.nanmean(counts_per_s[:, early_sl], axis=1)
    late_mean_spectrum = np.nanmean(counts_per_s[:, late_sl], axis=1)

    peak_early = float(np.nanmean(peak_metric[early_sl]))
    peak_late = float(np.nanmean(peak_metric[late_sl]))
    bg_early = float(np.nanmean(background_metric[early_sl]))
    bg_late = float(np.nanmean(background_metric[late_sl]))

    peak_change_pct = np.nan if peak_early == 0 else 100.0 * (peak_late - peak_early) / peak_early
    background_change_pct = np.nan if bg_early == 0 else 100.0 * (bg_late - bg_early) / bg_early

    cfg = dict(cfg)
    cfg["peak_window_cm1"] = [float(peak_window_cm1[0]), float(peak_window_cm1[1])]
    cfg["tail_window_cm1"] = [float(tail_window_cm1[0]), float(tail_window_cm1[1])]

    return PhotoactivityResult(
        config=cfg,
        wavelength_nm=wavelength_nm,
        wavenumber_cm1=wavenumber_cm1,
        counts=counts,
        counts_per_s=counts_per_s,
        time_s=time_s,
        time_norm=time_norm,
        frame_index=frame_index,
        peak_metric=peak_metric,
        background_metric=background_metric,
        mean_spectrum_counts_per_s=mean_spectrum,
        early_mean_counts_per_s=early_mean_spectrum,
        late_mean_counts_per_s=late_mean_spectrum,
        dominant_peak_center_nm=dominant_peak_center_nm,
        dominant_peak_center_cm1=dominant_peak_center_cm1,
        peak_change_pct=float(peak_change_pct),
        background_change_pct=float(background_change_pct),
    )


def _clean_axes(ax: plt.Axes) -> None:
    ax.grid(False)
    ax.minorticks_on()
    ax.tick_params(which="major", direction="in", length=5, width=1.0)
    ax.tick_params(which="minor", direction="in", length=3, width=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def _quantile_limits(z: np.ndarray, qlo: float = 0.01, qhi: float = 0.995) -> Tuple[float, float]:
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        return 0.0, 1.0
    vmin, vmax = np.quantile(finite, [qlo, qhi])
    if vmax <= vmin:
        vmax = vmin + 1.0
    return float(vmin), float(vmax)


def _plot_heatmap(ax: plt.Axes, result: PhotoactivityResult, title: str, vmin: float | None = None, vmax: float | None = None, xlim: Sequence[float] | None = None, show_windows: bool = True):
    z = result.counts_per_s.T
    if vmin is None or vmax is None:
        vmin, vmax = _quantile_limits(z)
    extent = [
        float(result.wavenumber_cm1[0]),
        float(result.wavenumber_cm1[-1]),
        float(result.frame_index[0]),
        float(result.frame_index[-1]),
    ]
    im = ax.imshow(
        z,
        aspect="auto",
        origin="lower",
        extent=extent,
        interpolation="nearest",
        vmin=vmin,
        vmax=vmax,
        cmap="viridis",
    )
    if show_windows:
        peak_lo, peak_hi = result.config["peak_window_cm1"]
        tail_lo, tail_hi = result.config["tail_window_cm1"]
        ax.axvspan(peak_lo, peak_hi, linestyle="--", linewidth=1.0, fill=False)
        ax.axvspan(tail_lo, tail_hi, linestyle=":", linewidth=1.0, fill=False)
    if xlim is not None:
        ax.set_xlim(float(xlim[0]), float(xlim[1]))
    ax.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Time (frame)")
    ax.set_title(title)
    _clean_axes(ax)
    return im


def plot_single_heatmap(result: PhotoactivityResult, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    im = _plot_heatmap(ax, result, f"{result.config['dataset_name']} heatmap")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Intensity (counts s$^{-1}$)")
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_exp2_two_pressure_heatmap(results: Sequence[PhotoactivityResult], outpath: Path) -> None:
    if len(results) != 2:
        raise ValueError("exp2 comparison expects exactly two datasets")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.3), sharey=True, constrained_layout=True)
    stacked = np.hstack([r.counts_per_s.T.ravel() for r in results])
    finite = stacked[np.isfinite(stacked)]
    vmin, vmax = (0.0, 1.0) if finite.size == 0 else np.quantile(finite, [0.01, 0.995])
    if vmax <= vmin:
        vmax = vmin + 1.0

    ims = []
    for ax, res in zip(axes, results):
        ims.append(_plot_heatmap(ax, res, f"{res.config['pressure_GPa']:.3f} GPa", vmin=float(vmin), vmax=float(vmax)))
    axes[1].set_ylabel("")
    cbar = fig.colorbar(ims[-1], ax=axes, fraction=0.03, pad=0.02)
    cbar.set_label("Intensity (counts s$^{-1}$)")
    fig.suptitle("exp2 two-pressure heatmap comparison", y=1.02)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_all_conditions_heatmap(results: Sequence[PhotoactivityResult], outpath: Path) -> None:
    fig, axes = plt.subplots(1, len(results), figsize=(6.3 * len(results), 5.4), sharey=True, constrained_layout=True)
    if len(results) == 1:
        axes = [axes]
    stacked = np.hstack([r.counts_per_s.T.ravel() for r in results])
    finite = stacked[np.isfinite(stacked)]
    vmin, vmax = (0.0, 1.0) if finite.size == 0 else np.quantile(finite, [0.01, 0.995])
    if vmax <= vmin:
        vmax = vmin + 1.0

    ims = []
    for ax, res in zip(axes, results):
        title = f"{res.config['experiment']} | {res.config['pressure_GPa']:.3f} GPa"
        ims.append(_plot_heatmap(ax, res, title, vmin=float(vmin), vmax=float(vmax)))
    for ax in axes[1:]:
        ax.set_ylabel("")
    cbar = fig.colorbar(ims[-1], ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("Intensity (counts s$^{-1}$)")
    fig.suptitle("all-conditions heatmap comparison", y=1.02)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_exp2_metric_comparison(results: Sequence[PhotoactivityResult], outpath: Path) -> None:
    if len(results) != 2:
        raise ValueError("exp2 comparison expects exactly two datasets")
    r1, r2 = sorted(results, key=lambda r: r.config['pressure_GPa'])
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.6))

    ax = axes[0]
    for res, color, tag in [(r1, 'C0', f"{r1.config['pressure_GPa']:.2f} GPa"), (r2, 'C2', f"{r2.config['pressure_GPa']:.2f} GPa")]:
        ax.plot(res.wavenumber_cm1, res.early_mean_counts_per_s, lw=2.0, color=color, label=f"{tag} early")
        ax.plot(res.wavenumber_cm1, res.late_mean_counts_per_s, lw=2.0, ls='--', color=color, label=f"{tag} late")
        ax.axvspan(res.config['peak_window_cm1'][0], res.config['peak_window_cm1'][1], alpha=0.08, color='grey')
        ax.axvspan(res.config['tail_window_cm1'][0], res.config['tail_window_cm1'][1], alpha=0.05, color='grey')
    ax.set_title("exp2 early/late averaged spectra")
    ax.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Exposure-normalized intensity (counts s$^{-1}$)")
    ax.legend(frameon=True)
    _clean_axes(ax)

    ax = axes[1]
    for res, color, ls_peak, ls_bg in [(r1, 'C0', '-', '--'), (r2, 'C2', '-', '--')]:
        early_sl = frame_slice(res.config['frame_windows']['early'])
        peak_norm = _norm_to_early(res.peak_metric, early_sl)
        bg_norm = _norm_to_early(res.background_metric, early_sl)
        ax.plot(res.time_s, peak_norm, lw=2.2, ls=ls_peak, color=color, label=f"{res.config['pressure_GPa']:.2f} GPa peak")
        ax.plot(res.time_s, bg_norm, lw=2.2, ls=ls_bg, color=color, label=f"{res.config['pressure_GPa']:.2f} GPa background")
    ax.axhline(1.0, color='0.3', lw=1.0)
    ax.set_title("exp2 metric evolution")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Metric / early-window mean")
    ax.legend(frameon=True, ncol=2)
    _clean_axes(ax)

    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_all_conditions_metric_panels(results: Sequence[PhotoactivityResult], outpath: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16.2, 10.8))
    res_sorted = sorted(results, key=lambda r: (r.config['experiment'], r.config['pressure_GPa']))

    for res in res_sorted:
        label = f"{res.config['experiment']} {res.config['pressure_GPa']:.2f} GPa"
        early_sl = frame_slice(res.config['frame_windows']['early'])
        peak_norm = _norm_to_early(res.peak_metric, early_sl)
        bg_norm = _norm_to_early(res.background_metric, early_sl)
        axes[0, 0].plot(res.time_s, peak_norm, lw=2.2, label=label)
        axes[0, 1].plot(res.time_s, bg_norm, lw=2.2, label=label)
        axes[1, 0].plot(res.time_norm, peak_norm, lw=2.2, label=label)
        axes[1, 1].plot(res.time_norm, bg_norm, lw=2.2, label=label)

    titles = [
        "Peak metric vs real time",
        "Background metric vs real time",
        "Peak metric vs normalized run progress",
        "Background metric vs normalized run progress",
    ]
    xlabels = ["Time (s)", "Time (s)", "t / t$_{max}$", "t / t$_{max}$"]
    for ax, title, xlabel in zip(axes.ravel(), titles, xlabels):
        ax.axhline(1.0, color='0.3', lw=1.0)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Metric / early-window mean")
        _clean_axes(ax)
    axes[0, 0].legend(frameon=True, loc='best')
    axes[0, 1].legend(frameon=True, loc='best')
    axes[1, 0].legend(frameon=True, loc='best')
    axes[1, 1].legend(frameon=True, loc='best')
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_change_bars(results: Sequence[PhotoactivityResult], outpath: Path) -> None:
    labels = [f"{r.config['experiment']} {r.config['pressure_GPa']:.2f} GPa" for r in results]
    peak_change = [r.peak_change_pct for r in results]
    bg_change = [r.background_change_pct for r in results]

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.2))
    axes[0].bar(labels, peak_change)
    axes[0].set_title("Early→late peak change")
    axes[0].set_ylabel("Change (%)")
    axes[0].tick_params(axis='x', rotation=18)
    _clean_axes(axes[0])

    axes[1].bar(labels, bg_change)
    axes[1].set_title("Early→late background change")
    axes[1].set_ylabel("Change (%)")
    axes[1].tick_params(axis='x', rotation=18)
    _clean_axes(axes[1])

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)


def build_summary_table(results: Sequence[PhotoactivityResult]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for res in results:
        early_sl = frame_slice(res.config['frame_windows']['early'])
        late_sl = frame_slice(res.config['frame_windows']['late'])
        rows.append({
            'dataset': res.config['dataset_name'],
            'experiment': res.config['experiment'],
            'pressure_GPa': res.config['pressure_GPa'],
            'excitation_wavelength_nm': res.config['excitation_wavelength_nm'],
            'exposure_s_per_frame': res.config['exposure_s_per_frame'],
            'n_frames': len(res.frame_index),
            'total_time_s': float(res.time_s[-1]),
            'peak_window_nm': format_window(res.config['peak_window_nm']),
            'peak_window_cm1': format_window(res.config['peak_window_cm1']),
            'tail_window_nm': format_window(res.config['tail_window_nm']),
            'tail_window_cm1': format_window(res.config['tail_window_cm1']),
            'dominant_peak_center_nm': res.dominant_peak_center_nm,
            'dominant_peak_center_cm1': res.dominant_peak_center_cm1,
            'peak_metric_early': float(np.nanmean(res.peak_metric[early_sl])),
            'peak_metric_late': float(np.nanmean(res.peak_metric[late_sl])),
            'peak_change_pct': res.peak_change_pct,
            'background_metric_early': float(np.nanmean(res.background_metric[early_sl])),
            'background_metric_late': float(np.nanmean(res.background_metric[late_sl])),
            'background_change_pct': res.background_change_pct,
        })
    return pd.DataFrame(rows)


def save_summary_workbook(results: Sequence[PhotoactivityResult], outpath: Path) -> None:
    summary_df = build_summary_table(results)
    readme_rows = [
        ["Photoactivity reanalysis summary", ""],
        ["", ""],
        ["Method summary", ""],
        ["1. Input structure", "first column = detected wavelength (nm); columns 2-601 = 600 time frames."],
        ["2. Axis conversion", "Detected wavelength was converted to wavenumber using wavenumber = 1e7*(1/lambda0 - 1/lambda), with lambda0 set by each experiment laser line."],
        ["3. Exposure normalization", "Counts were divided by exposure time per frame before comparison across experiments."],
        ["4. Peak metric", "Integrated area in a fixed narrow wavenumber window around the dominant Raman peak."],
        ["5. Background metric", "Mean intensity in a higher-wavenumber tail region without a sharp peak."],
        ["6. Early/late comparison", "Early = frames 20-70; late = frames 530-580."],
        ["7. Plot convention", "Single-dataset outputs are heatmaps with x = wavenumber and y = time (frame index). Comparison figures are exported to a dedicated comparison_plots folder."],
        ["8. Interpretation note", "This workflow tests time-dependent spectral evolution under laser illumination, not photocatalysis."],
    ]

    with pd.ExcelWriter(outpath, engine="openpyxl") as writer:
        pd.DataFrame(readme_rows).to_excel(writer, sheet_name="README", header=False, index=False)
        summary_df.to_excel(writer, sheet_name="Summary", index=False)

        for res in results:
            peak_norm = np.nanmean(res.peak_metric[frame_slice(res.config['frame_windows']['early'])])
            bg_norm = np.nanmean(res.background_metric[frame_slice(res.config['frame_windows']['early'])])
            metric_df = pd.DataFrame({
                'frame': res.frame_index,
                'time_s': res.time_s,
                'time_norm': res.time_norm,
                'peak_metric_counts_per_s_cm1': res.peak_metric,
                'background_metric_counts_per_s': res.background_metric,
                'peak_metric_norm_to_early': res.peak_metric / peak_norm,
                'background_metric_norm_to_early': res.background_metric / bg_norm,
            })
            spectra_df = pd.DataFrame({
                'wavelength_nm': res.wavelength_nm,
                'wavenumber_cm1': res.wavenumber_cm1,
                'mean_counts_per_s': res.mean_spectrum_counts_per_s,
                'early_mean_counts_per_s': res.early_mean_counts_per_s,
                'late_mean_counts_per_s': res.late_mean_counts_per_s,
            })
            short_name = (
                f"{res.config['experiment']}_{res.config['pressure_GPa']:.2f}".replace('.', '_')
                .replace('-', 'm')
            )
            metric_df.to_excel(writer, sheet_name=f"{short_name[:23]}_metrics", index=False)
            spectra_df.to_excel(writer, sheet_name=f"{short_name[:23]}_spectra", index=False)


def write_run_log(output_dir: Path, run_profile: Dict[str, Any], results: Sequence[PhotoactivityResult]) -> None:
    payload = {
        'run_profile': run_profile['name'],
        'description': run_profile.get('description', ''),
        'output_dir': str(output_dir.resolve()),
        'datasets': [
            {
                'dataset_name': res.config['dataset_name'],
                'excel_path': str(Path(res.config['excel_path']).resolve()),
                'pressure_GPa': res.config['pressure_GPa'],
                'excitation_wavelength_nm': res.config['excitation_wavelength_nm'],
                'exposure_s_per_frame': res.config['exposure_s_per_frame'],
                'peak_window_nm': res.config['peak_window_nm'],
                'peak_window_cm1': res.config['peak_window_cm1'],
                'tail_window_nm': res.config['tail_window_nm'],
                'tail_window_cm1': res.config['tail_window_cm1'],
                'peak_change_pct': res.peak_change_pct,
                'background_change_pct': res.background_change_pct,
            }
            for res in results
        ],
    }
    with open(output_dir / 'run_log.json', 'w', encoding='utf-8') as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)


def run_photoactivity_profile(run_profile: Dict[str, Any], dataset_configs: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    outdir = Path(run_profile['output_dir'])
    ensure_dir(outdir)

    results: List[PhotoactivityResult] = []
    for cfg in dataset_configs:
        excel_path = Path(cfg['excel_path'])
        if not excel_path.exists():
            raise FileNotFoundError(f"Input file not found: {excel_path}")
        res = analyze_single_dataset(cfg)
        results.append(res)
        plot_single_heatmap(res, outdir / f"{cfg['dataset_name']}_heatmap.png")

    comparison_dir = outdir / 'comparison_plots'
    if len(results) >= 2:
        ensure_dir(comparison_dir)
        exp2_results = [r for r in results if r.config['experiment'] == 'exp2']
        if len(exp2_results) == 2:
            plot_exp2_two_pressure_heatmap(exp2_results, comparison_dir / 'exp2_two_pressure_heatmap_comparison.png')
            plot_exp2_metric_comparison(exp2_results, comparison_dir / 'exp2_two_pressure_metric_comparison.png')

        plot_all_conditions_heatmap(results, comparison_dir / 'all_conditions_heatmap_comparison.png')
        plot_all_conditions_metric_panels(results, comparison_dir / 'all_conditions_metric_panels.png')
        plot_change_bars(results, comparison_dir / 'all_conditions_early_late_change_bars.png')

    summary_xlsx = outdir / 'photoactivity_reanalysis_summary.xlsx'
    summary_csv = outdir / 'summary_metrics.csv'
    save_summary_workbook(results, summary_xlsx)
    build_summary_table(results).to_csv(summary_csv, index=False)
    write_run_log(outdir, run_profile, results)

    return {
        'output_dir': str(outdir.resolve()),
        'n_datasets': len(results),
        'summary_xlsx': str(summary_xlsx.resolve()),
        'summary_csv': str(summary_csv.resolve()),
    }
