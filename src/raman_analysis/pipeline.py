import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import __version__
from .automatic import auto_fit
from .baseline import estimate_baseline
from .fitting import fit_peaks
from .io import load_spectrum
from .report import selection_html


def analyze(config, output: Path, *, base: Path = Path("."), save_figures: bool = True):
    if set(config) - {"input", "baseline", "fit", "data_kind", "metadata"}:
        raise ValueError("Unknown top-level configuration field")
    data_kind = config.get("data_kind", "unspecified")
    if data_kind not in {"synthetic", "measurement", "unspecified"}:
        raise ValueError("data_kind must be synthetic, measurement or unspecified")
    source = dict(config["input"])
    path = (base / source.pop("path")).resolve()
    spectrum = load_spectrum(path, **source)
    baseline_config = config.get("baseline", {"method": "none"})
    baseline = estimate_baseline(spectrum.x, spectrum.y, baseline_config)
    corrected = spectrum.y - baseline
    fit_config = dict(config["fit"])
    window = fit_config.pop("window", [float(spectrum.x[0]), float(spectrum.x[-1])])
    if len(window) != 2 or not np.isfinite(window).all() or window[0] >= window[1]:
        raise ValueError("Fit window must be a finite increasing interval")
    selection = (spectrum.x >= window[0]) & (spectrum.x <= window[1])
    if ("peaks" in fit_config) == ("auto" in fit_config):
        raise ValueError("Provide exactly one of fit.peaks or fit.auto")
    fit_function = auto_fit if "auto" in fit_config else fit_peaks
    fit = fit_function(
        spectrum.x[selection],
        corrected[selection],
        errors=spectrum.sigma[selection] if spectrum.sigma is not None else None,
        **fit_config,
    )
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    summary = {
        k: v
        for k, v in fit.items()
        if k not in {"fitted", "residuals", "components", "residual_background"}
    }
    summary.update(
        {
            "data_kind": data_kind,
            "program_version": __version__,
            "provenance": spectrum.provenance,
            "baseline": baseline_config,
            "fit_configuration": config["fit"],
            "metadata": config.get("metadata", {}),
            "python": platform.python_version(),
            "dependencies": {
                name: importlib.metadata.version(name)
                for name in ("numpy", "scipy", "matplotlib", "pandas", "openpyxl")
            },
        }
    )
    (output / "report.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    (output / "configuration.json").write_text(
        json.dumps(config, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    pd.DataFrame(
        fit["peaks"],
        columns=list(fit["peaks"][0])
        if fit["peaks"]
        else ["center_cm1", "fwhm_cm1", "amplitude", "area"],
    ).to_csv(output / "peaks.csv", index=False)
    if "model_selection" in fit:
        model_selection = fit["model_selection"]
        (output / "model-selection.html").write_text(
            selection_html(model_selection, fit["status"], fit["diagnostics"]), encoding="utf-8"
        )
        (output / "model-selection.json").write_text(
            json.dumps(model_selection, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        pd.DataFrame(
            [{k: v for k, v in c.items() if k != "attempts"} for c in model_selection["candidates"]]
        ).to_csv(output / "model-selection.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(
        {
            "shift_cm1": spectrum.x,
            "raw_intensity": spectrum.y,
            "estimated_baseline": baseline,
            "corrected_intensity": corrected,
        }
    ).to_csv(output / "processed.csv", index=False)
    fit_frame = {
        "shift_cm1": spectrum.x[selection],
        "fitted_corrected": fit["fitted"],
        "residual": fit["residuals"],
        "residual_background": fit["residual_background"],
    }
    if spectrum.sigma is not None:
        fit_frame["measurement_sigma"] = spectrum.sigma[selection]
    for i, component in enumerate(fit["components"]):
        fit_frame[f"gaussian_{i + 1}"] = component
    pd.DataFrame(fit_frame).to_csv(output / "fit.csv", index=False)
    if save_figures:
        plot_spectrum(spectrum, baseline, corrected, selection, fit, data_kind, output)
    manifest = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(output.iterdir())
        if p.is_file()
    }
    (output / "artifacts.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return summary


def plot_spectrum(spectrum, baseline, corrected, selection, fit, data_kind, output):
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(8.0, 8.0),
        sharex=True,
        layout="constrained",
        gridspec_kw={"height_ratios": [1, 1.3, 0.7]},
    )
    axes[0].plot(spectrum.x, spectrum.y, color="#435669", lw=0.9, label="Raw input")
    axes[0].plot(spectrum.x, baseline, color="#d47835", label="Estimated baseline")
    axes[1].plot(
        spectrum.x[selection],
        corrected[selection],
        color="#8294a1",
        lw=0.8,
        label="Baseline corrected",
    )
    axes[1].plot(
        spectrum.x[selection],
        fit["fitted"],
        color="#137b80",
        label="Gaussian model + linear residual background",
    )
    for i, curve in enumerate(fit["components"]):
        axes[1].plot(spectrum.x[selection], curve, "--", lw=1, label=f"Gaussian {i + 1}")
    axes[2].plot(spectrum.x[selection], fit["residuals"], color="#435669", lw=0.7)
    axes[2].axhline(0, color="#444444", lw=0.6)
    for ax in axes[:2]:
        ax.legend(fontsize=8, frameon=False)
        ax.set_ylabel("Intensity (a.u.)")
    axes[2].set_ylabel("Residual (a.u.)")
    axes[2].set_xlabel(r"Raman shift (cm$^{-1}$)")
    decision = "AUTO" if "model_selection" in fit else "MANUAL"
    fig.suptitle(
        f"{data_kind.upper()} | {decision} | {len(fit['peaks'])} Gaussian peaks | {fit['status']}",
        fontsize=12,
    )
    fig.savefig(output / "spectrum.png", dpi=180)
    fig.savefig(output / "spectrum.pdf")
    plt.close(fig)
