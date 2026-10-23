"""Generated, explicitly synthetic examples. No experimental profile is imported."""

from pathlib import Path

import numpy as np
import pandas as pd

from .batch import analyze_batch, write_json
from .fitting import gaussian


def demo_batch(output: Path, *, scenario="basic", automatic=False, no_figures=False, progress=None):
    if scenario not in {"basic", "crossing"}:
        raise ValueError("scenario must be basic or crossing")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    (output / "data").mkdir()
    x = np.linspace(100, 900, 641)
    rng = np.random.default_rng(20261023)
    if scenario == "basic":
        forward = np.linspace(0, 6, 6)
        series = {"sweep-a": np.insert(forward, 3, forward[2]), "sweep-b": np.linspace(6, 0, 4)}
    else:
        series = {"crossing": np.array([0.0, 1.0, 2.0, 4.0, 5.0, 6.0])}
    spectra, truth = [], []
    for name, pressures in series.items():
        for index, pressure in enumerate(pressures):
            centers = (
                [300 + 5 * pressure, 650 - 3 * pressure]
                if scenario == "basic"
                else [260 + 60 * pressure, 620 - 60 * pressure]
            )
            widths = [18, 26] if scenario == "basic" else [12, 20]
            amplitudes = [22, 15]
            y = (
                4
                + x * 0.002
                + sum(
                    gaussian(x, a, c, s)
                    for a, c, s in zip(amplitudes, centers, widths, strict=True)
                )
                + rng.normal(0, 0.08, len(x))
            )
            spectrum_id = f"{name}-{index + 1:03d}"
            filename = f"data/{spectrum_id}.csv"
            pd.DataFrame(
                {"shift_cm1": x, "intensity": y, "noise_sd": np.full(len(x), 0.08)}
            ).to_csv(output / filename, index=False)
            spectra.append(
                {
                    "id": spectrum_id,
                    "series": name,
                    "pressure_gpa": float(pressure),
                    "input": {"path": filename},
                    "analysis": {
                        "fit": {
                            "peaks": [
                                {"center": float(c + 2), "sigma": float(s * 0.9)}
                                for c, s in zip(centers, widths, strict=True)
                            ]
                        }
                    },
                }
            )
            truth.append(
                {
                    "spectrum_id": spectrum_id,
                    "series": name,
                    "pressure_gpa": float(pressure),
                    "centers_in_physical_branch_order_cm1": [float(c) for c in centers],
                }
            )
    manifest = {
        "schema_version": 1,
        "data_kind": "synthetic",
        "defaults": {"input": {"sigma": "noise_sd"}, "baseline": {"method": "none"}},
        "outputs": {"per_spectrum_figures": False, "tracking_figures": True},
        "tracking": {
            "prediction": "linear",
            "max_distance_cm1": 80 if scenario == "crossing" else 20,
            "ambiguity_margin_cm1": 1,
            "max_gap": 1,
        },
        "metadata": {
            "seed": 20261023,
            "scenario": scenario,
            "purpose": "synthetic software validation; not experimental data",
            "fit_initialization": "explicit demonstration priors near generated peaks; not blind peak discovery",
        },
        "spectra": spectra,
    }
    if automatic:
        manifest["defaults"]["fit"] = {"auto": {"max_peaks": 4}}
        manifest["metadata"]["fit_initialization"] = (
            "automatic proposals; generated peak parameters are not supplied to the fitter"
        )
        for spectrum in spectra:
            spectrum.pop("analysis")
    write_json(output / "manifest.json", manifest)
    write_json(output / "truth.json", {"data_kind": "synthetic", "spectra": truth})
    return analyze_batch(
        output / "manifest.json", output / "analysis", no_figures=no_figures, progress=progress
    )
