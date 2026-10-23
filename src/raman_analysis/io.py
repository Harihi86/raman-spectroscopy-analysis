"""Explicit tabular imports: no layout guessing, truncation or duplicate averaging."""

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Spectrum:
    x: np.ndarray
    y: np.ndarray
    sigma: np.ndarray | None
    provenance: dict


def load_spectrum(
    path: Path,
    *,
    x="shift_cm1",
    y="intensity",
    sigma=None,
    sheet=0,
    delimiter=None,
    missing="error",
    units="cm^-1",
):
    path = Path(path)
    if units != "cm^-1":
        raise ValueError(
            "The supported spectral-axis unit is cm^-1; convert wavelengths explicitly"
        )
    if missing not in {"error", "drop"}:
        raise ValueError("missing must be error or drop")
    if path.suffix.lower() == ".xlsx":
        frame = pd.read_excel(path, sheet_name=sheet)
    elif path.suffix.lower() in {".csv", ".tsv", ".txt"}:
        sep = (
            delimiter if delimiter is not None else ("\t" if path.suffix.lower() == ".tsv" else ",")
        )
        frame = pd.read_csv(path, sep=sep)
    else:
        raise ValueError("Supported formats: CSV, TSV, delimited TXT and XLSX, with a header row")

    def column(key):
        values = frame.iloc[:, key] if type(key) is int else frame[key]
        return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)

    axis, intensity = column(x), column(y)
    errors = column(sigma) if sigma is not None else None
    valid = np.isfinite(axis) & np.isfinite(intensity)
    if errors is not None:
        valid &= np.isfinite(errors) & (errors > 0)
    dropped = int(np.count_nonzero(~valid))
    if dropped and missing == "error":
        raise ValueError(
            f"{dropped} rows have missing, nonnumeric or invalid values; choose missing=drop explicitly"
        )
    axis, intensity = axis[valid], intensity[valid]
    if errors is not None:
        errors = errors[valid]
    if len(axis) < 8:
        raise ValueError("At least eight valid spectral points are required")
    if len(np.unique(axis)) != len(axis):
        raise ValueError("Duplicate axis values require an explicit aggregation decision")
    order = np.argsort(axis)
    sorted_input = bool(np.array_equal(order, np.arange(len(order))))
    axis, intensity = axis[order].copy(), intensity[order].copy()
    if errors is not None:
        errors = errors[order].copy()
    for array in (axis, intensity, errors):
        if array is not None:
            array.setflags(write=False)
    return Spectrum(
        axis,
        intensity,
        errors,
        {
            "input_name": path.name,
            "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "units": units,
            "x_column": x,
            "y_column": y,
            "sigma_column": sigma,
            "sheet": sheet if path.suffix.lower() == ".xlsx" else None,
            "missing_policy": missing,
            "dropped_rows": dropped,
            "input_already_sorted": sorted_input,
            "valid_points": len(axis),
            "duplicate_policy": "error",
        },
    )
