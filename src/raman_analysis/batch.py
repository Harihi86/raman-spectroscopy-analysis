"""Manifest-driven batch processing and portable, explicit tracking artifacts."""

import hashlib
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import __version__
from .manifest import read_schema, validate_json, validate_manifest
from .pipeline import analyze
from .report import batch_html
from .tracking import track_peaks


def write_json(path: Path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def write_artifacts(directory: Path):
    files = {
        p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(directory.rglob("*"))
        if p.is_file() and p != directory / "artifacts.json"
    }
    write_json(directory / "artifacts.json", files)


def export_tracking(result: dict, output: Path, *, figures=True):
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "tracking.json", result)
    columns = [
        "spectrum_id",
        "label",
        "data_kind",
        "series",
        "sequence",
        "pressure_gpa",
        "fit_status",
        "peak_index",
        "track_id",
        "event",
        "candidate_track_ids",
        "prediction_cm1",
        "center_cm1",
        "center_se_cm1",
        "fwhm_cm1",
        "amplitude",
        "area",
    ]
    rows = [
        {**r, "candidate_track_ids": json.dumps(r["candidate_track_ids"], ensure_ascii=False)}
        for r in result["assignments"]
    ]
    pd.DataFrame(rows, columns=columns).to_csv(
        output / "assignments.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(
        result["tracks"],
        columns=[
            "track_id",
            "series",
            "observations",
            "first_spectrum_id",
            "last_spectrum_id",
            "closed_by_gap_limit",
        ],
    ).to_csv(output / "tracks.csv", index=False, encoding="utf-8-sig")
    figures_created = []
    if figures:
        series_names = list(dict.fromkeys(r["series"] for r in result["assignments"]))
        for index, series in enumerate(series_names, start=1):
            frame = [r for r in result["assignments"] if r["series"] == series]
            fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
            ids = list(dict.fromkeys(r["track_id"] for r in frame if r["track_id"]))
            for track_id in ids:
                points = [r for r in frame if r["track_id"] == track_id]
                # GAP rows contain no measured center; NaNs break the line explicitly.
                x = [p["pressure_gpa"] if p["pressure_gpa"] is not None else np.nan for p in points]
                y = [
                    p["center_cm1"] if p["event"] in {"NEW", "MATCHED"} else np.nan for p in points
                ]
                error = [
                    p["center_se_cm1"] if p["center_se_cm1"] is not None else np.nan for p in points
                ]
                ax.errorbar(x, y, yerr=error, marker="o", ms=4, lw=1.2, capsize=2, label=track_id)
            excluded = [
                p
                for p in frame
                if p["track_id"] is None
                and p["center_cm1"] is not None
                and p["pressure_gpa"] is not None
            ]
            if excluded:
                ax.scatter(
                    [p["pressure_gpa"] for p in excluded],
                    [p["center_cm1"] for p in excluded],
                    marker="x",
                    color="#777777",
                    label="Ambiguous / excluded fit",
                )
            kinds = set(p["data_kind"] for p in frame)
            kind = next(iter(kinds)).upper() if len(kinds) == 1 else "MIXED DATA"
            ax.set(
                xlabel="Pressure (GPa)",
                ylabel=r"Peak center (cm$^{-1}$)",
                title=f"{kind} | {series} | {result['status']} | acquisition order",
            )
            ax.grid(alpha=0.2)
            if ids or excluded:
                ax.legend(fontsize=8, frameon=False)
            stem = f"series-{index:03d}"
            for extension in ("png", "pdf"):
                fig.savefig(output / f"{stem}.{extension}", dpi=180)
            plt.close(fig)
            figures_created.append({"series": series, "png": f"{stem}.png", "pdf": f"{stem}.pdf"})
    write_json(output / "figures.json", figures_created)
    write_artifacts(output)


def _batch_status(records, tracking):
    statuses = [r["status"] for r in records]
    if all(s == "PASS" for s in statuses) and tracking["status"] in {"PASS", "NOT_RUN"}:
        return "PASS"
    if any(s == "FAIL" for s in statuses):
        return "FAIL"
    if all(s == "BLOCKED" for s in statuses):
        return "BLOCKED"
    return "INCONCLUSIVE"


def analyze_batch(manifest_path: Path, output: Path, *, no_figures=False, progress=None) -> dict:
    manifest_path, output = Path(manifest_path), Path(output)
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes.decode("utf-8-sig"))
    normalized = validate_manifest(manifest)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "manifest.json", manifest)
    write_json(output / "resolved-spectra.json", normalized)
    state = {
        "format": "raman-analysis-batch",
        "schema_version": 1,
        "program_version": __version__,
        "state": "RUNNING",
        "status": "NOT_RUN",
        "manifest_name": manifest_path.name,
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "completed": 0,
        "total": len(normalized),
        "output": output.name,
    }
    write_json(output / "run.json", state)
    records = []
    per_spectrum_figures = manifest.get("outputs", {}).get("per_spectrum_figures", True)
    try:
        with (output / "spectra.jsonl").open("w", encoding="utf-8") as journal:
            for entry in normalized:
                directory = output / "spectra" / entry["id"]
                record = {
                    k: entry[k] for k in ("id", "label", "series", "sequence", "pressure_gpa")
                }
                record["data_kind"] = entry["analysis"]["data_kind"]
                record["output"] = directory.relative_to(output).as_posix()
                try:
                    report = analyze(
                        entry["analysis"],
                        directory,
                        base=manifest_path.resolve().parent,
                        save_figures=per_spectrum_figures and not no_figures,
                    )
                    record.update(
                        status=report["status"],
                        peaks=report["peaks"],
                        diagnostics=report["diagnostics"],
                        provenance=report["provenance"],
                    )
                    if "model_selection" in report:
                        decision = report["model_selection"]
                        record["model_selection"] = {
                            k: decision[k]
                            for k in (
                                "criterion",
                                "selected_peak_count",
                                "plausible_peak_counts",
                                "count_status",
                                "search_complete",
                                "selected_shape",
                                "shape_count_stability",
                            )
                        }
                except (ValueError, KeyError, TypeError, OSError, IndexError) as exc:
                    # One bad input is retained as BLOCKED, without discarding other spectra.
                    record.update(
                        status="BLOCKED",
                        peaks=[],
                        diagnostics=[str(exc)],
                        error_type=type(exc).__name__,
                    )
                    directory.mkdir(parents=True, exist_ok=True)
                    write_json(directory / "report.json", record)
                records.append(record)
                journal.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
                journal.flush()
                state["completed"] = len(records)
                write_json(output / "run.json", state)
                if progress is not None:
                    progress(len(records), len(normalized), record["id"], record["status"])
        tracking = track_peaks(records, manifest.get("tracking", {}))
        export_tracking(
            tracking,
            output / "tracking",
            figures=(manifest.get("outputs", {}).get("tracking_figures", True) and not no_figures),
        )
        columns = [
            "id",
            "label",
            "data_kind",
            "series",
            "sequence",
            "pressure_gpa",
            "status",
            "peak_count",
            "count_status",
            "plausible_peak_counts",
            "selected_shape",
            "shape_count_stability",
            "output",
        ]
        pd.DataFrame(
            [
                {
                    **r,
                    "peak_count": len(r["peaks"]),
                    "count_status": r.get("model_selection", {}).get("count_status", "MANUAL"),
                    "plausible_peak_counts": json.dumps(
                        r.get("model_selection", {}).get("plausible_peak_counts", [])
                    ),
                    "selected_shape": r.get("model_selection", {}).get("selected_shape", "MANUAL"),
                    "shape_count_stability": r.get("model_selection", {}).get(
                        "shape_count_stability", "NOT_ASSESSED"
                    ),
                }
                for r in records
            ],
            columns=columns,
        ).to_csv(output / "summary.csv", index=False, encoding="utf-8-sig")
        state.update(
            state="COMPLETE",
            status=_batch_status(records, tracking),
            fit_counts=dict(Counter(r["status"] for r in records)),
            tracking_status=tracking["status"],
            tracks=len(tracking["tracks"]),
        )
        report = {
            **state,
            "spectra": records,
            "tracking_settings": tracking["settings"],
            "tracking_report": "tracking/tracking.json",
            "metadata": manifest.get("metadata", {}),
        }
        (output / "report.html").write_text(batch_html(report, tracking, output), encoding="utf-8")
        write_json(output / "batch-report.json", report)
        write_json(output / "run.json", state)
        write_artifacts(output)
        return report
    except BaseException as exc:
        # A completed report is a reusable input to `track`; never leave one
        # behind if final export or artifact hashing failed.
        (output / "batch-report.json").unlink(missing_ok=True)
        state.update(
            state="INTERRUPTED" if isinstance(exc, KeyboardInterrupt) else "ABORTED",
            status="INCONCLUSIVE",
            error_type=type(exc).__name__,
        )
        write_json(output / "run.json", state)
        raise


def retrack(batch_report: Path, output: Path, *, settings_path: Path | None = None, figures=True):
    raw = Path(batch_report).read_bytes()
    report = json.loads(raw.decode("utf-8-sig"))
    if (
        report.get("format") != "raman-analysis-batch"
        or report.get("schema_version") != 1
        or report.get("state") != "COMPLETE"
    ):
        raise ValueError("Tracking requires a complete version-1 batch-report.json")
    settings = report.get("tracking_settings", {})
    if settings_path is not None:
        settings = json.loads(Path(settings_path).read_text(encoding="utf-8-sig"))
    validate_json(settings, read_schema("tracking-settings.schema.json"))
    result = track_peaks(report["spectra"], settings)
    result["source_batch_sha256"] = hashlib.sha256(raw).hexdigest()
    export_tracking(result, Path(output), figures=figures)
    return result
