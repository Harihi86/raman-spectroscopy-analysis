"""Pressure-aware assignment with explicit gaps, births and global ambiguity checks."""

import math
from collections import Counter
from dataclasses import asdict, dataclass, field

import numpy as np
from scipy.optimize import linear_sum_assignment


@dataclass(frozen=True)
class TrackingSettings:
    enabled: bool = True
    prediction: str = "linear"
    max_distance_cm1: float = 20.0
    ambiguity_margin_cm1: float = 1.0
    max_gap: int = 1
    width_weight: float = 0.0

    def __post_init__(self):
        if type(self.enabled) is not bool or self.prediction not in {"linear", "nearest"}:
            raise ValueError("Invalid tracking enabled/prediction setting")
        for name in ("max_distance_cm1", "ambiguity_margin_cm1", "width_weight"):
            value = getattr(self, name)
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError(f"tracking.{name} must be a nonnegative finite number")
        if self.max_distance_cm1 == 0:
            raise ValueError("tracking.max_distance_cm1 must be positive")
        if type(self.max_gap) is not int or not 0 <= self.max_gap <= 50:
            raise ValueError("tracking.max_gap must be an integer in [0, 50]")


@dataclass
class _Track:
    id: str
    series: str
    points: list[dict] = field(default_factory=list)
    last_step: int = -1
    closed: bool = False

    def predict(self, pressure: float, settings: TrackingSettings) -> float:
        last = self.points[-1]
        center = last["center_cm1"]
        if settings.prediction == "linear":
            # Replicates at one pressure do not introduce a division by zero.
            earlier = next(
                (
                    p
                    for p in reversed(self.points[:-1])
                    if not math.isclose(
                        p["pressure_gpa"], last["pressure_gpa"], rel_tol=1e-12, abs_tol=1e-12
                    )
                ),
                None,
            )
            if earlier is not None:
                slope = (center - earlier["center_cm1"]) / (
                    last["pressure_gpa"] - earlier["pressure_gpa"]
                )
                center += slope * (pressure - last["pressure_gpa"])
        if not math.isfinite(center):
            raise ValueError(
                "Linear prediction overflow: inspect pressure values or use nearest prediction"
            )
        return float(center)


def _finite(value, name, *, nullable=False, positive=False):
    if nullable and value is None:
        return
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")


def _validate_records(records, settings):
    ids, sequences = set(), set()
    for record in records:
        if not isinstance(record.get("id"), str) or not record["id"]:
            raise ValueError("Every tracking record needs an id")
        if record["id"] in ids:
            raise ValueError("Tracking spectrum ids must be unique")
        ids.add(record["id"])
        sequence = record.get("sequence")
        if type(sequence) is not int or sequence < 0 or sequence in sequences:
            raise ValueError("Tracking sequence numbers must be unique nonnegative integers")
        sequences.add(sequence)
        if not isinstance(record.get("series"), str) or not record["series"]:
            raise ValueError("Every tracking record needs an explicit series")
        if record.get("status") not in {"PASS", "FAIL", "INCONCLUSIVE", "BLOCKED"}:
            raise ValueError("Unknown fit status in tracking input")
        _finite(record.get("pressure_gpa"), "pressure_gpa", nullable=True)
        peaks = record.get("peaks", [])
        if not isinstance(peaks, list) or len(peaks) > 8:
            raise ValueError("Each tracking frame must contain at most eight peak observations")
        for peak in peaks:
            _finite(peak.get("center_cm1"), "center_cm1")
            if peak.get("center_se_cm1") is not None:
                _finite(peak["center_se_cm1"], "center_se_cm1")
                if peak["center_se_cm1"] < 0:
                    raise ValueError("center_se_cm1 cannot be negative")
            if settings.width_weight > 0 or peak.get("fwhm_cm1") is not None:
                _finite(peak.get("fwhm_cm1"), "fwhm_cm1", positive=True)


def _possible_assignments(cost, n_tracks, n_peaks, margin):
    """Enumerate globally competitive identities, including an unmatched observation.

    For each observation, force each feasible real edge or its birth edge and
    solve the remaining global assignment. This avoids row-wise nearest-neighbor
    confidence claims when observations compete for the same track.
    """
    rows, columns = linear_sum_assignment(cost)
    optimum = float(cost[rows, columns].sum())
    tolerance = 1e-9 * max(1.0, abs(optimum))
    possible = []
    for j in range(n_peaks):
        choices = []
        for i in [*range(n_tracks), n_tracks + j]:
            if not np.isfinite(cost[i, j]):
                continue
            forced = cost.copy()
            forced[i, :] = np.inf
            forced[:, j] = np.inf
            forced[i, j] = cost[i, j]
            try:
                rr, cc = linear_sum_assignment(forced)
            except ValueError:
                continue
            delta = float(forced[rr, cc].sum()) - optimum
            if delta <= margin + tolerance:
                choices.append((i if i < n_tracks else None, max(delta, 0.0)))
        possible.append(choices)
    return possible


def _row(
    record, *, peak_index=None, peak=None, track_id=None, event, candidates=None, prediction=None
):
    peak = peak or {}
    return {
        "spectrum_id": record["id"],
        "label": record.get("label", record["id"]),
        "data_kind": record.get("data_kind", "unspecified"),
        "series": record["series"],
        "sequence": record["sequence"],
        "pressure_gpa": record.get("pressure_gpa"),
        "fit_status": record["status"],
        "peak_index": peak_index,
        "track_id": track_id,
        "event": event,
        "candidate_track_ids": candidates or [],
        "prediction_cm1": prediction,
        "center_cm1": peak.get("center_cm1"),
        "center_se_cm1": peak.get("center_se_cm1"),
        "fwhm_cm1": peak.get("fwhm_cm1"),
        "amplitude": peak.get("amplitude"),
        "area": peak.get("area"),
    }


def track_peaks(records: list[dict], settings: dict | None = None) -> dict:
    """Track accepted observations in acquisition order within independent series.

    IDs denote algorithmic trajectories, not chemical assignments. Failed or
    inconclusive fits never update a trajectory. No missing peaks are invented.
    """
    config = TrackingSettings(**(settings or {}))
    _validate_records(records, config)
    if not config.enabled:
        return {
            "format": "raman-peak-tracking",
            "schema_version": 1,
            "status": "NOT_RUN",
            "settings": asdict(config),
            "assignments": [],
            "tracks": [],
            "diagnostics": [],
            "event_counts": {},
        }
    groups = {}
    for record in sorted(records, key=lambda r: r["sequence"]):
        groups.setdefault(record["series"], []).append(record)
    assignments, summaries, diagnostics = [], [], []
    for series, frames in groups.items():
        tracks = []
        for step, record in enumerate(frames):
            active = []
            for track in tracks:
                if step - track.last_step - 1 > config.max_gap:
                    track.closed = True
                if not track.closed:
                    active.append(track)
            pressure, peaks = record.get("pressure_gpa"), record.get("peaks", [])
            eligible = record["status"] == "PASS" and pressure is not None
            updated = set()
            if not eligible:
                reason = "EXCLUDED_FIT" if record["status"] != "PASS" else "MISSING_PRESSURE"
                diagnostics.append({"spectrum_id": record["id"], "code": reason})
                for j, peak in enumerate(peaks):
                    assignments.append(_row(record, peak_index=j, peak=peak, event=reason))
                if not peaks:
                    assignments.append(_row(record, event=reason))
            elif peaks:
                predictions = [t.predict(pressure, config) for t in active]
                n, m = len(active), len(peaks)
                if n:
                    cost = np.full((n + m, n + m), np.inf)
                    for i, track in enumerate(active):
                        for j, peak in enumerate(peaks):
                            distance = abs(predictions[i] - peak["center_cm1"])
                            if distance <= config.max_distance_cm1:
                                width_cost = 0.0
                                if config.width_weight:
                                    width_cost = config.width_weight * abs(
                                        track.points[-1]["fwhm_cm1"] - peak["fwhm_cm1"]
                                    )
                                cost[i, j] = distance + width_cost
                        cost[i, m + i] = config.max_distance_cm1 / 2
                    for j in range(m):
                        cost[n + j, j] = config.max_distance_cm1 / 2
                    cost[n:, m:] = 0
                    choices = _possible_assignments(cost, n, m, config.ambiguity_margin_cm1)
                else:
                    choices = [[(None, 0.0)] for _ in peaks]
                for j, (peak, options) in enumerate(zip(peaks, choices, strict=True)):
                    if len(options) != 1:
                        candidates = [active[i].id if i is not None else "NEW" for i, _ in options]
                        assignments.append(
                            _row(
                                record,
                                peak_index=j,
                                peak=peak,
                                event="AMBIGUOUS",
                                candidates=candidates,
                            )
                        )
                        diagnostics.append(
                            {
                                "spectrum_id": record["id"],
                                "peak_index": j,
                                "code": "AMBIGUOUS",
                                "candidates": candidates,
                            }
                        )
                        continue
                    index = options[0][0]
                    if index is None:
                        track = _Track(f"{series}-T{len(tracks) + 1:03d}", series)
                        tracks.append(track)
                        event, prediction = "NEW", None
                    else:
                        track = active[index]
                        event, prediction = "MATCHED", predictions[index]
                    point = _row(
                        record,
                        peak_index=j,
                        peak=peak,
                        track_id=track.id,
                        event=event,
                        prediction=prediction,
                    )
                    track.points.append(point)
                    track.last_step = step
                    updated.add(track.id)
                    assignments.append(point)
            for track in active:
                if track.id not in updated:
                    prediction = track.predict(pressure, config) if pressure is not None else None
                    assignments.append(
                        _row(record, track_id=track.id, event="GAP", prediction=prediction)
                    )
        for track in tracks:
            summaries.append(
                {
                    "track_id": track.id,
                    "series": series,
                    "observations": len(track.points),
                    "first_spectrum_id": track.points[0]["spectrum_id"],
                    "last_spectrum_id": track.points[-1]["spectrum_id"],
                    "closed_by_gap_limit": track.closed,
                }
            )
    assignments.sort(
        key=lambda r: (
            r["sequence"],
            r["peak_index"] is None,
            r["peak_index"] if r["peak_index"] is not None else 0,
            r["track_id"] or "",
        )
    )
    if not summaries:
        diagnostics.append(
            {"code": "NO_TRACKS", "message": "No accepted pressure-tagged observations"}
        )
    return {
        "format": "raman-peak-tracking",
        "schema_version": 1,
        "status": "INCONCLUSIVE" if diagnostics else "PASS",
        "settings": asdict(config),
        "assignments": assignments,
        "tracks": summaries,
        "diagnostics": diagnostics,
        "event_counts": dict(Counter(row["event"] for row in assignments)),
        "identity_scope": "algorithmic continuity labels; not chemical identification",
        "ambiguity_scope": "global assignment-cost margin heuristic; not a probability",
    }
