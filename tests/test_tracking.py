from itertools import permutations

import numpy as np
import pytest

from raman_analysis.tracking import _possible_assignments, track_peaks


def frame(index, pressure, centers, *, series="sample", status="PASS"):
    return {
        "id": f"s{index}",
        "sequence": index,
        "series": series,
        "pressure_gpa": pressure,
        "status": status,
        "data_kind": "synthetic",
        "peaks": [{"center_cm1": c, "center_se_cm1": 0.1, "fwhm_cm1": 10} for c in sorted(centers)],
    }


def observations(result, spectrum_id):
    return [
        r
        for r in result["assignments"]
        if r["spectrum_id"] == spectrum_id and r["peak_index"] is not None
    ]


def test_linear_prediction_preserves_identity_when_peak_order_crosses():
    pressures = [0, 1, 3, 4]
    records = [frame(i, p, [100 + 20 * p, 180 - 20 * p]) for i, p in enumerate(pressures)]
    result = track_peaks(records, {"max_distance_cm1": 25, "ambiguity_margin_cm1": 0.5})
    assert result["status"] == "PASS"
    final = {r["track_id"]: r["center_cm1"] for r in observations(result, "s3")}
    assert final == {"sample-T001": 180, "sample-T002": 100}
    assert observations(result, "s2")[0]["track_id"] == "sample-T002"


def test_exact_crossing_is_ambiguous_and_does_not_update_history():
    records = [frame(i, p, [100 + 20 * p, 180 - 20 * p]) for i, p in enumerate([0, 1, 2, 3])]
    result = track_peaks(records, {"max_distance_cm1": 25, "max_gap": 1})
    crossing = observations(result, "s2")
    assert result["status"] == "INCONCLUSIVE"
    assert all(r["event"] == "AMBIGUOUS" and r["track_id"] is None for r in crossing)
    assert all(set(r["candidate_track_ids"]) == {"sample-T001", "sample-T002"} for r in crossing)
    final = {r["track_id"]: r["center_cm1"] for r in observations(result, "s3")}
    assert final == {"sample-T001": 160, "sample-T002": 120}


def test_duplicate_and_decreasing_pressures_preserve_acquisition_identity():
    records = [frame(i, p, [10 + 2 * p]) for i, p in enumerate([3, 1, 1, 0, -1])]
    result = track_peaks(records, {"max_distance_cm1": 10})
    assert result["status"] == "PASS"
    assert len(result["tracks"]) == 1
    assert [r["pressure_gpa"] for r in result["assignments"]] == [3, 1, 1, 0, -1]


def test_interleaved_series_do_not_share_tracks_or_gap_counters():
    records = [
        frame(0, 0, [100], series="a"),
        frame(1, 0, [100], series="b"),
        frame(2, 1, [102], series="a"),
        frame(3, 1, [98], series="b"),
    ]
    result = track_peaks(records, {"max_gap": 0})
    assert result["status"] == "PASS"
    assert {r["track_id"] for r in result["tracks"]} == {"a-T001", "b-T001"}
    assert [r["spectrum_id"] for r in result["assignments"]] == ["s0", "s1", "s2", "s3"]


@pytest.mark.parametrize("status", ["FAIL", "INCONCLUSIVE", "BLOCKED"])
def test_untrusted_fits_never_advance_trajectories(status):
    result = track_peaks(
        [frame(0, 0, [100]), frame(1, 1, [500], status=status), frame(2, 2, [104])], {"max_gap": 1}
    )
    assert result["status"] == "INCONCLUSIVE"
    assert observations(result, "s1")[0]["event"] == "EXCLUDED_FIT"
    assert observations(result, "s2")[0]["track_id"] == "sample-T001"
    gap = [r for r in result["assignments"] if r["event"] == "GAP"]
    assert gap and all(r["center_cm1"] is None for r in gap)


def test_missing_pressure_stays_missing_and_gap_is_visible():
    result = track_peaks([frame(0, 0, [100]), frame(1, None, [900]), frame(2, 2, [104])])
    assert observations(result, "s1")[0]["pressure_gpa"] is None
    assert observations(result, "s1")[0]["event"] == "MISSING_PRESSURE"
    assert observations(result, "s2")[0]["track_id"] == "sample-T001"


def test_gap_limit_and_distance_gate_create_new_tracks():
    records = [frame(0, 0, [100]), frame(1, 1, []), frame(2, 2, [100])]
    result = track_peaks(records, {"max_gap": 0})
    assert observations(result, "s2")[0]["track_id"] == "sample-T002"
    far = track_peaks([frame(0, 0, [100]), frame(1, 1, [200])], {"max_distance_cm1": 5})
    assert observations(far, "s1")[0]["event"] == "NEW"
    assert len(far["tracks"]) == 2


def test_peak_birth_and_disappearance_are_not_forced_into_existing_tracks():
    result = track_peaks([frame(0, 0, [100]), frame(1, 1, [102, 300]), frame(2, 2, [104])])
    assert len(result["tracks"]) == 2
    assert observations(result, "s1")[1]["event"] == "NEW"
    assert any(
        r["event"] == "GAP" and r["track_id"] == "sample-T002" for r in result["assignments"]
    )


def test_global_ambiguity_margin_preserves_competing_assignments():
    records = [frame(0, 0, [100, 110]), frame(1, 1, [104, 106])]
    uncertain = track_peaks(records, {"ambiguity_margin_cm1": 5})
    assert all(r["event"] == "AMBIGUOUS" for r in observations(uncertain, "s1"))
    resolved = track_peaks(records, {"ambiguity_margin_cm1": 1})
    assert resolved["status"] == "PASS"


def test_near_gate_match_can_be_ambiguous_with_birth():
    result = track_peaks(
        [frame(0, 0, [100]), frame(1, 1, [119.5])],
        {"max_distance_cm1": 20, "ambiguity_margin_cm1": 1},
    )
    assert observations(result, "s1")[0]["event"] == "AMBIGUOUS"
    assert set(observations(result, "s1")[0]["candidate_track_ids"]) == {"sample-T001", "NEW"}


def test_assignment_candidates_match_exhaustive_global_oracle():
    # Literal 2-track / 2-observation cost matrix, including birth/gap alternatives.
    cost = np.array(
        [
            [4.0, 6.0, 10.0, np.inf],
            [6.0, 4.0, np.inf, 10.0],
            [10.0, np.inf, 0.0, 0.0],
            [np.inf, 10.0, 0.0, 0.0],
        ]
    )
    candidates = []
    for ordering in permutations(range(4)):
        total = sum(cost[i, j] for i, j in enumerate(ordering))
        if np.isfinite(total):
            candidates.append((total, ordering))
    optimum = min(c[0] for c in candidates)
    for margin in (0, 4, 13, 100):
        expected = []
        for j in range(2):
            possible = set()
            for total, ordering in candidates:
                if total <= optimum + margin:
                    i = ordering.index(j)
                    possible.add(i if i < 2 else None)
            expected.append(possible)
        actual = _possible_assignments(cost, 2, 2, margin)
        assert [{i for i, _ in choices} for choices in actual] == expected


def test_future_frames_do_not_change_previous_assignments():
    records = [frame(i, p, [100 + 2 * p]) for i, p in enumerate([0, 1, 2, 3])]
    prefix = track_peaks(records[:2])["assignments"]
    assert track_peaks(records)["assignments"][: len(prefix)] == prefix


@pytest.mark.parametrize(
    "settings",
    [
        {"max_gap": -1},
        {"max_distance_cm1": 0},
        {"width_weight": float("nan")},
        {"enabled": "false"},
        {"ambiguity_margin_cm1": True},
    ],
)
def test_invalid_tracking_settings_rejected(settings):
    with pytest.raises((ValueError, TypeError)):
        track_peaks([frame(0, 0, [100])], settings)


def test_disabled_tracking_is_explicit_and_missing_observations_are_inconclusive():
    assert track_peaks([frame(0, None, [100])], {"enabled": False})["status"] == "NOT_RUN"
    assert track_peaks([])["status"] == "INCONCLUSIVE"


def test_invalid_pressure_and_reused_identity_are_rejected():
    with pytest.raises(ValueError):
        track_peaks([frame(0, float("nan"), [100])])
    with pytest.raises(ValueError):
        track_peaks([frame(0, 0, [100]), frame(0, 1, [101])])


def test_explicit_width_cost_can_resolve_a_positional_tie():
    records = [frame(0, 0, [100, 110]), frame(1, 1, [105, 105])]
    for record in records:
        record["peaks"][1]["fwhm_cm1"] = 20
    unweighted = track_peaks(records, {"ambiguity_margin_cm1": 0})
    assert all(r["event"] == "AMBIGUOUS" for r in observations(unweighted, "s1"))
    weighted = track_peaks(records, {"ambiguity_margin_cm1": 0, "width_weight": 1})
    assert weighted["status"] == "PASS"
    assert [r["track_id"] for r in observations(weighted, "s1")] == ["sample-T001", "sample-T002"]
