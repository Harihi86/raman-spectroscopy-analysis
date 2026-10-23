import hashlib
import json

import numpy as np
import pandas as pd
import pytest
from jsonschema import Draft202012Validator

from raman_analysis.batch import analyze_batch, retrack
from raman_analysis.cli import main
from raman_analysis.examples import demo_batch
from raman_analysis.manifest import merge_settings, read_schema, validate_manifest
from raman_analysis.pipeline import analyze


def manifest_file(tmp_path, pressures=(0, 1, 1, 3), *, single_last=False):
    directory = tmp_path / "输入数据"
    directory.mkdir()
    entries = []
    x = np.linspace(0, 200, 401)
    rng = np.random.default_rng(75)
    for i, pressure in enumerate(pressures):
        first = 20 * np.exp(-0.5 * ((x - (50 + 2 * pressure)) / 6) ** 2)
        second = 14 * np.exp(-0.5 * ((x - (145 - pressure)) / 9) ** 2)
        y = 3 + 0.005 * x + first + rng.normal(0, 0.025, len(x))
        if not (single_last and i == len(pressures) - 1):
            y += second
        path = directory / f"光谱 {i}.csv"
        pd.DataFrame({"x": x, "counts": y, "error": np.full(len(x), 0.025)}).to_csv(
            path, index=False
        )
        entry = {
            "id": f"s{i}",
            "label": f"光谱 {i}",
            "pressure_gpa": pressure,
            "input": {"path": path.relative_to(tmp_path).as_posix()},
        }
        if single_last and i == len(pressures) - 1:
            entry["analysis"] = {"fit": {"peaks": [{"center": 50, "sigma": 6}]}}
        entries.append(entry)
    manifest = {
        "schema_version": 1,
        "data_kind": "synthetic",
        "defaults": {
            "input": {"x": "x", "y": "counts", "sigma": "error"},
            "baseline": {"method": "none"},
            "fit": {"peaks": [{"center": 48, "sigma": 5}, {"center": 143, "sigma": 8}]},
        },
        "tracking": {"max_distance_cm1": 20, "ambiguity_margin_cm1": 0.5},
        "outputs": {"per_spectrum_figures": False, "tracking_figures": False},
        "spectra": entries,
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False), encoding="utf-8")
    return path, manifest


def save(path, manifest):
    path.write_text(json.dumps(manifest, ensure_ascii=False), encoding="utf-8")


def test_batch_runs_from_another_directory_and_preserves_raw_inputs(tmp_path, monkeypatch):
    path, _ = manifest_file(tmp_path)
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.rglob("*.csv")}
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    output = tmp_path / "results"
    report = analyze_batch(path, output)
    assert report["state"] == "COMPLETE" and report["status"] == "PASS"
    assert report["total"] == 4 and report["tracks"] == 2
    assert [s["pressure_gpa"] for s in report["spectra"]] == [0, 1, 1, 3]
    for source, digest in hashes.items():
        assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
    assert (output / "summary.csv").read_bytes().startswith(b"\xef\xbb\xbf")
    assert "光谱 0" in (output / "report.html").read_text(encoding="utf-8")
    for name in [
        "batch-report.json",
        "run.json",
        "spectra.jsonl",
        "resolved-spectra.json",
        "tracking/tracking.json",
        "tracking/assignments.csv",
        "tracking/tracks.csv",
        "artifacts.json",
    ]:
        assert (output / name).stat().st_size > 0


def test_missing_input_is_isolated_and_never_counted_as_success(tmp_path):
    path, manifest = manifest_file(tmp_path)
    manifest["spectra"][1]["input"]["path"] = "not-found.csv"
    save(path, manifest)
    result = analyze_batch(path, tmp_path / "output")
    assert result["status"] == "INCONCLUSIVE"
    assert result["fit_counts"] == {"PASS": 3, "BLOCKED": 1}
    assert result["spectra"][1]["status"] == "BLOCKED"
    assert all(s["status"] == "PASS" for i, s in enumerate(result["spectra"]) if i != 1)


def test_corrupt_workbook_is_isolated_from_other_spectra(tmp_path):
    path, manifest = manifest_file(tmp_path)
    (tmp_path / "corrupt.xlsx").write_bytes(b"This is not an Excel archive")
    manifest["spectra"][1]["input"]["path"] = "corrupt.xlsx"
    save(path, manifest)
    result = analyze_batch(path, tmp_path / "output")
    assert result["fit_counts"] == {"PASS": 3, "BLOCKED": 1}
    assert "Cannot read XLSX" in result["spectra"][1]["diagnostics"][0]


def test_shared_excel_columns_match_independent_single_spectrum_runs(tmp_path):
    path, manifest = manifest_file(tmp_path)
    manifest["spectra"] = manifest["spectra"][:2]
    workbook = pd.read_csv(tmp_path / manifest["spectra"][0]["input"]["path"])
    for i, spectrum in enumerate(manifest["spectra"]):
        workbook[f"counts_{i}"] = pd.read_csv(tmp_path / spectrum["input"]["path"])["counts"]
        spectrum["input"] = {"path": "shared.xlsx", "sheet": "采集", "y": f"counts_{i}"}
    workbook.to_excel(tmp_path / "shared.xlsx", sheet_name="采集", index=False)
    save(path, manifest)
    result = analyze_batch(path, tmp_path / "batch")
    assert result["status"] == "PASS"
    for i, spectrum in enumerate(result["spectra"]):
        config = {
            "data_kind": "synthetic",
            "input": {
                "path": "shared.xlsx",
                "sheet": "采集",
                "x": "x",
                "y": f"counts_{i}",
                "sigma": "error",
            },
            "baseline": {"method": "none"},
            "fit": {"peaks": [{"center": 48, "sigma": 5}, {"center": 143, "sigma": 8}]},
        }
        single = analyze(config, tmp_path / f"single-{i}", base=tmp_path, save_figures=False)
        assert spectrum["peaks"] == single["peaks"]
        assert spectrum["provenance"] == single["provenance"]


def test_per_spectrum_peak_count_override_replaces_default_list(tmp_path):
    path, _ = manifest_file(tmp_path, single_last=True)
    result = analyze_batch(path, tmp_path / "output")
    assert result["status"] == "PASS"
    assert [len(s["peaks"]) for s in result["spectra"]] == [2, 2, 2, 1]


def test_batch_without_pressure_can_disable_tracking_explicitly(tmp_path):
    path, manifest = manifest_file(tmp_path)
    for spectrum in manifest["spectra"]:
        spectrum.pop("pressure_gpa")
    manifest["tracking"] = {"enabled": False}
    save(path, manifest)
    result = analyze_batch(path, tmp_path / "output")
    assert result["status"] == "PASS"
    assert result["tracking_status"] == "NOT_RUN" and result["tracks"] == 0


@pytest.mark.parametrize(
    "kind",
    [
        "duplicate-id",
        "case-duplicate",
        "escape-id",
        "reserved-id",
        "unknown-field",
        "nonfinite-pressure",
        "negative-column",
    ],
)
def test_invalid_manifest_fails_before_creating_output(tmp_path, kind):
    path, manifest = manifest_file(tmp_path)
    if kind == "duplicate-id":
        manifest["spectra"][1]["id"] = "s0"
    elif kind == "case-duplicate":
        manifest["spectra"][1]["id"] = "S0"
    elif kind == "escape-id":
        manifest["spectra"][0]["id"] = "../escape"
    elif kind == "reserved-id":
        manifest["spectra"][0]["id"] = "CON"
    elif kind == "unknown-field":
        manifest["tracking"]["max_distnace"] = 5
    elif kind == "negative-column":
        manifest["defaults"]["input"]["x"] = -1
    else:
        manifest["spectra"][0]["pressure_gpa"] = float("nan")
    save(path, manifest)
    with pytest.raises(ValueError):
        analyze_batch(path, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_changing_baseline_method_replaces_incompatible_parameters():
    merged = merge_settings(
        {"baseline": {"method": "asls", "lam": 1e6}}, {"baseline": {"method": "none"}}
    )
    assert merged["baseline"] == {"method": "none"}


def test_interruption_preserves_completed_records_and_never_reports_pass(tmp_path):
    path, _ = manifest_file(tmp_path)
    output = tmp_path / "output"

    def interrupt(*args):
        raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        analyze_batch(path, output, progress=interrupt)
    state = json.loads((output / "run.json").read_text(encoding="utf-8"))
    assert state["state"] == "INTERRUPTED" and state["status"] == "INCONCLUSIVE"
    assert state["completed"] == 1
    assert len((output / "spectra.jsonl").read_text(encoding="utf-8").splitlines()) == 1
    assert not (output / "batch-report.json").exists()
    with pytest.raises(FileExistsError):
        analyze_batch(path, output)


def test_retrack_uses_saved_fits_and_does_not_rerun_analysis(tmp_path, monkeypatch):
    path, _ = manifest_file(tmp_path)
    analyze_batch(path, tmp_path / "output")

    def forbidden(*args, **kwargs):
        raise AssertionError("Refitting was not requested")

    monkeypatch.setattr("raman_analysis.batch.analyze", forbidden)
    settings = tmp_path / "settings.json"
    settings.write_text('{"prediction":"nearest", "max_distance_cm1":20}')
    report_path = tmp_path / "output/batch-report.json"
    result = retrack(report_path, tmp_path / "retracked", settings_path=settings, figures=False)
    assert result["status"] == "PASS"
    assert result["source_batch_sha256"] == hashlib.sha256(report_path.read_bytes()).hexdigest()


def test_retrack_rejects_incomplete_report(tmp_path):
    path = tmp_path / "partial.json"
    path.write_text('{"format":"raman-analysis-batch", "schema_version":1, "state":"RUNNING"}')
    with pytest.raises(ValueError, match="complete"):
        retrack(path, tmp_path / "out", figures=False)


def test_chinese_help_and_json_summary(tmp_path, capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--lang", "zh", "--help"])
    assert exc.value.code == 0
    assert "跨压力峰追踪" in capsys.readouterr().out
    path, _ = manifest_file(tmp_path)
    assert main(["validate", str(path), "--check-files", "--lang", "zh", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["input_files_checked"] is True


def test_cli_partial_batch_has_nonzero_exit_and_readable_results(tmp_path, capsys):
    path, manifest = manifest_file(tmp_path)
    manifest["spectra"][0]["input"]["path"] = "missing.csv"
    save(path, manifest)
    assert main(["--lang", "zh", "batch", str(path), "--output", str(tmp_path / "output")]) == 1
    assert "INCONCLUSIVE" in capsys.readouterr().out


def test_html_escapes_user_labels(tmp_path):
    path, manifest = manifest_file(tmp_path)
    manifest["spectra"][0]["label"] = "<script>bad()</script>"
    save(path, manifest)
    analyze_batch(path, tmp_path / "output")
    html = (tmp_path / "output/report.html").read_text(encoding="utf-8")
    assert "<script>bad()" not in html and "&lt;script&gt;" in html


def test_bundled_schemas_remain_consistent(tmp_path):
    schema = read_schema()
    standalone = read_schema("tracking-settings.schema.json")
    Draft202012Validator.check_schema(schema)
    Draft202012Validator.check_schema(standalone)
    assert schema["$defs"]["tracking"] == {
        k: v for k, v in standalone.items() if k not in {"$schema", "title"}
    }
    _, manifest = manifest_file(tmp_path)
    assert len(validate_manifest(manifest)) == 4


@pytest.mark.parametrize(
    "scenario,expected_spectra,expected_tracks", [("basic", 11, 4), ("crossing", 6, 2)]
)
def test_end_to_end_synthetic_series_matches_independent_truth(
    tmp_path, scenario, expected_spectra, expected_tracks
):
    output = tmp_path / scenario
    report = demo_batch(output, scenario=scenario, no_figures=True)
    assert report["status"] == "PASS"
    assert report["total"] == expected_spectra and report["tracks"] == expected_tracks
    truth = json.loads((output / "truth.json").read_text(encoding="utf-8"))["spectra"]
    assignments = json.loads(
        (output / "analysis/tracking/tracking.json").read_text(encoding="utf-8")
    )["assignments"]
    for source in truth:
        rows = [
            r
            for r in assignments
            if r["spectrum_id"] == source["spectrum_id"] and r["event"] in {"NEW", "MATCHED"}
        ]
        rows.sort(key=lambda r: r["track_id"])
        np.testing.assert_allclose(
            [r["center_cm1"] for r in rows],
            source["centers_in_physical_branch_order_cm1"],
            atol=0.4,
        )


@pytest.mark.parametrize("export_name", ["batch_html", "write_artifacts"])
def test_final_export_failure_preserves_journal_without_completed_report(
    tmp_path, monkeypatch, export_name
):
    path, _ = manifest_file(tmp_path)
    output = tmp_path / "output"
    from raman_analysis import batch

    original = getattr(batch, export_name)

    def fail_export(*args, **kwargs):
        if export_name == "write_artifacts" and args[0] != output:
            return original(*args, **kwargs)
        raise ValueError("Export failed")

    monkeypatch.setattr("raman_analysis.batch." + export_name, fail_export)
    with pytest.raises(ValueError, match="Export failed"):
        analyze_batch(path, output)
    state = json.loads((output / "run.json").read_text(encoding="utf-8"))
    assert state["state"] == "ABORTED" and state["status"] == "INCONCLUSIVE"
    assert len((output / "spectra.jsonl").read_text(encoding="utf-8").splitlines()) == 4
    assert not (output / "batch-report.json").exists()
