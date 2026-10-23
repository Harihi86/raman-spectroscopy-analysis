"""Compare saved counts across declared baselines/windows without refitting.

Information criteria from different corrected observations are never ranked.
Agreement is separated from diagnostic support and physical-count accuracy.
"""

import hashlib
import json
from collections import Counter
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd


def _read_verified_report(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    report = json.loads(raw)
    if report.get("format") != "raman-analysis-batch" or report.get("state") != "COMPLETE":
        raise ValueError("Stability requires complete batch reports")
    for manifest in path.parent.rglob("artifacts.json"):
        for name, digest in json.loads(manifest.read_text(encoding="utf-8")).items():
            target = (manifest.parent / name).resolve()
            if not target.is_relative_to(path.parent):
                raise ValueError("Artifact path escapes the batch directory")
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError(f"Artifact hash mismatch: {name}")
    records = {}
    for record in report["spectra"]:
        if record["id"] in records:
            raise ValueError("Duplicate spectrum ID in saved report")
        folder = (path.parent / record["output"]).resolve()
        if not folder.is_relative_to(path.parent):
            raise ValueError("Spectrum output escapes the batch directory")
        detail = json.loads((folder / "report.json").read_text(encoding="utf-8"))
        if detail["status"] != record["status"]:
            raise ValueError("Saved spectrum status disagrees with the batch report")
        records[record["id"]] = (record, detail, folder)
    return path, hashlib.sha256(raw).hexdigest(), records


def compare_saved_batches(report_paths, output):
    """All IDs must match; blocked or missing models are retained as unresolved."""
    if len(report_paths) < 2:
        raise ValueError("Provide at least two independently configured batch reports")
    loaded = [_read_verified_report(p) for p in report_paths]
    if len({p for p, _, _ in loaded}) != len(loaded):
        raise ValueError("Repeated reports do not provide independent configuration variants")
    ids = list(loaded[0][2])
    if not ids:
        raise ValueError("Compared batches must contain at least one spectrum")
    if any(set(records) != set(ids) for _, _, records in loaded):
        raise ValueError("Compared batches must contain exactly the same spectrum IDs")
    variants = [
        {"variant": i + 1, "report_name": path.parent.name, "report_sha256": digest}
        for i, (path, digest, _) in enumerate(loaded)
    ]
    outcomes, rows = [], []
    for spectrum_id in ids:
        observations, identity, first_raw, settings, kind = [], None, None, None, None
        for i, (_, _, records) in enumerate(loaded):
            record, detail, folder = records[spectrum_id]
            kind = kind or record["data_kind"]
            if kind != record["data_kind"]:
                raise ValueError("Compared data kinds differ")
            if record["status"] == "BLOCKED":
                observations.append(
                    {
                        "variant": i + 1,
                        "status": "BLOCKED",
                        "selected_peak_count": None,
                        "families": [],
                        "reason": record["diagnostics"],
                    }
                )
                continue
            provenance = detail["provenance"]
            current_identity = {
                k: provenance[k]
                for k in (
                    "input_sha256",
                    "header",
                    "units",
                    "x_column",
                    "y_column",
                    "sigma_column",
                    "sheet",
                    "missing_policy",
                    "dropped_rows",
                )
            }
            if identity is not None and identity != current_identity:
                raise ValueError(
                    f"Spectrum {spectrum_id}: input provenance differs across variants"
                )
            identity = current_identity
            raw = pd.read_csv(folder / "processed.csv")[["shift_cm1", "raw_intensity"]].to_numpy()
            if first_raw is not None and (
                first_raw.shape != raw.shape
                or not np.allclose(first_raw, raw, rtol=1e-12, atol=1e-12)
            ):
                raise ValueError(f"Spectrum {spectrum_id}: raw observations differ across variants")
            first_raw = raw
            decision = detail.get("model_selection")
            if not decision or "shape_comparison" not in decision:
                raise ValueError(
                    "Stability needs automatic family-comparison outputs from version 0.4 or later"
                )
            # Only baseline and declared fit window may differ. No threshold,
            # shape-set, iteration budget, or uncertainty-policy drift is hidden.
            controls = dict(detail["fit_configuration"])
            controls.pop("window", None)
            current_settings = (decision["settings"], controls)
            if settings is not None and settings != current_settings:
                raise ValueError("Automatic settings or fit budgets differ across variants")
            settings = current_settings
            families = [
                {
                    k: family[k]
                    for k in (
                        "shape",
                        "selected_peak_count",
                        "plausible_peak_counts",
                        "status",
                        "search_complete",
                    )
                }
                for family in decision["shape_comparison"]
            ]
            observations.append(
                {
                    "variant": i + 1,
                    "status": detail["status"],
                    "baseline": detail["baseline"],
                    "window": detail["fit_configuration"].get("window"),
                    "baseline_solver_converged": detail["baseline_diagnostics"].get("converged"),
                    "selected_peak_count": decision["selected_peak_count"],
                    "selected_shape": decision["selected_shape"],
                    "plausible_peak_counts": decision["plausible_peak_counts"],
                    "families": families,
                }
            )
        counts = sorted(
            {o["selected_peak_count"] for o in observations if o["selected_peak_count"] is not None}
            | {f["selected_peak_count"] for o in observations for f in o["families"]}
        )
        incomplete = any(o["status"] == "BLOCKED" for o in observations)
        supported = not incomplete and all(
            o["status"] == "PASS"
            and all(
                f["status"] == "PASS"
                and f["search_complete"]
                and len(f["plausible_peak_counts"]) == 1
                for f in o["families"]
            )
            for o in observations
        )
        stability = (
            "UNRESOLVED"
            if incomplete
            else (
                "SENSITIVE"
                if len(counts) > 1
                else ("SUPPORTED_STABLE" if supported else "AGREEMENT_ONLY")
            )
        )
        outcome = {
            "id": spectrum_id,
            "data_kind": kind,
            "stability": stability,
            "observed_best_counts": counts,
            "variants": observations,
        }
        outcomes.append(outcome)
        for observation in observations:
            rows.append(
                {
                    "id": spectrum_id,
                    "stability": stability,
                    "variant": observation["variant"],
                    "selected_peak_count": observation["selected_peak_count"],
                    "selected_shape": observation.get("selected_shape"),
                    "fit_status": observation["status"],
                    "family_best_counts": json.dumps(
                        {f["shape"]: f["selected_peak_count"] for f in observation["families"]},
                        sort_keys=True,
                    ),
                    "baseline": json.dumps(observation.get("baseline")),
                    "window": json.dumps(observation.get("window")),
                }
            )
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    totals = dict(Counter(o["stability"] for o in outcomes))
    result = {
        "format": "raman-count-stability",
        "schema_version": 1,
        "state": "COMPLETE",
        "status": "PASS"
        if all(o["stability"] == "SUPPORTED_STABLE" for o in outcomes)
        else "INCONCLUSIVE",
        "spectra": len(ids),
        "variant_count": len(loaded),
        "stability_counts": totals,
        "variants": variants,
        "records": outcomes,
        "physical_peak_count_accuracy": "NOT_ASSESSED",
        "scope": "Sensitivity of configured component counts; no cross-baseline information-criterion ranking, no true-count labels, no propagated baseline/model uncertainty",
    }
    (output / "stability.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8"
    )
    pd.DataFrame(rows).to_csv(output / "stability.csv", index=False, encoding="utf-8-sig")
    table = "".join(
        "<tr><td>"
        + escape(o["id"])
        + "</td><td>"
        + o["stability"]
        + "</td><td>"
        + escape(str(o["observed_best_counts"]))
        + "</td></tr>"
        for o in outcomes
    )
    variant_list = "".join(
        "<li>" + escape(f"{v['variant']}: {v['report_name']} ({v['report_sha256']})") + "</li>"
        for v in variants
    )
    (output / "report.html").write_text(
        f"""<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>峰数稳定性 / Count stability</title><style>body{{font:16px/1.7 system-ui;max-width:1050px;margin:32px auto;padding:0 20px;color:#233742}}table{{width:100%;border-collapse:collapse}}td,th{{padding:8px;text-align:left;border-bottom:1px solid #ddd}}li{{overflow-wrap:anywhere}}.scroll{{overflow:auto}}</style>
<h1>峰数稳定性 / Count stability</h1><p>{result["status"]} · {len(ids)} 条光谱 · {len(loaded)} 个配置变体</p>
<p>{escape(str(totals))}</p><p>SUPPORTED_STABLE：所有变体及每种峰型的最佳峰数一致，且全部诊断通过。AGREEMENT_ONLY：数目相同，但诊断仍有问题。SENSITIVE：峰数发生改变。UNRESOLVED：存在无法运行的记录。</p>
<p>数目一致不等于物理峰数准确。不同基线下的 BIC/AICc 没有相互排名；基线与峰型分歧不会被隐藏成零误差。</p>
<p><a href="stability.csv">逐谱逐配置 CSV</a> · <a href="stability.json">完整判断与来源 JSON</a></p>
<div class="scroll"><table><tr><th>光谱 / Spectrum</th><th>稳定性 / Stability</th><th>出现的最佳峰数 / Best counts</th></tr>{table}</table></div>
<h2>比较来源 / Sources</h2><ol>{variant_list}</ol></html>""",
        encoding="utf-8",
    )
    hashes = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file()
    }
    (output / "artifacts.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf-8")
    return result
