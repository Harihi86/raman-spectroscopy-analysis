"""Offline bilingual result index; all displayed configuration text is escaped."""

import json
from html import escape
from pathlib import Path


def selection_html(selection: dict, status: str, diagnostics: list[str]) -> str:
    """A readable offline score table, separate from machine-readable traces."""
    rows = []
    for candidate in selection["candidates"]:
        count = candidate["peak_count"]
        adopted = count == selection["selected_peak_count"]
        cells = [str(count) + (" ★" if adopted else ""), candidate["fit_status"]]
        for name in ("bic", "aicc", "delta"):
            value = candidate.get(name)
            cells.append("未得到 / unavailable" if value is None else f"{value:.4f}")
        cells.append(str(len(candidate["attempts"])))
        rows.append("<tr>" + "".join("<td>" + escape(v) + "</td>" for v in cells) + "</tr>")
    warnings = "".join("<li>" + escape(d) + "</li>" for d in diagnostics)
    settings = escape(json.dumps(selection["settings"], indent=2, ensure_ascii=False))
    return f"""<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>自动峰数判断 / Peak-count selection</title>
<style>body{{font:16px/1.7 system-ui,sans-serif;max-width:1000px;margin:32px auto;padding:0 20px;color:#233742}}
table{{border-collapse:collapse;width:100%}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #ddd}}
th{{background:#edf3f5}}a{{color:#12696d}}pre{{white-space:pre-wrap;background:#f5f7f8;padding:16px}}
.note{{border-left:4px solid #b97825;padding-left:16px}}.scroll{{overflow:auto}}</style>
<h1>自动峰数判断 / Peak-count selection</h1>
<p><strong>{escape(status)}</strong> · 最佳候选 / Best count: {selection["selected_peak_count"]}<br>
竞争峰数 / Competing counts: {escape(str(selection["plausible_peak_counts"]))} · 判据 / Criterion: {escape(selection["criterion"])}</p>
<p class="note">★ 仅表示当前模型与基线下代价最小的候选；INCONCLUSIVE 表示判断仍有未解决问题。代价差不是概率。<br>
The best candidate does not establish a physical mode count. Score differences are not probabilities.</p>
<div class="scroll"><table><tr><th>峰数 / Count</th><th>拟合状态 / Fit</th><th>BIC</th><th>AICc</th><th>代价差 / Delta</th><th>尝试 / Starts</th></tr>{"".join(rows)}</table></div>
<h2>需要复核 / Review required</h2><ul>{warnings or "<li>配置的诊断通过 / Configured diagnostics passed.</li>"}</ul>
<p>请结合原始谱、残差、基线敏感性、分辨率和物理背景判断。零峰不是“没有振动模式”的证明；达到峰数上限不代表已确认峰数。</p>
<p><a href="model-selection.csv">峰数比较 CSV</a> · <a href="model-selection.json">完整搜索记录 JSON</a> ·
<a href="processed.csv">原始强度与基线 CSV</a> · <a href="fit.csv">拟合和残差 CSV</a> · <a href="report.json">完整拟合报告 JSON</a></p>
<details><summary>自动选峰参数 / Settings</summary><pre>{settings}</pre></details>
<p>所有候选拟合相同的原始采样点；平滑只用于初值。局部协方差不包含基线和峰数选择误差。<br>
All models fit the same original observations; smoothing is for initialization only. Local errors exclude baseline and selection uncertainty.</p>
</html>"""


def batch_html(report: dict, tracking: dict, output: Path) -> str:
    def link(path, label):
        return f'<a href="{escape(path, quote=True)}">{escape(label)}</a>'

    rows = []
    for spectrum in report["spectra"]:
        pressure = spectrum["pressure_gpa"]
        pressure_text = "未提供 / missing" if pressure is None else f"{pressure:g}"
        artifacts = []
        for name, label in [
            ("report.json", "详情 / JSON"),
            ("processed.csv", "光谱 / CSV"),
            ("peaks.csv", "峰参数 / peaks"),
            ("fit.csv", "曲线与残差 / residuals"),
            ("model-selection.html", "峰数比较 / model selection"),
            ("spectrum.png", "拟合图 / figure"),
        ]:
            relative = f"{spectrum['output']}/{name}"
            if (output / relative).is_file():
                artifacts.append(link(relative, label))
        diagnostics = "; ".join(spectrum["diagnostics"]) or "—"
        decision = spectrum.get("model_selection")
        count_text = str(len(spectrum["peaks"]))
        if decision:
            count_text += " / " + ", ".join(str(n) for n in decision["plausible_peak_counts"])
        values = [
            spectrum["id"],
            spectrum["label"],
            spectrum["series"],
            pressure_text,
            spectrum["status"],
            count_text,
            diagnostics,
        ]
        rows.append(
            "<tr>"
            + "".join(f"<td>{escape(str(v))}</td>" for v in values)
            + "<td>"
            + " · ".join(artifacts)
            + "</td></tr>"
        )
    diagnostics = "".join(
        "<li>" + escape(json.dumps(d, ensure_ascii=False)) + "</li>"
        for d in tracking["diagnostics"]
    )
    figures = json.loads((output / "tracking/figures.json").read_text(encoding="utf-8"))
    gallery = "".join(
        f'<figure><img src="tracking/{escape(f["png"], quote=True)}" alt="{escape(f["series"], quote=True)}">'
        f"<figcaption>{escape(f['series'])} — "
        + link("tracking/" + f["pdf"], "矢量 PDF / vector PDF")
        + "</figcaption></figure>"
        for f in figures
    )
    kinds = sorted({s["data_kind"] for s in report["spectra"]})
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Raman 批处理结果 / Batch analysis</title>
<style>
body{{font:16px/1.6 system-ui,sans-serif;color:#233742;background:#f5f7f8;margin:0}}
main{{max-width:1200px;margin:auto;padding:32px}}h1,h2{{line-height:1.3}}
.panel{{background:white;border:1px solid #d4dfe2;border-radius:10px;padding:22px;margin:18px 0}}
.status{{font-size:24px;font-weight:700;color:#12696d}}a{{color:#12696d}}.scroll{{overflow:auto}}
.status.INCONCLUSIVE,.status.BLOCKED{{color:#965a0b}}.status.FAIL{{color:#a32828}}
table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{text-align:left;padding:10px;border-bottom:1px solid #ddd;vertical-align:top}}
th{{background:#edf3f5}}img{{max-width:100%;height:auto}}figure{{margin:14px 0}}code{{background:#edf3f5;padding:2px 5px}}
small{{color:#546a75}}li{{overflow-wrap:anywhere}}.note{{border-left:4px solid #b97825;padding-left:14px}}
</style></head><body><main>
<h1>Raman 批处理分析 <small>Batch analysis</small></h1>
<section class="panel"><div class="status {escape(report["status"])}">{escape(report["status"])}</div>
<p>数据类型 / Data kind: <strong>{escape(", ".join(kinds))}</strong><br>
光谱 / Spectra: {report["total"]} · 轨迹 / Tracks: {report["tracks"]} · 峰追踪 / Tracking: {escape(report["tracking_status"])}</p>
<p class="note">轨迹标签表示算法连续性，不直接代表化学或振动模式归属。合成数据不属于实验观测。<br>
Tracks are algorithmic continuity labels. Synthetic spectra are not experimental observations.</p>
<p>{link("summary.csv", "批处理汇总 / Summary CSV")} · {link("batch-report.json", "完整报告 / Full JSON")} ·
{link("manifest.json", "运行清单 / Manifest")} · {link("tracking/assignments.csv", "全部配对和缺口 / Assignments CSV")} ·
{link("tracking/tracks.csv", "轨迹表 / Tracks CSV")}</p></section>
<section class="panel"><h2>逐谱结果 / Per-spectrum results</h2>
<p>PASS：配置的数值检查通过；FAIL：拟合未通过；INCONCLUSIVE：结果存在未解决问题；BLOCKED：输入或配置阻碍分析。</p>
<p>自动模式的峰数列显示“最佳候选 / 竞争峰数”。最佳候选不等于已确认的物理峰数；详见峰数比较和残差。<br>Automatic counts show best / competing counts; rankings do not establish physical mode counts.</p>
<div class="scroll"><table><thead><tr><th>ID</th><th>名称 / Label</th><th>序列 / Series</th><th>GPa</th><th>状态 / Status</th><th>峰数 / Counts</th><th>诊断 / Diagnostics</th><th>结果 / Files</th></tr></thead>
<tbody>{"".join(rows)}</tbody></table></div></section>
<section class="panel"><h2>压力轨迹 / Pressure trajectories</h2>
<p>圆点为采用的拟合峰位；灰色叉号为歧义或排除结果；缺失帧处断线。连接线只辅助阅读，不增加测量点。误差棒是固定模型下的局部拟合误差。<br>
Markers are fitted peak estimates. Lines guide the eye; gaps and unassigned estimates remain explicit.</p>{gallery or "<p>本次未生成图 / Figures disabled.</p>"}</section>
<section class="panel"><h2>需复核事项 / Review required</h2><ul>{diagnostics or "<li>无追踪告警 / No tracking warnings.</li>"}</ul>
<p>追踪使用全局代价及歧义阈值，代价差不等于概率。请结合原始光谱、残差和实验分辨率复核。<br>
Assignment margins are heuristic costs, not probabilities. Inspect spectra, residuals and resolution.</p></section>
<small>Raman Spectroscopy Analysis {escape(report["program_version"])} · UTF-8 · 离线报告 / Offline report</small>
</main></body></html>"""
