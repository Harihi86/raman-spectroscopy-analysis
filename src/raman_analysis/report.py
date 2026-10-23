"""Offline bilingual result index; all displayed configuration text is escaped."""

import json
from html import escape
from pathlib import Path


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
            ("spectrum.png", "拟合图 / figure"),
        ]:
            relative = f"{spectrum['output']}/{name}"
            if (output / relative).is_file():
                artifacts.append(link(relative, label))
        diagnostics = "; ".join(spectrum["diagnostics"]) or "—"
        values = [
            spectrum["id"],
            spectrum["label"],
            spectrum["series"],
            pressure_text,
            spectrum["status"],
            str(len(spectrum["peaks"])),
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
table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{text-align:left;padding:10px;border-bottom:1px solid #ddd;vertical-align:top}}
th{{background:#edf3f5}}img{{max-width:100%;height:auto}}figure{{margin:14px 0}}code{{background:#edf3f5;padding:2px 5px}}
small{{color:#546a75}}li{{overflow-wrap:anywhere}}.note{{border-left:4px solid #b97825;padding-left:14px}}
</style></head><body><main>
<h1>Raman 批处理分析 <small>Batch analysis</small></h1>
<section class="panel"><div class="status">{escape(report["status"])}</div>
<p>数据类型 / Data kind: <strong>{escape(", ".join(kinds))}</strong><br>
光谱 / Spectra: {report["total"]} · 轨迹 / Tracks: {report["tracks"]} · 峰追踪 / Tracking: {escape(report["tracking_status"])}</p>
<p class="note">轨迹标签表示算法连续性，不直接代表化学或振动模式归属。合成数据不属于实验观测。<br>
Tracks are algorithmic continuity labels. Synthetic spectra are not experimental observations.</p>
<p>{link("summary.csv", "批处理汇总 / Summary CSV")} · {link("batch-report.json", "完整报告 / Full JSON")} ·
{link("manifest.json", "运行清单 / Manifest")} · {link("tracking/assignments.csv", "全部配对和缺口 / Assignments CSV")} ·
{link("tracking/tracks.csv", "轨迹表 / Tracks CSV")}</p></section>
<section class="panel"><h2>逐谱结果 / Per-spectrum results</h2>
<p>PASS：配置的数值检查通过；FAIL：拟合未通过；INCONCLUSIVE：结果存在未解决问题；BLOCKED：输入或配置阻碍分析。</p>
<div class="scroll"><table><thead><tr><th>ID</th><th>名称 / Label</th><th>序列 / Series</th><th>GPa</th><th>状态 / Status</th><th>峰数</th><th>诊断 / Diagnostics</th><th>结果 / Files</th></tr></thead>
<tbody>{"".join(rows)}</tbody></table></div></section>
<section class="panel"><h2>压力轨迹 / Pressure trajectories</h2>
<p>圆点为采用的拟合峰位；灰色叉号为歧义或排除结果；缺失帧处断线。连接线只辅助阅读，不增加测量点。误差棒是固定模型下的局部拟合误差。<br>
Markers are fitted peak estimates. Lines guide the eye; gaps and unassigned estimates remain explicit.</p>{gallery or "<p>本次未生成图 / Figures disabled.</p>"}</section>
<section class="panel"><h2>需复核事项 / Review required</h2><ul>{diagnostics or "<li>无追踪告警 / No tracking warnings.</li>"}</ul>
<p>追踪使用全局代价及歧义阈值，代价差不等于概率。请结合原始光谱、残差和实验分辨率复核。<br>
Assignment margins are heuristic costs, not probabilities. Inspect spectra, residuals and resolution.</p></section>
<small>Raman Spectroscopy Analysis {escape(report["program_version"])} · UTF-8 · 离线报告 / Offline report</small>
</main></body></html>"""
