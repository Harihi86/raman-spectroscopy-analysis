import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .batch import analyze_batch, retrack
from .examples import demo_batch
from .fitting import gaussian
from .io import load_spectrum
from .manifest import validate_manifest
from .pipeline import analyze
from .stability import compare_saved_batches


def demo(output: Path):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    x = np.linspace(100, 1100, 1001)
    rng = np.random.default_rng(20261023)
    y = 3 + 0.003 * x + gaussian(x, 18, 350, 22) + gaussian(x, 12, 720, 35)
    y += rng.normal(0, 0.15, len(x))
    pd.DataFrame({"shift_cm1": x, "intensity": y, "noise_sd": np.full(len(x), 0.15)}).to_csv(
        output / "synthetic.csv", index=False
    )
    config = {
        "data_kind": "synthetic",
        "input": {"path": "synthetic.csv", "sigma": "noise_sd"},
        "baseline": {"method": "asls", "lam": 1e9, "p": 0.01, "iterations": 30},
        "fit": {"peaks": [{"center": 340, "sigma": 20}, {"center": 710, "sigma": 30}]},
        "metadata": {
            "seed": 20261023,
            "purpose": "software demonstration; no experimental observations",
        },
    }
    (output / "config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    (output / "truth.json").write_text(
        json.dumps(
            {
                "kind": "synthetic",
                "centers_cm1": [350, 720],
                "sigma_cm1": [22, 35],
                "amplitudes": [18, 12],
                "noise_sd": 0.15,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return analyze(config, output / "analysis", base=output)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--lang", choices=["en", "zh"], default="en")
    pre.add_argument("--json", action="store_true")
    options, remaining = pre.parse_known_args(argv)
    zh = options.lang == "zh"

    def text(english, chinese):
        return chinese if zh else english

    parser = argparse.ArgumentParser(
        description=text(
            "Reproducible Raman spectra, batches and pressure-peak tracking",
            "可复现的 Raman 单谱分析、批处理与跨压力峰追踪",
        )
    )
    parser.add_argument(
        "--lang",
        choices=["en", "zh"],
        default=options.lang,
        help=text("Help and console language", "帮助及终端显示语言"),
    )
    parser.add_argument(
        "--json",
        action="store_true",
        default=options.json,
        help=text("Write a machine-readable summary to stdout", "在标准输出中显示 JSON 摘要"),
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="operation", required=True)
    example = sub.add_parser(
        "demo", help=text("Run a synthetic single-spectrum example", "运行合成单谱示例")
    )
    run = sub.add_parser("run", help=text("Analyze one spectrum", "分析单条光谱"))
    run.add_argument(
        "config", type=Path, help=text("Single-spectrum JSON configuration", "单谱 JSON 配置")
    )
    batch = sub.add_parser(
        "batch",
        help=text("Fit a manifest and track pressure-dependent peaks", "批量拟合并追踪跨压力峰"),
    )
    batch.add_argument(
        "manifest", type=Path, help=text("Version-1 batch manifest", "版本 1 批处理清单")
    )
    demonstration = sub.add_parser(
        "demo-batch", help=text("Run synthetic pressure-series examples", "运行合成压力序列示例")
    )
    demonstration.add_argument(
        "--automatic",
        action="store_true",
        help=text("Select peak seeds and counts from the signal", "从信号自动选峰和判断峰数"),
    )
    demonstration.add_argument(
        "--scenario",
        choices=["basic", "crossing"],
        default="basic",
        help=text("Separate sweeps or crossing peaks", "独立压力序列或交叉峰场景"),
    )
    track = sub.add_parser(
        "track", help=text("Retrack saved fits without refitting", "复用已保存拟合结果，重新追踪")
    )
    track.add_argument(
        "report", type=Path, help=text("Completed batch-report.json", "已完成的 batch-report.json")
    )
    track.add_argument(
        "--settings", type=Path, help=text("Optional tracking settings JSON", "可选的追踪参数 JSON")
    )
    validation = sub.add_parser(
        "validate", help=text("Validate a batch manifest before execution", "运行前校验批处理清单")
    )
    validation.add_argument("manifest", type=Path)
    stability = sub.add_parser(
        "stability",
        help=text(
            "Compare saved peak counts across configurations",
            "比较不同基线或窗口下保存的峰数稳定性",
        ),
    )
    stability.add_argument(
        "reports",
        nargs="+",
        type=Path,
        help=text(
            "Complete batch reports for the same spectra", "同一批光谱的多个完整 batch-report.json"
        ),
    )
    validation.add_argument(
        "--check-files",
        action="store_true",
        help=text(
            "Also load and validate every spectrum input", "同时读取并校验各条光谱的输入文件"
        ),
    )
    for command in (example, run, batch, demonstration, track, stability):
        command.add_argument(
            "--output",
            required=True,
            type=Path,
            help=text(
                "New output directory; existing directories are preserved",
                "新的输出目录；不会覆盖已有目录",
            ),
        )
    for command in (batch, demonstration, track):
        command.add_argument(
            "--no-figures",
            action="store_true",
            help=text("Save numerical results without figures", "仅保存数值结果，不绘图"),
        )
    args = parser.parse_args(remaining)

    def progress(done, total, spectrum_id, status):
        print(f"[{done}/{total}] {spectrum_id}: {status}", file=sys.stderr)

    try:
        if args.operation == "demo":
            report = demo(args.output)
        elif args.operation == "demo-batch":
            report = demo_batch(
                args.output,
                scenario=args.scenario,
                automatic=args.automatic,
                no_figures=args.no_figures,
                progress=progress,
            )
        elif args.operation == "batch":
            report = analyze_batch(
                args.manifest, args.output, no_figures=args.no_figures, progress=progress
            )
        elif args.operation == "track":
            report = retrack(
                args.report, args.output, settings_path=args.settings, figures=not args.no_figures
            )
        elif args.operation == "stability":
            report = compare_saved_batches(args.reports, args.output)
        elif args.operation == "validate":
            manifest = json.loads(args.manifest.read_text(encoding="utf-8-sig"))
            spectra = validate_manifest(manifest)
            if args.check_files:
                for spectrum in spectra:
                    source = dict(spectrum["analysis"]["input"])
                    path = args.manifest.resolve().parent / source.pop("path")
                    load_spectrum(path, **source)
            report = {
                "status": "PASS",
                "spectra_count": len(spectra),
                "series_count": len({s["series"] for s in spectra}),
                "input_files_checked": args.check_files,
            }
        else:
            config = json.loads(args.config.read_text(encoding="utf-8-sig"))
            report = analyze(config, args.output, base=args.config.resolve().parent)
    except FileExistsError:
        parser.exit(
            2,
            text(
                "Output already exists; select a new directory.\n",
                "输出目录已存在，请换用新的目录。\n",
            ),
        )
    except (ValueError, KeyError, TypeError, OSError, IndexError) as exc:
        parser.exit(
            2,
            text(
                f"Operation could not run: {exc}\n",
                f"无法执行：{exc}\n请检查清单、列名、参数和输入文件；详见中文使用指南。\n",
            ),
        )
    except KeyboardInterrupt:
        parser.exit(
            130,
            text(
                "Interrupted. Batch progress is preserved in run.json and spectra.jsonl.\n",
                "操作已中断。批任务进度保存在 run.json 和 spectra.jsonl 中。\n",
            ),
        )
    if args.operation in {"batch", "demo-batch"}:
        summary = {
            k: report[k]
            for k in ("status", "state", "total", "fit_counts", "tracking_status", "tracks")
        }
    elif args.operation == "track":
        summary = {
            "status": report["status"],
            "tracks": len(report["tracks"]),
            "event_counts": report["event_counts"],
        }
    elif args.operation == "validate":
        summary = report
    elif args.operation == "stability":
        summary = {
            k: report[k]
            for k in (
                "status",
                "spectra",
                "variant_count",
                "stability_counts",
                "physical_peak_count_accuracy",
            )
        }
    else:
        summary = {k: report[k] for k in ("status", "data_kind", "peaks", "diagnostics")}
    if options.json:
        print(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False))
    else:
        print(text("Result", "结果") + ": " + summary["status"])
        if "total" in summary:
            print(
                text("Spectra", "光谱数")
                + f": {summary['total']}; "
                + text("Fit outcomes", "拟合状态统计")
                + f": {summary['fit_counts']}"
            )
            print(
                text("Tracking", "峰追踪")
                + f": {summary['tracking_status']}; "
                + text("Tracks", "轨迹数")
                + f": {summary['tracks']}"
            )
        elif "spectra_count" in summary:
            print(
                text("Validated spectra", "已校验光谱")
                + f": {summary['spectra_count']}; "
                + text("Independent series", "独立序列")
                + f": {summary['series_count']}"
            )
        elif "stability_counts" in summary:
            print(
                text("Spectra", "光谱数")
                + f": {summary['spectra']}; "
                + text("Variants", "配置变体")
                + f": {summary['variant_count']}"
            )
            print(text("Count stability", "峰数稳定性") + ": " + str(summary["stability_counts"]))
        if hasattr(args, "output"):
            print(text("Output", "结果目录") + ": " + str(args.output))
            if args.operation in {"batch", "demo-batch", "stability"}:
                directory = (
                    args.output / "analysis" if args.operation == "demo-batch" else args.output
                )
                print(
                    text("Open result overview", "打开结果总览")
                    + ": "
                    + str(directory / "report.html")
                )
    return 0 if report["status"] in {"PASS", "NOT_RUN"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
