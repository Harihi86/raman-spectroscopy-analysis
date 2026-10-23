import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .fitting import gaussian
from .pipeline import analyze


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


def main():
    parser = argparse.ArgumentParser(
        description="Reproducible tabular Raman analysis with explicit assumptions"
    )
    sub = parser.add_subparsers(dest="operation", required=True)
    example = sub.add_parser("demo")
    example.add_argument("--output", required=True, type=Path)
    run = sub.add_parser("run")
    run.add_argument("config", type=Path)
    run.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        if args.operation == "demo":
            report = demo(args.output)
        else:
            config = json.loads(args.config.read_text(encoding="utf-8"))
            report = analyze(config, args.output, base=args.config.resolve().parent)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.exit(2, f"Analysis could not run: {exc}\n")
    print(
        json.dumps(
            {
                "status": report["status"],
                "data_kind": report["data_kind"],
                "peaks": report["peaks"],
                "diagnostics": report["diagnostics"],
            },
            indent=2,
        )
    )
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
