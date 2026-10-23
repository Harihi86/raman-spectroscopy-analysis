#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

from photoactivity_processing import run_photoactivity_profile
from photoactivity_profiles import DATASETS, RUN_PROFILES, deep_copy_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the photoactivity time-series pipeline.")
    parser.add_argument("--profile", type=str, default=None, help="Profile key from RUN_PROFILES.")
    parser.add_argument("--list-profiles", action="store_true", help="List all available profiles and exit.")
    parser.add_argument("--excel", type=str, default=None, help="Override excel_path for single-dataset profiles only.")
    parser.add_argument("--output", type=str, default=None, help="Override output_dir.")
    return parser.parse_args()


def write_cli_log(output_dir: Path, run_profile: Dict[str, Any], cli_args: Dict[str, Any]) -> None:
    payload = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "profile": run_profile.get("name"),
        "output_dir": str(output_dir.resolve()),
        "cli_args": cli_args,
    }
    with open(output_dir / "cli_log.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)


def main() -> None:
    args = parse_args()

    if args.list_profiles:
        print("Available profiles:")
        for key, prof in RUN_PROFILES.items():
            print(f"  {key}: {prof.get('description', '')}")
        return

    if not args.profile:
        raise SystemExit("Please provide --profile, or use --list-profiles first.")
    if args.profile not in RUN_PROFILES:
        raise SystemExit(f"Unknown profile: {args.profile}")

    run_profile: Dict[str, Any] = deep_copy_config(RUN_PROFILES[args.profile])
    dataset_keys: List[str] = list(run_profile["dataset_keys"])
    dataset_configs: List[Dict[str, Any]] = [deep_copy_config(DATASETS[k]) for k in dataset_keys]

    if args.excel is not None:
        if len(dataset_configs) != 1:
            raise SystemExit("--excel override is only supported for single-dataset profiles.")
        dataset_configs[0]["excel_path"] = args.excel

    if args.output is not None:
        run_profile["output_dir"] = args.output

    print(f"Running profile: {run_profile['name']}")
    print(f"Output: {run_profile['output_dir']}")
    for cfg in dataset_configs:
        print(f"  - {cfg['dataset_name']}: {cfg['excel_path']}")

    result = run_photoactivity_profile(run_profile, dataset_configs)
    output_dir = Path(result["output_dir"])
    write_cli_log(output_dir, run_profile, vars(args))

    print("Done.")
    print(f"Datasets processed: {result['n_datasets']}")
    print(f"Output directory: {result['output_dir']}")
    print(f"Summary workbook: {result['summary_xlsx']}")
    print(f"Summary csv: {result['summary_csv']}")


if __name__ == "__main__":
    main()
