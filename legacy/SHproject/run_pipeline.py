#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

from core_processing import deep_copy_config, process_profile, run_comparison_group
from profiles import PROFILES, COMPARISON_GROUPS


def load_pressures_from_text(path: str) -> List[float]:
    vals: List[float] = []
    with open(path, "r", encoding="utf-8") as f:
        text = f.read().strip()
    if not text:
        return vals
    if path.lower().endswith(".json"):
        arr = json.loads(text)
        return [float(v) for v in arr]
    for line in text.replace(",", " ").split():
        vals.append(float(line))
    return vals


def infer_default_output_dir(cfg: Dict[str, Any]) -> Path:
    group = str(cfg.get("output_group", "misc")).strip().replace("\\", "/")
    leaf = str(cfg.get("output_dir", cfg.get("name", "run_output"))).strip()
    return Path("outputs") / Path(group) / leaf


def sanitize_for_json(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [sanitize_for_json(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    return obj


def write_run_log(output_dir: Path, final_cfg: Dict[str, Any], cli_args: Dict[str, Any], result: Dict[str, Any]) -> None:
    payload = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "profile": final_cfg.get("name"),
        "excel_path": str(Path(final_cfg["excel_path"]).resolve()) if final_cfg.get("excel_path") else None,
        "output_dir": str(output_dir.resolve()),
        "n_spectra": result.get("n_spectra"),
        "cli_args": sanitize_for_json(cli_args),
        "final_config": sanitize_for_json(final_cfg),
    }
    with open(output_dir / "run_log.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)



def parse_excel_map(values: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in values:
        if '=' not in item:
            raise SystemExit(f'Invalid --excel-map entry: {item}. Expected PROFILE=path.xlsx')
        key, value = item.split('=', 1)
        key = key.strip()
        value = value.strip()
        if not key or not value:
            raise SystemExit(f'Invalid --excel-map entry: {item}. Expected PROFILE=path.xlsx')
        out[key] = value
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Raman pipeline with a selected profile or comparison group.")
    parser.add_argument("--profile", type=str, default=None, help="Profile key from PROFILES.")
    parser.add_argument("--compare-group", type=str, default=None, help="Comparison group key from COMPARISON_GROUPS.")
    parser.add_argument("--list-profiles", action="store_true", help="List all available profiles and exit.")
    parser.add_argument("--list-compare-groups", action="store_true", help="List all comparison groups and exit.")
    parser.add_argument("--excel", type=str, default=None, help="Override excel_path for a single profile run.")
    parser.add_argument("--excel-map", action="append", default=[], help="Comparison input mapping: PROFILE=path.xlsx . Repeat as needed.")
    parser.add_argument("--sheet", type=str, default=None, help="Override sheet name or sheet index.")
    parser.add_argument("--output", type=str, default=None, help="Override output_dir completely.")
    parser.add_argument("--compare-output", type=str, default=None, help="Override output directory for comparison group runs.")
    parser.add_argument("--pressures-file", type=str, default=None, help="Text or JSON file containing pressure values.")
    parser.add_argument("--set-fit-mode", type=str, default=None, help="Override fit mode: single | double | auto | none.")
    parser.add_argument("--no-progress", action="store_true", help="Disable terminal progress bars.")
    return parser.parse_args()



def main() -> None:
    args = parse_args()

    if args.list_profiles:
        print("Available profiles:")
        for key, prof in PROFILES.items():
            print(f"  {key}: {prof.get('description', '')}")
        if not args.list_compare_groups:
            return

    if args.list_compare_groups:
        print("Available comparison groups:")
        for key, grp in COMPARISON_GROUPS.items():
            print(f"  {key}: {grp.get('description', '')}")
        if not args.profile and not args.compare_group:
            return

    if args.compare_group:
        if args.compare_group not in COMPARISON_GROUPS:
            raise SystemExit(f"Unknown comparison group: {args.compare_group}")
        excel_map = parse_excel_map(args.excel_map or [])
        group_cfg = deep_copy_config(COMPARISON_GROUPS[args.compare_group])
        group_cfg['progress'] = not args.no_progress
        if args.compare_output is not None:
            outdir = Path(args.compare_output)
        elif args.output is not None:
            outdir = Path(args.output)
        else:
            slug = str(group_cfg.get('output_dir', args.compare_group.lower()))
            outdir = Path('outputs') / 'comparisons' / slug
        print(f"Running comparison group: {args.compare_group}")
        print(f"Output: {outdir}")
        result = run_comparison_group(group_cfg, PROFILES, excel_map, outdir)
        print("Done.")
        print(f"Comparison output directory: {result['output_dir']}")
        print(f"Summary workbook: {result['summary_path']}")
        return

    if not args.profile:
        raise SystemExit("Please provide --profile or --compare-group, or use listing flags first.")
    if args.profile not in PROFILES:
        raise SystemExit(f"Unknown profile: {args.profile}")

    cfg: Dict[str, Any] = deep_copy_config(PROFILES[args.profile])
    cfg['progress'] = not args.no_progress

    if args.excel is not None:
        cfg["excel_path"] = args.excel
    if args.sheet is not None:
        try:
            cfg["sheet_name"] = int(args.sheet)
        except ValueError:
            cfg["sheet_name"] = args.sheet
    if args.pressures_file is not None:
        cfg["pressure_list"] = load_pressures_from_text(args.pressures_file)
    if args.set_fit_mode is not None:
        cfg.setdefault("peak_tracking", {})["mode"] = args.set_fit_mode

    if args.output is not None:
        cfg["output_dir"] = args.output
    else:
        cfg["output_dir"] = str(infer_default_output_dir(cfg))

    print(f"Running profile: {cfg['name']}")
    print(f"Excel: {cfg['excel_path']}")
    print(f"Output: {cfg['output_dir']}")

    result = process_profile(cfg)
    outdir = Path(result["output_dir"])
    write_run_log(outdir, cfg, vars(args), result)

    print("Done.")
    print(f"Profile: {cfg['name']}")
    print(f"Read {result['n_spectra']} spectra")
    print(f"Output directory: {result['output_dir']}")
    print(f"Run log: {str((outdir / 'run_log.json').resolve())}")


if __name__ == "__main__":
    main()
