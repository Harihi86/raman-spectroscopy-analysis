"""Strict, versioned manifests. Configuration is data, never executable Python."""

import copy
import json
from importlib.resources import files

from jsonschema import Draft202012Validator

from .tracking import TrackingSettings


def read_schema(name="batch-manifest.schema.json") -> dict:
    return json.loads(files("raman_analysis").joinpath("schemas", name).read_text(encoding="utf-8"))


def validate_json(data, schema):
    json.dumps(data, allow_nan=False)
    errors = sorted(Draft202012Validator(schema).iter_errors(data), key=lambda e: str(e.json_path))
    if errors:
        details = "; ".join(f"{e.json_path}: {e.message}" for e in errors[:4])
        raise ValueError(f"Configuration validation failed: {details}")


def merge_settings(defaults, overrides):
    result = copy.deepcopy(defaults)
    for name, value in overrides.items():
        if isinstance(value, dict) and isinstance(result.get(name), dict):
            if name == "baseline" and value.get("method") != result[name].get("method"):
                result[name] = copy.deepcopy(value)
            else:
                result[name] = merge_settings(result[name], value)
        else:
            result[name] = copy.deepcopy(value)
    return result


def validate_manifest(manifest) -> list[dict]:
    schema = read_schema()
    validate_json(manifest, schema)
    defaults = manifest.get("defaults", {})
    if "path" in defaults.get("input", {}):
        raise ValueError("Specify input.path for each spectrum, not in defaults")
    normalized, ids = [], set()
    reserved = {
        "CON",
        "PRN",
        "AUX",
        "NUL",
        *[f"COM{i}" for i in range(1, 10)],
        *[f"LPT{i}" for i in range(1, 10)],
    }
    for sequence, entry in enumerate(manifest["spectra"]):
        spectrum_id = entry["id"]
        if spectrum_id.casefold() in ids or spectrum_id.upper() in reserved:
            raise ValueError("Spectrum IDs must be unique ignoring case and safe on Windows")
        ids.add(spectrum_id.casefold())
        analysis = merge_settings(defaults, entry.get("analysis", {}))
        analysis["input"] = merge_settings(analysis.get("input", {}), entry["input"])
        analysis.setdefault("data_kind", manifest.get("data_kind", "unspecified"))
        analysis.setdefault("baseline", {"method": "none"})
        analysis["metadata"] = merge_settings(
            analysis.get("metadata", {}), entry.get("metadata", {})
        )
        if "fit" not in analysis or "peaks" not in analysis["fit"]:
            raise ValueError(f"Spectrum {spectrum_id}: provide fit.peaks in defaults or analysis")
        resolved_schema = {**schema, "$ref": "#/$defs/analysis"}
        # Remove root-only constraints; retain local references for validation.
        resolved_schema = {
            k: v for k, v in resolved_schema.items() if k in {"$schema", "$defs", "$ref"}
        }
        validate_json(analysis, resolved_schema)
        normalized.append(
            {
                "id": spectrum_id,
                "series": entry.get("series", "default"),
                "label": entry.get("label", spectrum_id),
                "sequence": sequence,
                "pressure_gpa": entry.get("pressure_gpa"),
                "analysis": analysis,
            }
        )
    TrackingSettings(**manifest.get("tracking", {}))
    return normalized
