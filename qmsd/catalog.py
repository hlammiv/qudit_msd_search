"""Portable search exports and validated local explorer catalog imports."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Iterable


SEARCH_SCHEMA = "qmsd.search-results.v1"
CATALOG_SCHEMA = "qmsd.catalog-record.v1"
_ROOT = Path(__file__).resolve().parents[1]


def catalog_directory(override: str | Path | None = None) -> Path:
    """Return the local catalog directory, allowing an override for automation/tests."""
    if override is not None:
        return Path(override).expanduser().resolve()
    configured = os.environ.get("QMSD_CATALOG_DIR")
    return (Path(configured).expanduser().resolve() if configured
            else _ROOT / "qmsd" / "data" / "catalog")


def code_payload(code) -> dict:
    """Serialize a ``Code`` without weakening its certification metadata."""
    return {
        "label": code.label, "p": code.p, "m": code.m, "r": code.r, "w": code.w,
        "n": code.n, "k": code.k, "d": code.d, "A_d": code.A_d,
        "gamma": code.gamma, "puncture_columns_1indexed": (
            list(code.puncture_columns) if code.puncture_columns is not None else None),
        "full_rank": code.full_rank, "distance_certified": bool(code.d_certified),
        "A_d_exact": code.A_d_exact,
    }


def write_search_export(path: str | Path, codes: Iterable, search: dict) -> Path:
    """Write a versioned search-result bundle and return its resolved path."""
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {"schema": SEARCH_SCHEMA, "search": search,
               "codes": [code_payload(code) for code in codes]}
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return target


def _stable_id(code: dict) -> str:
    identity = {key: code.get(key) for key in (
        "p", "m", "r", "w", "n", "k", "puncture_columns_1indexed")}
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True,
        separators=(",", ":")).encode()).hexdigest()[:12]
    return f"imported-p{code['p']}-n{code['n']}-k{code['k']}-{digest}"


def _validate_code(code: dict, requested_status: str | None) -> dict:
    required = ("p", "n", "k", "d", "distance_certified")
    missing = [key for key in required if key not in code]
    if missing:
        raise ValueError(f"code record is missing: {', '.join(missing)}")
    for key in ("p", "n", "k"):
        if not isinstance(code[key], int) or code[key] <= 0:
            raise ValueError(f"{key} must be a positive integer")
    columns = code.get("puncture_columns_1indexed")
    if columns is not None:
        if not isinstance(columns, list) or any(not isinstance(c, int) or c <= 0 for c in columns):
            raise ValueError("puncture_columns_1indexed must contain positive integers")
        if len(columns) != len(set(columns)):
            raise ValueError("puncture columns contain duplicates")
        if code.get("m") is None or len(columns) != code["k"]:
            raise ValueError("explicit records require m and exactly k puncture columns")
        if any(c > code["p"] ** code["m"] for c in columns):
            raise ValueError("puncture column exceeds p**m")
    inferred = ("confirmed" if code.get("distance_certified") and code.get("d") is not None
                and code.get("full_rank") is not False else "candidate")
    status = requested_status or inferred
    if status == "confirmed" and inferred != "confirmed":
        raise ValueError("cannot import an uncertified or rank-deficient result as confirmed")
    return {**code, "status": status}


def import_search_export(source: str | Path, destination: str | Path | None = None,
                         status: str | None = None) -> dict:
    """Validate and import all codes in a search bundle into the explorer catalog."""
    source_path = Path(source).expanduser().resolve()
    payload = json.loads(source_path.read_text())
    if payload.get("schema") != SEARCH_SCHEMA or not isinstance(payload.get("codes"), list):
        raise ValueError(f"expected a {SEARCH_SCHEMA} bundle")
    target_dir = catalog_directory(destination)
    target_dir.mkdir(parents=True, exist_ok=True)
    imported, unchanged = [], []
    provenance = payload.get("search", {})
    for raw in payload["codes"]:
        code = _validate_code(raw, status)
        artifact_id = _stable_id(code)
        family = "imported explicit search" if code.get("puncture_columns_1indexed") else (
            "imported analytic Manhattan")
        record = {
            "schema": CATALOG_SCHEMA, "artifact_id": artifact_id,
            "label": code.get("label") or f"[[{code['n']},{code['k']},{code['d']}]]_{code['p']}",
            "p": code["p"], "m": code.get("m"), "r": code.get("r"), "w": code.get("w"),
            "n": code["n"], "k": code["k"], "d": code.get("d"), "A_d": code.get("A_d"),
            "gamma_value": code.get("gamma"),
            "puncture_columns": code.get("puncture_columns_1indexed"),
            "family": family, "provenance": "qmsd CLI search import",
            "status": code["status"],
            "distance_evidence": ("Certified by the qmsd search engine" if code.get("distance_certified")
                                  else "Not certified by the qmsd search engine"),
            "Ad_evidence": ("Exact during qmsd search" if code.get("A_d_exact")
                            else "Search value; exactness not established" if code.get("A_d") is not None
                            else "Not computed"),
            "source": source_path.name,
            "note": "Imported search parameters: " + json.dumps(provenance, sort_keys=True),
        }
        target = target_dir / f"{artifact_id}.json"
        encoded = json.dumps(record, indent=2, sort_keys=True) + "\n"
        if target.exists():
            existing = json.loads(target.read_text())
            identity_keys = ("p", "m", "r", "w", "n", "k", "d", "puncture_columns")
            if any(existing.get(key) != record.get(key) for key in identity_keys):
                raise ValueError(f"artifact ID collision at {target}")
            unchanged.append(artifact_id)
        else:
            target.write_text(encoded)
            imported.append(artifact_id)
    return {"directory": str(target_dir), "imported": imported, "unchanged": unchanged}
