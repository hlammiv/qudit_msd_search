"""Portable search exports and the validated local explorer catalog.

Catalog records are JSON files in ``qmsd/data/catalog`` (schema ``qmsd.catalog-record.v2``, see
``CATALOG.md``).  Besides the v1 summary fields (p, n, k, d, A_d, status, free-text evidence) a v2
record carries four structured blocks so that every code can be *rebuilt* and its claims *checked*
by a program:

* ``construction``  how to rebuild G0/H1: ``punctured_rm`` (puncture columns of RM(r,m)),
                    ``reed_solomon`` (m = 1), ``weighted_rm`` (charge polynomial + outputs in X),
                    ``manhattan`` (weight-w puncturing), or ``none``.
* ``certificate``   lower/upper distance bounds, the method that established each, an explicit
                    minimum-weight word when there is one, and the stages/machine/date.  The
                    evidence status is *derived* from it (``status_from_certificate``).
* ``performance``   which performance data exist: exact curve, rigorous bounds, leading-order
                    census, acceptance Monte Carlo (pointers, not copies).
* ``geometry``      optional search-family descriptors (residual checks, class-block rank, plane
                    distance, maximum outputs per line).

v1 files are upgraded in memory on load (``upgrade_record``) and on disk by ``qmsd catalog migrate``.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
from pathlib import Path
from typing import Iterable


SEARCH_SCHEMA = "qmsd.search-results.v1"
CATALOG_SCHEMA = "qmsd.catalog-record.v2"
CATALOG_SCHEMA_V1 = "qmsd.catalog-record.v1"
ACCEPTED_SCHEMAS = (CATALOG_SCHEMA, CATALOG_SCHEMA_V1)
CONSTRUCTION_FAMILIES = ("punctured_rm", "reed_solomon", "weighted_rm", "manhattan", "none")
_ROOT = Path(__file__).resolve().parents[1]


def catalog_directory(override: str | Path | None = None) -> Path:
    """Return the local catalog directory, allowing an override for automation/tests."""
    if override is not None:
        return Path(override).expanduser().resolve()
    configured = os.environ.get("QMSD_CATALOG_DIR")
    return (Path(configured).expanduser().resolve() if configured
            else _ROOT / "qmsd" / "data" / "catalog")


# ---------------------------------------------------------------------------------------------
# structured blocks
# ---------------------------------------------------------------------------------------------

def construction_block(family: str, p: int, m: int | None = None, r: int | None = None, *,
                       puncture_columns_1indexed: list[int] | None = None, w: int | None = None,
                       charge: dict | None = None, outputs_in_X: list[int] | None = None,
                       class_unit: int | None = None) -> dict:
    """A rebuild recipe.  ``charge`` for weighted_rm is {"kind": "quadric", "diag": [...],
    "lin": [...], "const": c}; ``outputs_in_X`` are 0-based indices into the support X."""
    if family not in CONSTRUCTION_FAMILIES:
        raise ValueError(f"unknown construction family {family!r}")
    block = {"family": family, "p": int(p), "m": m, "r": r}
    if family in ("punctured_rm", "reed_solomon"):
        block["puncture_columns_1indexed"] = list(puncture_columns_1indexed or [])
    elif family == "manhattan":
        block["w"] = w
    elif family == "weighted_rm":
        block["charge"] = charge
        block["outputs_in_X"] = list(outputs_in_X or [])
        block["class_unit"] = class_unit if class_unit is not None else 1
    return block


def certificate_block(lower_bound: int | None, upper_bound: int | None, *,
                      lower_method: str = "", upper_method: str = "",
                      witness: dict | None = None, stages: list | None = None,
                      machine: str = "", date: str | None = None, refuted: bool = False) -> dict:
    """Distance certificate.  ``witness`` = {"columns": [...G0 columns...], "coefficients": [...],
    "weight": w} is an explicit code word (proves d <= w).  ``exact`` iff the bounds meet."""
    exact = lower_bound is not None and upper_bound is not None and lower_bound == upper_bound
    return {"lower_bound": lower_bound, "upper_bound": upper_bound, "exact": bool(exact),
            "lower_method": lower_method, "upper_method": upper_method,
            "witness": witness, "stages": stages or [], "machine": machine,
            "date": date or _dt.date.today().isoformat(), "refuted": bool(refuted)}


def status_from_certificate(cert: dict | None, fallback: str = "candidate") -> str:
    """confirmed = exact distance; partial = a rigorous lower bound or a closed-form distance with
    the exact value not independently certified; candidate = no rigorous bound; refuted."""
    if not cert:
        return fallback
    if cert.get("refuted"):
        return "refuted"
    if cert.get("exact"):
        return "confirmed"
    if cert.get("lower_bound") is not None:
        return "partial"
    return "candidate"


def performance_block(exact: bool = False, bounds: bool = False, census: dict | None = None,
                      acceptance_mc: dict | None = None, notes: str = "") -> dict:
    """What performance data exist for the record (pointers; the numbers live in
    qmsd/data/performance and runs/)."""
    return {"exact_curve": bool(exact), "rigorous_bounds": bool(bounds), "census": census,
            "acceptance_mc": acceptance_mc, "notes": notes}


def _infer_construction(code: dict) -> dict | None:
    cols = code.get("puncture_columns_1indexed")
    if code.get("construction"):
        return code["construction"]
    if cols:
        fam = "reed_solomon" if code.get("m") == 1 else "punctured_rm"
        return construction_block(fam, code["p"], code.get("m"), code.get("r"), puncture_columns_1indexed=cols)
    if code.get("w") is not None:
        return construction_block("manhattan", code["p"], code.get("m"), code.get("r"), w=code["w"])
    return None


def _infer_certificate(code: dict) -> dict | None:
    if code.get("certificate"):
        return code["certificate"]
    d = code.get("d")
    if d is None:
        return None
    if code.get("distance_certified"):
        return certificate_block(d, d, lower_method="qmsd search engine (MITM min_dependent_columns)",
                                 upper_method="qmsd search engine (explicit dependent set)")
    return None


def upgrade_record(data: dict) -> dict:
    """Upgrade a v1 catalog record dict to v2 in memory (idempotent on v2)."""
    if data.get("schema") == CATALOG_SCHEMA and all(k in data for k in ("construction", "certificate")):
        return data
    out = dict(data)
    out["schema"] = CATALOG_SCHEMA
    if not out.get("construction"):
        out["construction"] = _infer_construction({
            "p": out["p"], "m": out.get("m"), "r": out.get("r"), "w": out.get("w"),
            "puncture_columns_1indexed": out.get("puncture_columns")}) or {"family": "none", "p": out["p"], "m": out.get("m"), "r": out.get("r")}
    if "certificate" not in out:
        d = out.get("d"); st = out.get("status")
        if st == "refuted":
            out["certificate"] = certificate_block(None, None, refuted=True, lower_method=out.get("distance_evidence", ""))
        elif st == "confirmed" and d is not None:
            out["certificate"] = certificate_block(d, d, lower_method=out.get("distance_evidence", ""), upper_method=out.get("distance_evidence", ""))
        elif st == "partial" and d is not None:
            out["certificate"] = certificate_block(d, d, lower_method=out.get("distance_evidence", ""), upper_method=out.get("distance_evidence", ""))
            out["certificate"]["exact"] = False  # closed-form / not independently certified
        else:
            out["certificate"] = None
    out.setdefault("performance", performance_block())
    out.setdefault("geometry", None)
    return out


# ---------------------------------------------------------------------------------------------
# search exports and imports
# ---------------------------------------------------------------------------------------------

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
    cons = code.get("construction") or {}
    if cons.get("family") == "weighted_rm":
        identity["weighted"] = {"charge": cons.get("charge"), "outputs_in_X": cons.get("outputs_in_X")}
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
    cons = code.get("construction")
    if cons is not None and cons.get("family") not in CONSTRUCTION_FAMILIES:
        raise ValueError(f"unknown construction family {cons.get('family')!r}")
    if cons is not None and cons.get("family") == "weighted_rm":
        if not cons.get("charge") or len(cons.get("outputs_in_X") or []) != code["k"]:
            raise ValueError("weighted_rm construction needs a charge and exactly k outputs_in_X")
    cert = _infer_certificate(code)
    inferred = status_from_certificate(cert)
    if inferred == "confirmed" and code.get("full_rank") is False:
        inferred = "candidate"
    status = requested_status or inferred
    if status == "confirmed" and inferred != "confirmed":
        raise ValueError("cannot import an uncertified or rank-deficient result as confirmed")
    return {**code, "status": status, "certificate": cert, "construction": _infer_construction(code)}


def import_search_export(source: str | Path, destination: str | Path | None = None,
                         status: str | None = None) -> dict:
    """Validate and import all codes in a search bundle into the explorer catalog (v2 records).

    Bundle codes may carry optional ``construction``, ``certificate``, ``performance`` and
    ``geometry`` blocks (see CATALOG.md); otherwise they are inferred from the v1 fields."""
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
        cons = code.get("construction") or {}
        family = {"punctured_rm": "imported explicit search", "reed_solomon": "imported Reed-Solomon (m = 1)",
                  "weighted_rm": "imported weighted support", "manhattan": "imported analytic Manhattan"}.get(cons.get("family"), "imported search")
        cert = code.get("certificate")
        record = {
            "schema": CATALOG_SCHEMA, "artifact_id": artifact_id,
            "label": code.get("label") or f"[[{code['n']},{code['k']},{code['d']}]]_{code['p']}",
            "p": code["p"], "m": code.get("m"), "r": code.get("r"), "w": code.get("w"),
            "n": code["n"], "k": code["k"], "d": code.get("d"), "A_d": code.get("A_d"),
            "gamma_value": code.get("gamma"),
            "puncture_columns": code.get("puncture_columns_1indexed"),
            "family": family, "provenance": "qmsd CLI search import",
            "status": code["status"],
            "distance_evidence": (code.get("distance_evidence") or (cert and cert.get("lower_method")) or
                                  ("Certified by the qmsd search engine" if code.get("distance_certified")
                                   else "Not certified by the qmsd search engine")),
            "Ad_evidence": ("Exact during qmsd search" if code.get("A_d_exact")
                            else "Search value; exactness not established" if code.get("A_d") is not None
                            else "Not computed"),
            "source": source_path.name,
            "note": code.get("note") or ("Imported search parameters: " + json.dumps(provenance, sort_keys=True)),
            "construction": cons or {"family": "none", "p": code["p"], "m": code.get("m"), "r": code.get("r")},
            "certificate": cert,
            "performance": code.get("performance") or performance_block(),
            "geometry": code.get("geometry"),
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


def migrate_catalog(directory: str | Path | None = None) -> dict:
    """Rewrite every v1 record in the catalog directory as v2 (in place; idempotent)."""
    target_dir = catalog_directory(directory)
    upgraded, kept = [], []
    for path in sorted(target_dir.glob("*.json")):
        data = json.loads(path.read_text())
        if data.get("schema") not in ACCEPTED_SCHEMAS:
            raise ValueError(f"unsupported catalog schema in {path.name!r}")
        new = upgrade_record(data)
        if new != data:
            path.write_text(json.dumps(new, indent=2, sort_keys=True) + "\n")
            upgraded.append(path.name)
        else:
            kept.append(path.name)
    return {"directory": str(target_dir), "upgraded": upgraded, "unchanged": kept}
