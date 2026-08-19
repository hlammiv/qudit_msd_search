"""Curated, read-only research results for the Streamlit explorer.

The files in the repository predate a common result schema.  This module is the
deliberately small normalization boundary: scientific claims are registered
explicitly, so an old JSON file can never silently become a confirmed result.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from functools import lru_cache
import numpy as np

from pathlib import Path
from typing import Iterable, Literal

from .codes import code_from_manhattan
from .catalog import CATALOG_SCHEMA, catalog_directory
from .distillation import cost, delta_out_avg, nbar_T
from .field import GFp
from .oracle import load_oracle
from .puncture import column_to_point
from .reedmuller import points
from .triorthogonal import build_triorthogonal_code

EvidenceStatus = Literal["confirmed", "partial", "candidate", "refuted"]

_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class ResultRecord:
    artifact_id: str
    label: str
    p: int
    n: int
    k: int
    d: int | None
    m: int | None = None
    r: int | None = None
    w: int | None = None
    A_d: int | None = None
    gamma_value: float | None = None
    puncture_columns: tuple[int, ...] | None = None
    family: str = "explicit puncture"
    provenance: str = "repository result"
    status: EvidenceStatus = "confirmed"
    distance_evidence: str = "Not recorded"
    Ad_evidence: str = "Not recorded"
    source: str = ""
    note: str = ""

    @property
    def gamma(self) -> float | None:
        if self.gamma_value is not None:
            return self.gamma_value
        if self.d is None or self.d <= 1 or self.n <= 0 or self.k <= 0:
            return None
        return math.log(self.n / self.k) / math.log(self.d)

    @property
    def rate(self) -> float:
        return self.k / self.n

    @property
    def sublogarithmic(self) -> bool:
        return self.gamma is not None and self.gamma < 1

    @property
    def distance_certified(self) -> bool:
        return self.status == "confirmed" and self.d is not None


# Paper Table 3 counts, keyed independently of puncture-set identity.
_TABLE3 = (
    (3, 4, 80, 1, 5, 2080), (3, 4, 79, 2, 4, 130),
    (3, 4, 72, 9, 3, 648), (3, 5, 230, 13, 6, 572),
    (3, 5, 215, 28, 5, 1104), (3, 5, 206, 37, 4, 880),
    (3, 5, 200, 43, 3, 1700), (3, 6, 690, 39, 5, 1128),
    (3, 6, 667, 62, 4, 3972), (5, 2, 24, 1, 3, 96),
    (5, 2, 20, 5, 2, 760), (5, 3, 124, 1, 4, 124),
    (5, 3, 112, 13, 3, 512), (5, 4, 519, 106, 5, 2180),
)

_TABLE2 = (
    (2, 58, 14), (3, 32, 16), (5, 16, 16), (7, 13, 20),
    (11, 7, 19), (13, 7, 23), (17, 4, 17), (19, 4, 19), (23, 1, 5),
)


def _paper_records() -> list[ResultRecord]:
    oracle = {(o.p, o.m, o.n, o.k, o.d): o for o in load_oracle()}
    rows: list[ResultRecord] = []
    for p, m, n, k, d, ad in _TABLE3:
        oc = oracle.get((p, m, n, k, d))
        cols = oc.puncture_columns_1indexed if oc else None
        rows.append(ResultRecord(
            artifact_id=f"paper-table3-p{p}-m{m}-n{n}-k{k}",
            label=f"[[{n},{k},{d}]]_{p}", p=p, m=m, r=(m * (p - 1) - 1) // 3, n=n, k=k, d=d, A_d=ad,
            puncture_columns=cols, family="paper explicit search",
            provenance="arXiv:2510.10852 Table 3",
            distance_evidence="Exact; reproduced by the qmsd distance certifier",
            Ad_evidence="Paper Table 3 value; reproduced where computationally feasible",
            source="qmsd/data/puncture_locations.json" if cols else "paper Table 3 ground truth",
        ))
    return rows


def _manhattan_records() -> list[ResultRecord]:
    rows = []
    for p, m, w in _TABLE2:
        c = code_from_manhattan(p, m, w)
        rows.append(ResultRecord(
            artifact_id=f"manhattan-table2-p{p}-m{m}-w{w}", label=c.label,
            p=p, m=m, r=c.r, w=w, n=c.n, k=c.k, d=c.d,
            family="analytic Manhattan", provenance="arXiv:2510.10852 Table 2",
            distance_evidence="Exact closed-form Manhattan-family distance",
            Ad_evidence="Not computed", source="qmsd.codes.code_from_manhattan",
        ))
    return rows


def _read_json(name: str) -> dict:
    return json.loads((_ROOT / name).read_text())


def _artifact_records() -> list[ResultRecord]:
    rows: list[ResultRecord] = []

    def add(name: str, artifact_id: str, *, status: EvidenceStatus = "confirmed",
            family: str = "new explicit search", label_key: str = "code",
            ad_key: str | None = None, distance_evidence: str = "Exact, repository-certified",
            ad_evidence: str = "Not computed", note_key: str = "note") -> None:
        data = _read_json(name)
        d = data.get("d")
        prime = int(data.get("p", 3))
        label = data.get(label_key) or data.get("code_label") or (
            f"[[{data['n']},{data['k']},{d if d is not None else '?'}]]_{prime}")
        cols = data.get("puncture_columns_1indexed", data.get("cols"))
        ad = data.get(ad_key) if ad_key else data.get("A_d")
        rows.append(ResultRecord(
            artifact_id=artifact_id, label=label, p=prime,
            m=data.get("m", 7 if name == "cap_qutrit_code.json" else None), r=data.get("r"),
            n=int(data["n"]), k=int(data["k"]), d=d, A_d=ad,
            gamma_value=data.get("gamma", data.get("gamma_upper")),
            puncture_columns=tuple(cols) if cols else None, family=family,
            provenance="This repository", status=status,
            distance_evidence=distance_evidence, Ad_evidence=ad_evidence,
            source=name, note=str(data.get(note_key, "")),
        ))

    add("p17_d6_code.json", "new-p17-d6", ad_key="A_6",
        distance_evidence="d=6 exact by MITM lower bound and structured line witness",
        ad_evidence="A_6 exact by balanced streamed enumeration", note_key="insight")
    add("p19_lock.json", "new-p19-d5",
        distance_evidence="d=5 exact by MITM lower bound and structured line witness")
    add("qutrit_Ad572.json", "optimized-qutrit-ad572",
        family="A_d-optimized puncture set", ad_evidence="A_d exact by two independent methods")
    add("p5_Ad_code.json", "optimized-p5-m3-ad396",
        family="A_d-optimized puncture set", ad_evidence="A_d exact by MacWilliams enumeration")
    add("flagship_p5_d5_Ad.json", "optimized-p5-flagship-ad1904", ad_key="A_5",
        family="A_d-optimized puncture set",
        distance_evidence="d=5 exact by meet-in-the-middle certification",
        ad_evidence="A_5 exact by streamed weight-5 enumeration")
    add("cap_qutrit_code.json", "refuted-qutrit-cap", status="refuted",
        family="refuted cap construction", label_key="code_label",
        distance_evidence="Refuted: later validation found true distance d=1",
        ad_evidence="Not applicable",
        note_key="note")
    add("p17_d7_candidate.json", "candidate-p17-d7", status="candidate",
        family="line-bound candidate", label_key="missing",
        distance_evidence="Only the line-derived bound was recorded; full-span d=6 obstruction found")
    add("p19_d7_candidate.json", "candidate-p19-d7", status="candidate",
        family="line-bound candidate", label_key="missing",
        distance_evidence="Only the line-derived bound was recorded; full-span d=5 obstruction found")
    return rows


def _catalog_fingerprint() -> tuple[tuple[str, int, int], ...]:
    directory = catalog_directory()
    if not directory.exists():
        return ()
    return tuple((path.name, path.stat().st_mtime_ns, path.stat().st_size)
                 for path in sorted(directory.glob("*.json")))


def _imported_records(directory: Path,
                      fingerprint: tuple[tuple[str, int, int], ...]) -> list[ResultRecord]:
    rows = []
    for name, _, _ in fingerprint:
        data = json.loads((directory / name).read_text())
        if data.get("schema") != CATALOG_SCHEMA:
            raise ValueError(f"unsupported catalog schema in {name!r}")
        rows.append(ResultRecord(
            artifact_id=data["artifact_id"], label=data["label"], p=int(data["p"]),
            m=data.get("m"), r=data.get("r"), w=data.get("w"), n=int(data["n"]),
            k=int(data["k"]), d=data.get("d"), A_d=data.get("A_d"),
            gamma_value=data.get("gamma_value"),
            puncture_columns=(tuple(data["puncture_columns"])
                              if data.get("puncture_columns") else None),
            family=data.get("family", "imported search"),
            provenance=data.get("provenance", "qmsd catalog import"),
            status=data.get("status", "candidate"),
            distance_evidence=data.get("distance_evidence", "Not recorded"),
            Ad_evidence=data.get("Ad_evidence", "Not recorded"),
            source=data.get("source", name), note=data.get("note", ""),
        ))
    return rows


@lru_cache(maxsize=16)
def _load_result_catalog(directory_name: str,
                         fingerprint: tuple[tuple[str, int, int], ...]) -> tuple[ResultRecord, ...]:
    """Return the curated catalog, with stable artifact IDs and explicit evidence states."""
    records = (_paper_records() + _manhattan_records() + _artifact_records()
               + _imported_records(Path(directory_name), fingerprint))
    ids = [r.artifact_id for r in records]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate result artifact_id in curated catalog")
    return tuple(records)


def load_result_catalog() -> tuple[ResultRecord, ...]:
    """Return built-in and imported records, refreshing when catalog files change."""
    directory = catalog_directory()
    return _load_result_catalog(str(directory), _catalog_fingerprint())


def distillation_series(record: ResultRecord, delta_inputs: Iterable[float]) -> list[dict]:
    """Compute cheap, display-only distillation metrics without recertifying a code."""
    out = []
    for delta in delta_inputs:
        nbar = nbar_T(record.n, record.k, record.p, float(delta))
        row = {"delta_in": float(delta), "accepted_outputs": nbar,
               "cost": cost(record.n, nbar) if nbar > 0 else float("inf"),
               "delta_out": None}
        if (nbar > 0 and record.status == "confirmed" and record.d
                and record.A_d is not None):
            row["delta_out"] = delta_out_avg(
                record.n, record.k, record.d, record.A_d, record.p, float(delta))
        out.append(row)
    return out


def pareto_front(records: Iterable[ResultRecord], x: str, y: str) -> set[str]:
    """IDs minimizing both named numeric properties; uncertain claims are excluded."""
    candidates = []
    for record in records:
        if record.status != "confirmed":
            continue
        xv, yv = getattr(record, x), getattr(record, y)
        xv = xv() if callable(xv) else xv
        yv = yv() if callable(yv) else yv
        if xv is not None and yv is not None and math.isfinite(float(xv)) and math.isfinite(float(yv)):
            candidates.append((record, float(xv), float(yv)))
    return {a.artifact_id for a, ax, ay in candidates
            if not any((bx <= ax and by <= ay) and (bx < ax or by < ay)
                       for b, bx, by in candidates if b is not a)}


def puncture_points(record: ResultRecord) -> list[tuple[int, ...]]:
    if record.m is None or not record.puncture_columns:
        return []
    return [column_to_point(c, record.m, record.p) for c in record.puncture_columns]


@lru_cache(maxsize=32)
def affine_line_profile(record: ResultRecord, max_ambient: int = 700) -> dict | None:
    """Count selected punctures on every affine line in ``F_p^m``."""
    if (record.m is None or not record.puncture_columns
            or record.p ** record.m > max_ambient):
        return None
    p, m = record.p, record.m
    ambient = np.asarray(points(m, p), dtype=np.int64)
    punctured = set(record.puncture_columns)
    directions = []
    for vector in ambient[1:]:
        first = int(np.flatnonzero(vector)[0])
        if vector[first] == 1 and np.all(vector[:first] == 0):
            directions.append(vector)
    lines: set[tuple[int, ...]] = set()
    for direction in directions:
        for origin in ambient:
            line = []
            for scalar in range(p):
                point = tuple(((origin + scalar * direction) % p).tolist())
                column = 1 + sum(value * p ** j for j, value in enumerate(point))
                line.append(column)
            lines.add(tuple(sorted(line)))
    occupancies = np.asarray(
        [sum(column in punctured for column in line) for line in lines], dtype=int)
    max_occupancy = int(occupancies.max(initial=0))
    witness = next((line for line in sorted(lines)
                    if sum(column in punctured for column in line) == max_occupancy), ())
    values, counts = np.unique(occupancies, return_counts=True)
    return {
        "histogram": tuple((int(value), int(count)) for value, count in zip(values, counts)),
        "max_occupancy": max_occupancy,
        "max_line_columns": tuple(witness),
        "line_count": len(lines),
    }


def record_dict(record: ResultRecord) -> dict:
    data = asdict(record)
    data["gamma"] = record.gamma
    data["rate"] = record.rate
    data["sublogarithmic"] = record.sublogarithmic
    return data


@lru_cache(maxsize=32)
def code_structure(record: ResultRecord, max_ambient: int = 2500) -> dict | None:
    """Build exact CSS incidence data and logical-X quotient representatives.

    The result describes the algebraic code, not a hardware coupling graph.  It is
    intentionally gated because analytic codes can have astronomical ambient spaces.
    """
    if (record.puncture_columns is None or record.m is None or record.r is None
            or record.p ** record.m > max_ambient):
        return None
    built = build_triorthogonal_code(
        record.p, record.m, record.r, record.puncture_columns)
    x_stab = np.asarray(built["X_stab"], dtype=np.int64) % record.p
    z_stab = np.asarray(built["Z_stab"], dtype=np.int64) % record.p
    gp = np.asarray(built["Gp"], dtype=np.int64) % record.p

    # Select a basis of Gp/G0.  Each added row is an X-logical representative;
    # its support/coefficients show how one logical qudit is carried physically.
    gf = GFp(record.p)
    basis = [row.copy() for row in x_stab]
    rank = int(np.linalg.matrix_rank(gf(x_stab))) if len(x_stab) else 0
    logical_x = []
    for row in gp:
        trial = np.vstack([*basis, row]) if basis else row.reshape(1, -1)
        new_rank = int(np.linalg.matrix_rank(gf(trial)))
        if new_rank > rank:
            logical_x.append(row.copy())
            basis.append(row.copy())
            rank = new_rank
        if len(logical_x) == record.k:
            break

    punctured = set(record.puncture_columns)
    physical_columns = tuple(
        c for c in range(1, record.p ** record.m + 1) if c not in punctured)
    return {
        "X_stab": x_stab,
        "Z_stab": z_stab,
        "logical_X": np.asarray(logical_x, dtype=np.int64),
        "physical_columns": physical_columns,
        "full_rank": built["full_rank"],
    }
