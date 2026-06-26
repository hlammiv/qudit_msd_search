"""Loader for the paper's published search codes -- the correctness ORACLE.

Reads qmsd/data/puncture_locations.json (the machine-validated puncture-column
sets for the 10 Table-3 codes of arXiv:2510.10852 that have more than two punctures)
and exposes them as immutable OracleCode records. These are the deterministic
ground truth the discovery pipeline must reproduce exactly: given p, m and the
exact puncture columns, rebuilding the triorthogonal code must yield the listed
[[n,k,d]] (and, ideally, A_d). See tests/test_oracle.py.

The 1- and 2-puncture Table-3 codes are NOT here (their distance is
location-independent by Lemma 2, so they need no explicit column list); they live
in tests/ground_truth.py instead.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

_JSON = Path(__file__).resolve().parent / "data" / "puncture_locations.json"


@dataclass(frozen=True)
class OracleCode:
    label: str                              # e.g. "[[519,106,5]]_5"
    p: int                                  # prime qudit dimension
    n: int                                  # surviving columns  (= p**m - k)
    k: int                                  # logical qudits     (= #punctures)
    d: int                                  # distance
    m: int                                  # RM variables       (p**m = n + k)
    p_to_m: int                             # p**m
    num_punctures: int                      # == k
    puncture_columns_1indexed: tuple[int, ...]

    @property
    def r_max(self) -> int:
        """Largest RM degree with 3r < m(p-1) (triorthogonality; NOTES Thm 1).

        Uses m(p-1), NOT the paper's misprinted p(m-1) (see NOTES.md, section 4).
        """
        return (self.m * (self.p - 1) - 1) // 3

    @property
    def r_tilde(self) -> int:
        """Dual degree controlling the quantum distance: m(p-1) - r_max - 1."""
        return self.m * (self.p - 1) - self.r_max - 1


@lru_cache(maxsize=1)
def load_oracle(path: str | Path | None = None) -> tuple[OracleCode, ...]:
    data = json.loads(Path(path or _JSON).read_text())
    codes = []
    for c in data["codes"]:
        oc = OracleCode(
            label=c["label"],
            p=c["p"],
            n=c["n"],
            k=c["k"],
            d=c["d"],
            m=c["m"],
            p_to_m=c["p_to_m"],
            num_punctures=c["num_punctures"],
            puncture_columns_1indexed=tuple(c["puncture_columns_1indexed"]),
        )
        # Defensive invariants (the JSON was validated at extraction time).
        assert oc.p_to_m == oc.p ** oc.m == oc.n + oc.k, oc.label
        assert oc.num_punctures == oc.k == len(oc.puncture_columns_1indexed), oc.label
        codes.append(oc)
    return tuple(codes)
