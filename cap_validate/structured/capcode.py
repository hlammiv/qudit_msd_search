"""Load the m=7 cap qutrit code and expose the geometry needed by Method 2.

The code under validation is

    built = build_triorthogonal_code(3, 7, 4, S, G=rm_generator(4,7,3))
    G0    = built["X_stab"]                       # 55 x 1968 over F_3

with  S = the 219 cap puncture columns (1-indexed).  Its quantum distance is

    d = min_{c in RM_3(9,7), c != 0}  |supp(c) \\ S|

because G0^perp = RM_3(9,7) punctured at S  (RM_3(4,7)^perp = RM_3(rtilde,7),
rtilde = m(p-1)-r-1 = 7*2-4-1 = 9).  Method 2 lower-bounds d by enumerating the
*structured* low-weight codewords of RM_3(9,7) (weight 18..45) and minimising
|supp(c) \\ S| over them.

This module provides the F_3^7 point lattice and the 219 cap points as integer
coordinate rows; everything downstream is pure geometry over those points.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

P = 3
M = 7
R = 4                     # RM_3(4,7) is the triorthogonal generator (3r=12 < m(p-1)=14)
N = P ** M               # 2187
RTILDE = M * (P - 1) - R - 1   # 9  -> dual code is RM_3(9,7)
D_RM_MIN = 18            # min weight of RM_3(9,7) (Schwartz-Zippel): 3^3 * (3-1)//3 = 18

REPO = Path(__file__).resolve().parents[2]
CAP_JSON = REPO / "cap_qutrit_code.json"


def points() -> np.ndarray:
    """All p^m points of F_3^7 as rows of base-3 digits (x_1 least significant).

    Row c (0-indexed) is the point whose 1-indexed column in the RM generator is c+1
    -- matches qmsd.reedmuller.points / the puncture-column convention.
    """
    pts = np.empty((N, M), dtype=np.int8)
    for j in range(M):
        pts[:, j] = (np.arange(N) // (P ** j)) % P
    return pts


def load_cap(path: Path | str | None = None) -> dict:
    """Load the cap S: returns dict with 1-indexed columns, 0-indexed indices, points."""
    path = Path(path) if path is not None else CAP_JSON
    d = json.loads(path.read_text())
    cols1 = list(d["puncture_columns_1indexed"])
    idx0 = np.asarray([c - 1 for c in cols1], dtype=np.int64)
    assert len(set(cols1)) == len(cols1), "duplicate cap columns"
    assert all(1 <= c <= N for c in cols1), "cap column out of range"
    pts = points()
    return {
        "cols1": cols1,
        "idx0": idx0,                       # 0-indexed point ids, shape (219,)
        "cap_pts": pts[idx0].astype(np.int64),  # (219, 7) coords over F_3
        "size": len(cols1),
        "all_pts": pts,
    }


def verify_build() -> dict:
    """Sanity: rebuild G0 via qmsd and confirm shape/rank, and confirm the dual is
    RM_3(9,7) by checking the weight-18 generator membership.  Returns a small report."""
    import qmsd.reedmuller as rm
    import qmsd.triorthogonal as tri

    cap = load_cap()
    G = rm.rm_generator(R, M, P)
    built = tri.build_triorthogonal_code(P, M, R, cap["cols1"], G=G)
    G0 = np.asarray(built["X_stab"]).astype(np.int64) % P
    return {
        "G0_shape": tuple(G0.shape),
        "full_rank": bool(built["full_rank"]),
        "n": int(built["n"]),
        "k": int(built["k"]),
        "rtilde": RTILDE,
        "dual_code": f"RM_{P}({RTILDE},{M})",
        "d_rm_min": D_RM_MIN,
    }
