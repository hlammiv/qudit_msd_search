"""Lenore sweep: the ring punctured-RM(2,9) T-like family [[512-k, k, d]]_4.

For each puncture set P (plane-spread/Sidon sets: no 4 points in a 2-flat => no 5 in a 3-flat => d_Z <= 4
is attainable; random information sets up to k = dim RM(2,9) = 46), build the Z_4 lift (M_X = shortened
RM(2,9), logical X = quotient reps, M_Z = (M_X + L)^perp), check hosting of the carrier (automatic by the
quadruply-even theorem; re-verified), the logical-escape parity, certify d <= 4 with the Z_4 MITM, record
gamma and the max 3-flat occupancy of P (the geometric distance predictor: d_Z = 8 - max occupancy when
the min-weight dual words dominate).

    PYTHONPATH=. python -m primepower_msd.rm_ring_sweep --jobs 28 --out runs/diverse/rm29_sweep.json
"""

from __future__ import annotations

import argparse
import itertools as it
import json
import math
import random
import time

import numpy as np

from . import lifts, rgt
from . import phase_tables as PT
from .z4_distance import code_distance


def max_3flat_occupancy(P: list[int], m: int = 9) -> int:
    """max over 3-flats F of |F ∩ P| (F spanned by 4 affinely independent points of P)."""
    pts = np.array([[(p >> i) & 1 for i in range(m)] for p in P], dtype=np.int64)
    Pset = set(P)
    best = min(4, len(P))
    for idx in it.combinations(range(len(P)), 4):
        a, b, c, d = pts[list(idx)]
        # affine span of {a,b,c,d}: a + span{b-a, c-a, d-a} (over F_2, -a = a)
        dirs = [(b ^ a), (c ^ a), (d ^ a)]
        flat = set()
        for coeffs in it.product(range(2), repeat=3):
            v = a.copy()
            for cf, dv in zip(coeffs, dirs):
                if cf:
                    v = v ^ dv
            flat.add(int(sum(int(v[i]) << i for i in range(m))))
        if len(flat) < 8:
            continue                      # dependent 4-set: lies in a 2-flat (occupancy counted elsewhere)
        occ = len(flat & Pset)
        if occ > best:
            best = occ
    return best


def one(P: list[int], label: str, d_max: int = 4) -> dict:
    t0 = time.time()
    q = np.array(PT.CARRIER_Q)
    cand = lifts.rm_puncture_candidate(9, P, 2)
    if cand is None:
        return {"label": label, "P": P, "status": "rank fail"}
    c = lifts.lift_css_from_x(cand["MX"], cand["Xlog"], cand["n"], label=label)
    h = rgt.hosts(c, rgt.TransversalLayer.uniform(c.n, q), "support")
    par = rgt.axis_parity_report(c)
    try:
        d, cert, det = code_distance(c, d_max=d_max, sides=("Z",))   # X-logicals have weight >= 2^7 - k
    except MemoryError:
        d, cert, det = None, False, {}
    g = math.log(c.n / c.k) / math.log(d) if (d and d > 1 and c.k > 0) else None
    return {"label": label, "P": [int(x) for x in P], "n": c.n, "k": c.k, "k_genuine": c.k_genuine, "soft": c.soft,
            "free": c.free_quotient_ok, "hosted": bool(h.hosted), "escape_possible": par["logical_escape_possible"],
            "distance": d, "certified": cert, "detail": det, "gamma": g,
            "max_3flat_occupancy": max_3flat_occupancy(P) if len(P) <= 46 else None,
            "seconds": round(time.time() - t0, 1), "status": "ok"}


def puncture_sets(seed: int, n_each: int = 4):
    from qmsd.sampling import random_plane_spread, all_points, points_to_columns
    rng = random.Random(seed)
    allpts = all_points(9, 2)
    jobs = []
    for k in (12, 16, 20, 24, 28):
        for s in range(n_each):
            pts = random_plane_spread(9, 2, k, rng, allpts)
            if pts is None:
                continue
            P = sorted(int(cc) - 1 for cc in points_to_columns(pts, 2))
            jobs.append((P, f"plane-spread k={k} s={s}"))
    for k in (30, 36, 40, 44, 46):
        for s in range(n_each):
            for _ in range(20):
                P = sorted(rng.sample(range(512), k))
                if lifts.rm_puncture_candidate(9, P, 2) is not None:
                    jobs.append((P, f"random k={k} s={s}")); break
    return jobs


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--each", type=int, default=4)
    ap.add_argument("--out", default="runs/diverse/rm29_sweep.json")
    a = ap.parse_args(argv)
    from joblib import Parallel, delayed
    jobs = puncture_sets(a.seed, a.each)
    print(f"{len(jobs)} puncture sets", flush=True)
    t0 = time.time()
    res = Parallel(n_jobs=a.jobs, backend="loky")(delayed(one)(P, lab) for (P, lab) in jobs)
    import os
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as f:
        json.dump({"seconds": round(time.time() - t0, 1), "results": res}, f, indent=1, default=str)
    ok = [r for r in res if r.get("status") == "ok" and r.get("distance")]
    for r in sorted(ok, key=lambda r: r["gamma"] or 99)[:15]:
        print(f"{r['label']:26s} [[{r['n']},{r['k']},{r['distance']}]] gamma={r['gamma']:.3f} occ3={r['max_3flat_occupancy']} hosted={r['hosted']} esc={r['escape_possible']}")


if __name__ == "__main__":
    main()
