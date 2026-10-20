"""Top-down exhaustive search for failing flats (see rm_puncture_cert.FlatBFS) — certificate driver.

    python -m primepower_msd.rm_cert_bfs --points runs/diverse/sidon82_F2_13.json --codim 5 --jobs 30 \
        --out runs/diverse/cert_sidon82.json

If no flat of codimension `codim` fails (X_G = 0 for all of them), then every Z-logical S = A ∪ R with
dim aff(R) <= m - codim has A ⊆ (a flat of dim m - codim containing R) and the case-2 (local flat)
argument finishes the certificate.  For the RM(3,13) constructions the local lemma covers flats of
dim <= 8, so --codim 5 suffices (Sidon-82: no failing flat already at codim 4 => d_Z = 10 certified).
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time

import numpy as np

from .rm_puncture_cert import FlatBFS, canon_forms

_B: FlatBFS | None = None


def _init(points):
    global _B
    _B = FlatBFS(np.array(points, dtype=np.uint8), 3)


def _has_const(key, m):
    v = 1 << m
    basis = {x.bit_length() - 1: x for x in key}
    while v:
        h = v.bit_length() - 1
        if h not in basis:
            return False
        v ^= basis[h]
    return True


def _children(parent):
    forms = list(parent)
    _, basis, _ = _B.check(forms)
    m = _B.m
    seen = set()
    out = []
    c = len(forms)
    for g in _B.val:
        key = canon_forms(forms + [g])
        if len(key) != c + 1 or key in seen:
            continue
        seen.add(key)
        if _has_const(key, m):
            continue
        fails, _, _ = _B.check(forms + [g], basis)
        if fails:
            out.append(key)
    return out, len(seen)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", required=True)
    ap.add_argument("--codim", type=int, required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    pts = json.load(open(a.points))["points"]
    _init(pts)
    t0 = time.time()
    level = []
    for f in _B.val:
        fails, _, _ = _B.check([f])
        if fails:
            level.append((f,))
    log = [{"codim": 1, "failing": len(level), "seconds": round(time.time() - t0, 1)}]
    print(log[-1], flush=True)
    with mp.Pool(a.jobs, initializer=_init, initargs=(pts,)) as pool:
        for c in range(2, a.codim + 1):
            t = time.time()
            nxt = set()
            tested = 0
            for kids, nt in pool.imap_unordered(_children, level, chunksize=4):
                nxt.update(kids)
                tested += nt
            level = sorted(nxt)
            log.append({"codim": c, "children_tested_with_dups": tested, "failing": len(level),
                        "seconds": round(time.time() - t, 1)})
            print(log[-1], flush=True)
            if not level:
                break
    res = {"points_file": a.points, "target_codim": a.codim, "levels": log,
           "failing_at_target": [list(k) for k in level] if len(log) == a.codim and log[-1]["codim"] == a.codim else [],
           "certified_no_failing_flat": (log[-1]["failing"] == 0)}
    json.dump(res, open(a.out, "w"), indent=1)
    print("certified_no_failing_flat:", res["certified_no_failing_flat"])


if __name__ == "__main__":
    main()
