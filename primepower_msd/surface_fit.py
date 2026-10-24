"""Threshold and logical-error scaling from runs/surface/zd_d{d}.jsonl.

Table of p_fail(D, p) per basis, threshold estimate (crossing of neighbouring D pairs, linear interpolation in log p),
and a fit of the footprint model p_L(D) = A (p / p_th)^((D+1)/2) to the sub-threshold points (fails >= 5, p_fail < 0.1).

python -m primepower_msd.surface_fit d [basis]
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict

import numpy as np


def load(d):
    rows = [json.loads(l) for l in open(f"runs/surface/zd_d{d}.jsonl")]
    agg = defaultdict(lambda: [0, 0])
    for r in rows:
        a = agg[(r["basis"], r["D"], r["p"])]; a[0] += r["fails"]; a[1] += r["shots"]
    return {k: (f, s) for k, (f, s) in agg.items()}


def crossings(pts, basis):
    Ds = sorted({D for b, D, p in pts if b == basis})
    out = []
    for D1, D2 in zip(Ds[:-1], Ds[1:]):
        ps = sorted({p for b, D, p in pts if b == basis and D == D1} & {p for b, D, p in pts if b == basis and D == D2})
        diff = []
        for p in ps:
            f1, s1 = pts[(basis, D1, p)]; f2, s2 = pts[(basis, D2, p)]
            if f1 and f2:
                diff.append((p, math.log(f2 / s2) - math.log(f1 / s1)))
        for (pa, da), (pb, db) in zip(diff[:-1], diff[1:]):
            if da < 0 <= db:
                la, lb = math.log(pa), math.log(pb)
                out.append((D1, D2, math.exp(la + (lb - la) * (-da) / (db - da))))
    return out


def fit(pts, basis):
    X, y, w = [], [], []
    for (b, D, p), (f, s) in pts.items():
        if b != basis or f < 5 or f / s >= 0.1:
            continue
        X.append([1.0, (D + 1) / 2, (D + 1) / 2 * math.log(p)]); y.append(math.log(f / s)); w.append(math.sqrt(f))
    if len(X) < 4:
        return None
    X, y, w = np.array(X), np.array(y), np.array(w)
    # log P = logA + ((D+1)/2) (log p - log p_th): two parameters, column 2 has coefficient 1 by construction
    A_ = np.stack([X[:, 0], X[:, 1]], axis=1)
    rhs = y - X[:, 2]
    sol, *_ = np.linalg.lstsq(A_ * w[:, None], rhs * w, rcond=None)
    logA, m = sol                       # m = -log p_th
    resid = rhs - A_ @ sol
    return dict(A=math.exp(logA), p_th=math.exp(-m), n=len(y), rms=float(np.sqrt(np.mean(resid ** 2))))


def main():
    d = int(sys.argv[1]); bases = [sys.argv[2]] if len(sys.argv) > 2 else ["Z", "X"]
    pts = load(d)
    for basis in bases:
        Ds = sorted({D for b, D, p in pts if b == basis}); ps = sorted({p for b, D, p in pts if b == basis})
        print(f"== d={d} basis {basis}: p_fail after D rounds (fails/shots)")
        print("p        " + "".join(f"D={D:<14}" for D in Ds))
        for p in ps:
            row = f"{p:<9.4g}"
            for D in Ds:
                f, s = pts.get((basis, D, p), (0, 0))
                row += f"{f / s if s else float('nan'):<8.2e}{f:>5}  " if s else " " * 15
            print(row)
        print("crossings:", [(a, b, f"{x:.2e}") for a, b, x in crossings(pts, basis)])
        print("fit p_L = A (p/p_th)^((D+1)/2):", fit(pts, basis))


if __name__ == "__main__":
    main()
