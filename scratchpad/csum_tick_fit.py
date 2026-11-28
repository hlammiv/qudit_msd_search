"""Joint fit of the tick-idle (f=1) logical CSUM (Pauli-frame 'seq' decoder) and two-patch memory ('base'):
log p_fail = log A_k + ((D+1)/2)(log p - log p_th), shared p_th, separate A for seq / base, over sub-threshold points
(>= 5 events, p_fail < 0.1) of runs/surface/csum_d4_tick1.jsonl.  python scratchpad/csum_tick_fit.py [Dmin]
Basis 'avg' = Z+X fails over Z+X shots (per-basis mean, the CSUM_Z4 / BASE2_Z4 convention); 'X' = the host-census input."""
import json, math, sys
from collections import defaultdict
import numpy as np

def load(path="runs/surface/csum_d4_tick1.jsonl"):
    agg = defaultdict(lambda: [0, 0])
    for l in open(path):
        r = json.loads(l)
        k = "seq" if r["decoder"] == "seq" else ("csum" if r["csum"] else "base")
        for b in (r["basis"], "avg"):
            a = agg[(k, b, r["D"], r["p"])]; a[0] += r["fails"]; a[1] += r["shots"]
    return agg

def joint(agg, basis, Dmin=5, keys=("seq", "base")):
    X, y, w = [], [], []
    for (k, b, D, p), (f, s) in agg.items():
        if b != basis or k not in keys or D < Dmin or f < 5 or f / s >= 0.1:
            continue
        n = (D + 1) / 2
        X.append([1.0 * (k == kk) for kk in keys] + [n]); y.append(math.log(f / s) - n * math.log(p)); w.append(math.sqrt(f))
    X, y, w = map(np.array, (X, y, w))
    sol, *_ = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)
    r = y - X @ sol
    return {**{kk: round(math.exp(sol[i]), 4) for i, kk in enumerate(keys)}, "p_th": float(f"{math.exp(-sol[-1]):.3g}"),
            "n": len(y), "rms": round(float(np.sqrt(np.mean(r ** 2))), 3)}

if __name__ == "__main__":
    Dmin = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    agg = load()
    for b in ("avg", "Z", "X"):
        print(b, f"D>={Dmin}", joint(agg, b, Dmin))
