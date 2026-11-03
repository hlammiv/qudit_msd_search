"""Build and verify the planar witness found by levelset_lines.py: the dual word u = (product over the
excluded parameter values) on (p-b) parallel lines of one plane, including the full lines of O.
Checks (no big matrices): sum_x u(x) x^e = 0 mod p for every monomial of degree <= r (u in RM(r,4)^perp,
so u restricted off O is a kernel word of G0), and weight off O = lower bound.
usage: levelset_witness.py CODE.json LINES.json OUT.json"""
import sys, json, itertools, numpy as np
from qmsd.reedmuller import points
src = json.load(open(sys.argv[1])); b = src.get("best", src); L = json.load(open(sys.argv[2])); f = L["found"]
p, m, r = b["p"], 4, b["r"]; O = set(c - 1 for c in b["puncture_columns_1indexed"])
P = np.asarray(points(m, p), dtype=np.int64); index = {tuple(int(t) for t in x): i for i, x in enumerate(P)}
v = np.array(f["direction"]); j = int(np.nonzero(v)[0][0])
def lift(q): q = list(q); return np.array(q[:j] + [0] + q[j:])
base, d3 = np.array(f["quotient_base"]), np.array(f["quotient_line_dir"])
full_t = sorted({next(t for t in range(p) if tuple((base + t * d3) % p) == tuple(q)) for q in f["quotient_points"]})
a1, b1 = divmod(m * (p - 1) - r - 1, p - 1); nl_word = p - b1
plane = {t: [index[tuple(int(c) for c in (lift((base + t * d3) % p) + s * v) % p)] for s in range(p)] for t in range(p)}
occ = {t: sum(1 for x in plane[t] if x in O) for t in range(p)}
chosen = full_t + sorted([t for t in range(p) if t not in full_t], key=lambda t: -occ[t])[: nl_word - len(full_t)]
excluded = [t for t in range(p) if t not in chosen]
u = {}
for t in chosen:
    val = 1
    for e in excluded: val = (val * (t - e)) % p
    for x in plane[t]: u[x] = val
idx = np.array(list(u)); vals = np.array([u[x] for x in idx]); Px = P[idx]
bad = 0
POW = np.array([[pow(x, k, p) for k in range(p)] for x in range(p)], dtype=np.int64)   # exact modular powers (no int64 overflow)
for e in itertools.product(range(p), repeat=m):
    if sum(e) > r: continue
    term = vals.copy()
    for i in range(m):
        term = (term * POW[Px[:, i], e[i]]) % p
    if int(term.sum() % p): bad += 1
w_off = sum(1 for x in u if x not in O)
res = {"p": p, "word_weight": len(u), "points_in_O": len(u) - w_off, "kernel_word_weight": w_off, "lower_bound": L["lower_bound"],
       "orthogonal_to_RM_r": bad == 0, "exact": bad == 0 and w_off == L["lower_bound"],
       "witness_points_off_O_0based_grid": sorted(int(x) for x in u if x not in O), "coefficients": [int(u[x]) for x in sorted(u) if x not in O]}
json.dump(res, open(sys.argv[3], "w")); print({k: v for k, v in res.items() if not k.startswith("witness") and k != "coefficients"}, flush=True)
