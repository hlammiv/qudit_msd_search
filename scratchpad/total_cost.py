"""Combined cost: R * (raw magic states per output) + (two-qudit Clifford gates per output), minimised over
protocols (chains of rounds), for a range of R = cost of one magic state in Clifford gates.
Per round attempt a code costs n raw states (inputs) and g gates (one-shot axis-by-axis count for grid
RM/RS codes, systematic encoder count for the cyclic / doubled codes).  Bookkeeping per output of round i:
  M_i = a_i M_{i-1},  Y_i = a_i Y_{i-1} + g_i / (k_i P_i),  a_i = n_i / (k_i P_i).
Curves: exact (standard library, RS family), rigorous census-free (CorrectingBoundCurve t = 0) for the
mid and large search / level-set codes, leading-order surrogates from censuses for the p = 7 sphere codes
(flagged).  usage: total_cost.py OUT.json"""
import sys, json, math, itertools, numpy as np, mpmath
from qmsd.protocol import CodeCurve, curve, _rm_k1, CYCLIC
from qmsd.cyclic import construction_A, construction_B
from qmsd.doubling import paper_codes
from qmsd.rs_family import rs_profile, rs_parameters, r_max_rs
from qmsd.objective import LeadingOrderCurve
from qmsd.rigorous_curve import CorrectingBoundCurve
from scratchpad.fast_encoder import stage_gates, add_gates, line_matrix
from scratchpad.fast_matrices import _rref

def systematic_gates(G0, H1, p):
    G0 = np.ascontiguousarray(np.asarray(G0, dtype=np.int64) % p); H1 = np.asarray(H1, dtype=np.int64) % p
    Gs, piv = _rref(G0, p)
    H = H1.copy()
    for i, c in enumerate(piv): H = (H - np.outer(H[:, c], Gs[i])) % p
    return int((Gs != 0).sum() + (H != 0).sum())

def library_with_gates(p):
    pool = []
    G0, H1 = _rm_k1(p, 1); c = curve(f"[[{p-1},1,?]]_{p} RS1", G0, H1, p)
    if c: pool.append((c, systematic_gates(G0, H1, p), "exact"))
    for label, cons, n, S in CYCLIC.get(p, []):
        if n > 60: continue
        G0, H1 = (construction_A if cons == "A" else construction_B)(p, n, S); c = curve(label, G0, H1, p)
        if c: pool.append((c, systematic_gates(G0, H1, p), "exact"))
    if p >= 7 and p not in CYCLIC:
        nu, rr = (p + 1) // 2, p // 6
        G0, H1 = construction_B(p, p + 1, set(range(nu - rr, nu + rr + 1))); c = curve(f"[[{p},1,{2*rr+1}]]_{p} cyc-B", G0, H1, p)
        if c: pool.append((c, systematic_gates(G0, H1, p), "exact"))
    for label, G0, H1, coeffs in paper_codes(p):
        c = curve(label, G0, H1, p, coeffs=coeffs)
        if c: pool.append((c, systematic_gates(G0, H1, p), "exact"))
    for kk in range(1, r_max_rs(p) + 1):     # RS family: one-shot on the line F_p, RM(r,1)
        nn, _, dd = rs_parameters(p, kk)
        if dd >= 2:
            g = stage_gates(p, 1, r_max_rs(p))[0]
            pool.append((CodeCurve(f"[[{nn},{kk},{dd}]]_{p} RS", p, nn, kk, rs_profile(p, kk)), g, "exact"))
    return pool

def best(pool, eps0, target, R, depth=5):
    bestv = [math.inf, None]
    def dfs(eps, M, Y, seq):
        if eps <= target:
            v = R * M + Y
            if v < bestv[0]: bestv[0] = v; bestv[1] = (list(seq), M, Y)
            return
        if len(seq) >= depth or R * M + Y >= bestv[0]: return
        for c, g, kind in pool:
            pa = float(c.pacc_lo(eps) if hasattr(c, "pacc_lo") and not getattr(c, "exact", True) else c.pacc(eps))
            if pa <= 1e-12: continue
            e2 = float(c.f_hi(eps) if not getattr(c, "exact", True) else c.f(eps))
            if not (e2 < eps * 0.9): continue
            a = c.n / (c.k * pa)
            dfs(e2, a * M, a * Y + g / (c.k * pa), seq + [c.name])
    dfs(eps0, 1.0, 0.0, [])
    return bestv

out = {}
for p in (7, 11, 17):
    pool = library_with_gates(p)
    big = []
    if p == 7:
        g7 = stage_gates(7, 4, 7)[0]
        big = [(LeadingOrderCurve("[[1746,270,7]]_7 (sphere, surrogate)", 7, 1746, 270, 7, 193.2), g7, "surrogate"),
               (LeadingOrderCurve("[[1756,260,7]]_7 (sphere, surrogate)", 7, 1756, 260, 7, 1.292), g7, "surrogate")]
    if p == 11:
        big = [(CorrectingBoundCurve("[[13332,1309,22]]_11 (level set, rigorous)", 11, 13332, 1309, 22, 0), stage_gates(11, 4, 13)[0], "rigorous")]
    if p == 17:
        big = [(CorrectingBoundCurve("[[237,52,6]]_17 (m=2, rigorous)", 17, 237, 52, 6, 0), stage_gates(17, 2, 10)[0], "rigorous"),
               (CorrectingBoundCurve("[[78455,5066,68]]_17 (level set, rigorous)", 17, 78455, 5066, 68, 0), stage_gates(17, 4, 21)[0], "rigorous")]
    out[p] = {"pool": [(c.name, c.n, c.k, g, kind) for c, g, kind in pool + big], "rows": []}
    print(f"\n=== p = {p}: {len(pool)} small codes + {len(big)} large; gates per block: " + ", ".join(f"{c.name.split(' ')[0]} {g}" for c, g, _ in pool[:4] + big), flush=True)
    for eps0, tgt in [(1e-3, 1e-15), (1e-3, 1e-30)]:
        for R in (1, 10, 100, 1000, 10000):
            v0 = best(pool, eps0, tgt, R); v1 = best(pool + big, eps0, tgt, R)
            row = {"eps0": eps0, "target": tgt, "R": R, "small_only": v0[0], "with_large": v1[0], "saving": v0[0] / v1[0] if v1[0] else None,
                   "protocol": v1[1][0] if v1[1] else None, "magic_per_output": v1[1][1] if v1[1] else None, "gates_per_output": v1[1][2] if v1[1] else None}
            out[p]["rows"].append(row)
            print(f"  {eps0:g}->{tgt:g} R={R:>5}: small codes {v0[0]:10.1f} | with large {v1[0]:10.1f} ({row['saving']:.2f}x) {[s.split(' ')[0] for s in row['protocol']]} M={row['magic_per_output']:.1f} G={row['gates_per_output']:.0f}", flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1, default=str)
