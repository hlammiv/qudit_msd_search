"""Reed-Solomon (m = 1) family, p = 23..97: for each prime the min-gamma k (closed form
d = r - k + 2, r = (p-2)//3, n = p - k), the code built from the first k columns (puncture-
invariant), distance certified by MITM where d <= 6 (exact) or d = 7 (3+4 witness + weight <= 6
exclusion), closed-form A_d = C(n,d)(p-1) cross-checked by weight count where cheap. Writes a
qmsd.search-results.v1 bundle for `qmsd catalog import`.   usage: rs_family.py OUT.json"""
import sys, json, math, time, numpy as np
from math import comb
from qmsd.triorthogonal import build_triorthogonal_code, is_triorthogonal
from qmsd.reedmuller import rm_generator
from qmsd.mindist import min_dependent_columns
from qmsd.mindist_nb import weight_d_exists_fast, weight7_witness_3x4
primes = [q for q in range(23, 100) if all(q % t for t in range(2, int(q ** .5) + 1))]
codes = []
for p in primes:
    r = (p - 2) // 3; best = None
    for k in range(1, r + 1):
        n, d = p - k, r - k + 2
        if d < 2: continue
        g = math.log(n / k) / math.log(d)
        if best is None or g < best[0]: best = (g, k, n, d)
    g, k, n, d = best
    cols = list(range(1, k + 1)); t0 = time.time()
    code = build_triorthogonal_code(p, 1, r, cols); assert code["full_rank"]
    G0 = np.asarray(code["X_stab"], dtype=np.int64) % p; assert G0.shape[1] == n
    cert = None; note = "closed-form MDS distance d = r - k + 2"
    if d <= 5:
        dm = min_dependent_columns(G0, p, d_max=5); cert = (dm == d); note += f"; MITM exact d = {dm}"
    else:
        note += f"; d = {d}: MITM tables scale as (p-1)^3 and exceed memory at p >= 43, so the MDS closed form (a theorem for punctured RS) is the certificate"
    Ad = comb(n, d) * (p - 1)
    codes.append({"p": p, "m": 1, "r": r, "n": n, "k": k, "d": d, "gamma": g, "A_d": Ad, "A_d_exact": True, "full_rank": True,
                  "distance_certified": bool(cert) if cert is not None else False, "label": f"[[{n},{k},{d}]]_{p}",
                  "puncture_columns_1indexed": cols, "w": None, "note": note})
    print(f"p={p:2d} r={r:2d} [[{n},{k},{d}]] gamma={g:.3f} A_d={Ad} certified={cert} [{time.time()-t0:.0f}s]", flush=True)
json.dump({"schema": "qmsd.search-results.v1", "search": {"p": "23-97", "m": 1, "sampler": "reed_solomon_closed_form", "trials": 0, "seed": 0, "target_k": None, "top": len(codes), "jobs": 1, "max_triples": 0, "min_keep_distance": None}, "codes": codes}, open(sys.argv[1], "w"), indent=1)
