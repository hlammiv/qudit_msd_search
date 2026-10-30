"""Exact performance of every Reed-Solomon magic code [[p-k, k, r-k+2]]_p (r = r_max, 1 <= k <= r) for
all primes up to 97, from the closed-form MDS enumerators (qmsd.rs_family): leading order, exact
threshold (per-output and block), acceptance at 1e-2 / 1e-3, gamma.  Output: runs/coverage/rs_family_all.json"""
import json, math, time, sys
from qmsd.rs_family import rs_profile, rs_parameters, r_max_rs
from qmsd.performance import leading_order, threshold
primes = [q for q in range(5, 100) if all(q % t for t in range(2, int(q ** .5) + 1))]
rows = []; t0 = time.time()
for p in primes:
    r = r_max_rs(p)
    for k in range(1, r + 1):
        n, _, d = rs_parameters(p, k)
        prof = rs_profile(p, k)
        la, lb = leading_order(prof, "avg"), leading_order(prof, "block")
        ta, tb = threshold(prof, "avg"), threshold(prof, "block")
        rows.append({"p": p, "r": r, "n": n, "k": k, "d": d, "dim_A": r + 1 - k, "gamma": math.log(n / k) / math.log(d) if d > 1 else None,
                     "c_eps_avg": la.c_eps, "L_d_avg": str(la.L_d), "c_eps_block": lb.c_eps, "threshold_avg": ta.eps, "threshold_block": tb.eps,
                     "P_acc_1e-2": float(prof.p_accept(1e-2)), "P_acc_1e-3": float(prof.p_accept(1e-3)),
                     "eps_out_avg_1e-2": float(prof.eps_out(1e-2, "avg")), "eps_out_avg_1e-3": float(prof.eps_out(1e-3, "avg"))})
        print(f"p={p:2d} [[{n},{k},{d}]] gamma={rows[-1]['gamma'] if rows[-1]['gamma'] is None else round(rows[-1]['gamma'],3)} c={la.c_eps:.3g} eps*={ta.eps if ta.eps is None else round(ta.eps,4)} Pacc(1e-2)={rows[-1]['P_acc_1e-2']:.3f} [{time.time()-t0:.0f}s]", flush=True)
json.dump({"description": "Reed-Solomon (m = 1) magic codes, exact closed-form MDS performance (qmsd.rs_family); eps = paper convention", "rows": rows}, open("runs/coverage/rs_family_all.json", "w"), indent=1)
print("rows", len(rows))
