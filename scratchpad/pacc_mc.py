"""Acceptance probability P_acc(eps) by Monte Carlo for the large new codes (exact enumeration is
out of reach: p^dim(A) with dim A = 13-61), plus the leading approximation (1-eps)^n. Only the
acceptance is sampled (cheap); the output error at these n needs the leading-order census."""
import json, math, numpy as np
from qmsd.results import load_result_catalog
from qmsd.performance import record_matrices
from qmsd.weighted import weighted_support, quadric, weighted_code
from qmsd.mc_distill import simulate
codes = {}
cat = list(load_result_catalog())
def rec(label): return next(r for r in cat if r.label == label)
for lab in ("[[237,52,6]]_17", "[[293,68,5]]_19", "[[296,65,6]]_19", "[[424,105,6]]_23"):
    G0, H1 = record_matrices(rec(lab)); codes[lab] = (int(lab.split("_")[-1]), G0, H1)
GX, X, mu = weighted_support(7, 4, quadric(7, 4, (1, 1, 1, 1)), 7)
for lab, f in (("[[1756,260,7]]_7", "runs/weighted/p7_m4_sphere_1756_260_cols.json"), ("[[1746,270,>=7]]_7", "runs/certify/p7_k270_seed407.json")):
    cols = json.load(open(f))["cols_in_X"]; G0, H1, _ = weighted_code(GX, mu, cols, 7); codes[lab] = (7, G0, H1)
out = {}
for lab, (p, G0, H1) in codes.items():
    n = G0.shape[1]; out[lab] = {"n": n, "k": int(H1.shape[0]), "dim_A": int(G0.shape[0]), "P_acc": {}}
    for eps in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        shots = 200_000 if n > 1000 else 400_000
        mc = simulate(G0, H1, p, eps, shots=shots, seed=1, batch=20_000)
        out[lab]["P_acc"][str(eps)] = {"mc": mc["P_acc"], "err": mc["P_acc_err"], "no_error": (1 - eps) ** n, "eps_out_block_mc": mc["eps_out_block"], "accepted": mc["accepted"]}
        print(f"{lab:20s} n={n:4d} dimA={G0.shape[0]:3d} eps={eps:.0e}: P_acc={mc['P_acc']:.4f}±{mc['P_acc_err']:.4f}  (1-eps)^n={(1-eps)**n:.4f}  block eps_out(MC)={mc['eps_out_block']:.2e} acc={mc['accepted']}", flush=True)
json.dump(out, open("runs/certify/pacc_mc.json", "w"), indent=1)
