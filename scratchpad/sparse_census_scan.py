"""Scan sparse normal bases / encoders of the [[7,1,3]] s = 1 gadget: CSUMs, depth, Z-type exposures (tau_F = 0, 15),
inner failure coefficient C_3 and the census floor at p = 1e-3, D = 21.  Output: runs/meas_distill/sparse_census_d{d}.jsonl"""
import json, sys
import numpy as np
from joblib import Parallel, delayed
from primepower_msd import host_census as h
from primepower_msd.meas_sim import leading_order

def one(d, s):
    try:
        P, C, i = h.gadget(d) if s < 0 else h.gadget(d, True, seed=s, tries=1)
    except Exception as e:
        return dict(seed=s, error=repr(e))
    r = dict(seed=s, d=d, N=[int(x) for x in np.asarray(P.N).ravel()], csums=sum(op[0] == "CSUM" for op in C.ops))
    for tf in (0, 15):
        F, T = h.faults_logical(C, i, d, tf); N, NZ, W = h.census(C, i, d, F)
        r[f"T{tf}"] = T; r[f"gate{tf}"] = NZ["gate"]; r[f"rounds{tf}"] = NZ["idle"] + NZ["fourier"] + NZ["meas"] + NZ["prep"]
        nz = dict(gate=NZ["gate"], idle=NZ["idle"], fourier=NZ["fourier"], meas=NZ["meas"], prep=NZ["prep"])
        r[f"floor{tf}_D21"] = h.floor(1e-3, 21, nz, tf)
    r["C3"] = leading_order(P, 3)["C"]
    return r

d, n = int(sys.argv[1]), int(sys.argv[2])
out = f"runs/meas_distill/sparse_census_d{d}.jsonl"
for r in Parallel(n_jobs=8, return_as="generator_unordered")(delayed(one)(d, s) for s in range(-1, n)):
    with open(out, "a") as f:
        f.write(json.dumps(r) + "\n")
