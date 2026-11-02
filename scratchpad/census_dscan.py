"""Census exposures (Z-type, per output) of the dense and best sparse [[7,1,3]] gadget at d = 4..32, plus the two-wire
CSUM failure census (host_census.faults_pairs, 'uniform' and 'copy') -> runs/meas_distill/census_dscan.jsonl.

python scratchpad/census_dscan.py 4,8,16,32
"""
import json
import sys

from primepower_msd import host_census as H


def row(d, sparse):
    P, C, inputs = H.gadget(d, sparse, seed=0, tries=1) if sparse else H.gadget(d)
    F, T_end, = H.faults_logical(C, inputs, d)[:2]
    N, NZ, W = H.census(C, inputs, d, F)
    r = dict(d=d, gadget="sparse-best" if sparse else "dense", n_csum=sum(op[0] == "CSUM" for op in C.ops), depth=T_end,
             gate=NZ["gate"] / P.n_in, rounds=(NZ["idle"] + NZ["fourier"] + NZ["meas"] + NZ["prep"]) / P.n_in)
    for mode in ("uniform", "copy"):
        N2, NZ2, W2 = H.census(C, inputs, d, H.faults_pairs(C, d, mode))
        r[f"two_wire_{mode}"] = NZ2["gate2"] / P.n_in
    return r


if __name__ == "__main__":
    for d in [int(x) for x in sys.argv[1].split(",")]:
        for sparse in (False, True):
            r = row(d, sparse)
            with open("runs/meas_distill/census_dscan.jsonl", "a") as f:
                f.write(json.dumps(r) + "\n")
            print(r, flush=True)
