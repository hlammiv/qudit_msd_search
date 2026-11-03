"""Joint (both-sector) gap cut of an escape variant on trivial frames: Z- and X-basis runs with the same seed share their
noise (Escape.run(chunk_seeded=True)), so keep = one escape's postselection and each run's gap is its sector's gap.
Reports, per (cut_Z, cut_X), the joint keep rate and failures (either readout wrong).
    PYTHONPATH=. python scratchpad/escape_joint_cut.py variant shots p [seed] -> runs/surface/escape_joint_cut.jsonl"""
import json
import sys
import time

import numpy as np
from joblib import Parallel, delayed

from primepower_msd.graft_run import VARIANTS
from primepower_msd.graft_zd import Escape


def one(variant, basis, shots, p, seed):
    ex = Escape(4, basis, **VARIANTS[variant])
    r = ex.run(p, shots, seed=seed, decoder="peel-gap", inject=np.zeros((shots, 19, 2), dtype=np.int64), per_shot=True,
               chunk_seeded=True)
    return r["shot_keep"], r["shot_flip"], np.nan_to_num(r["shot_gap"], nan=0.0)


def main():
    variant, shots, p = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    t = time.time()
    (kz, fz, gz), (kx, fx, gx) = Parallel(n_jobs=2)(delayed(one)(variant, b, shots, p, seed) for b in "ZX")
    keep = kz & kx
    bad = (fz != 0) | (fx != 0)
    cuts = (0, 10, 11, 12, 13, 14, 15, 16, 18)
    grid = {f"{cz},{cx}": [int((keep & (gz >= cz) & (gx >= cx)).sum()), int((keep & (gz >= cz) & (gx >= cx) & bad).sum())]
            for cz in cuts for cx in cuts}
    out = dict(variant=variant, p=p, shots=shots, seed=seed, keep=int(keep.sum()), keep_disagree=float((kz != kx).mean()),
               grid=grid, gap_max=[float(gz.max()), float(gx.max())], seconds=round(time.time() - t, 1))
    with open("runs/surface/escape_joint_cut.jsonl", "a") as f:
        f.write(json.dumps(out) + "\n")
    best = sorted(((v[0], k, v[1]) for k, v in grid.items()), reverse=True)
    print(variant, "keep", keep.mean(), "gap max", out["gap_max"])
    for kept, k, fl in best[:40]:
        print("  cut(Z,X)", k, "keep %.3f" % (kept / shots), "fails", fl)


if __name__ == "__main__":
    main()
