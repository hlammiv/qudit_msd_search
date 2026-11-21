"""Per-2-adic-level complementary gaps of the f = 1 closed-window escape on CLEAN frames (MEAS_DISTILL 5k, step 3):
the escape's own failures at high gap cut are order-2 (top level).  Records per shot, for both bases (same seed per slice):
keep, decoded flip (mod d), syndrome-window gap (min over levels) and the per-level gaps.  Checkpointed per slice.
    PYTHONPATH=.:scratchpad python scratchpad/escape_level_gaps.py p shots jobs out_dir [priors.npz]
"""
import os
import sys
import time

import numpy as np
from joblib import Parallel, delayed

from cult_escape_e2e import D
from primepower_msd.graft_run import VARIANTS
from primepower_msd.graft_zd import Escape

VARIANT = os.environ.get("ELG_VARIANT", "one+gap-tick1-closed")
SLICE = int(os.environ.get("ELG_SLICE", 10_000))
SEED_OFF = int(os.environ.get("ELG_SEED", 500000))
_cache = {}


def run_slice(basis, n, p, seed, priors):
    if basis not in _cache:
        pri = None
        if priors:
            z = np.load(priors); pri = (z["frames"].astype(np.int64), z["probs"])
        _cache[basis] = Escape(D, basis, **VARIANTS[VARIANT], priors=pri)
    r = _cache[basis].run(p, n, seed=seed, decoder="peel-gap", inject=np.zeros((n, 19, 2), dtype=np.int64), per_shot=True,
                          chunk_seeded=True, flag_cert=True)
    return dict(keep=r["shot_keep"], flip=r["shot_flip"].astype(np.int8), gap=np.nan_to_num(r["shot_gap_syn"], nan=0).astype(np.float32),
                lev=np.nan_to_num(r["shot_gap_syn_levels"], nan=-1, posinf=99).astype(np.float32))


def main():
    p, shots, jobs, out = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    priors = sys.argv[5] if len(sys.argv) > 5 else None
    os.makedirs(out, exist_ok=True)
    todo = [(k, b) for k in range(-(-shots // SLICE)) for b in "ZX" if not os.path.exists(f"{out}/t{k}_{b}.npz")]
    print(f"{VARIANT} d{D} p{p:g}: {len(todo)} slices to do", flush=True)
    t0 = time.time(); B = 2 * jobs; tot = len(todo)
    for i in range(0, tot, B):
        batch = todo[i:i + B]
        res = Parallel(n_jobs=jobs)(delayed(run_slice)(b, min(SLICE, shots - k * SLICE), p, 9000 + SEED_OFF + k, priors) for k, b in batch)
        for (k, b), r in zip(batch, res):
            np.savez(f"{out}/t{k}_{b}.tmp.npz", **r); os.replace(f"{out}/t{k}_{b}.tmp.npz", f"{out}/t{k}_{b}.npz")
        el = (time.time() - t0) / 60; dn = i + len(batch)
        print(f"progress {dn}/{tot} slices, {el:.1f} min, eta {el / dn * (tot - dn):.1f} min", flush=True)


if __name__ == "__main__":
    main()
