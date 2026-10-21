"""Monte Carlo of one distillation round under twirled depolarizing noise (cross-check).

Samples Z(j) with each j_l nonzero w.p. ``eps`` (value uniform on F_p^*), post-selects on
G0 j = 0, and records per-output (H1 j != 0) and block (any) logical errors.  This is the
independent sampler that the exact ``qmsd.performance`` engine is validated against, and the
analytic-side baseline for noisy-circuit (sdim) studies.
"""
from __future__ import annotations

import math

import numpy as np


def simulate(G0, H1, p, eps, shots, seed=0, batch=200_000) -> dict:
    """Estimate P_acc, eps_out (per-output average) and eps_out (block) with 1-sigma errors."""
    G0 = np.asarray(G0, dtype=np.int64) % p
    H1 = np.asarray(H1, dtype=np.int64) % p
    k, n = H1.shape
    rng = np.random.default_rng(seed)
    acc = hits = blk = 0
    done = 0
    while done < shots:
        b = min(batch, shots - done)
        j = np.where(rng.random((b, n)) < eps, rng.integers(1, p, size=(b, n)), 0)
        ok = ~np.any((j @ G0.T) % p, axis=1)
        logical = ((j[ok] @ H1.T) % p) != 0
        acc += int(ok.sum())
        hits += int(logical.sum())
        blk += int(logical.any(axis=1).sum())
        done += b

    def frac(x, m):
        f = x / m if m else float("nan")
        return f, math.sqrt(f * (1 - f) / m) if m else float("nan")

    pa, pa_err = frac(acc, shots)
    eb, eb_err = frac(blk, acc)
    ea = hits / (k * acc) if acc else float("nan")
    # per-output indicators within a shot are correlated; bound the error by the block-style binomial
    ea_err = math.sqrt(max(ea * (1 - ea), 0.0) / acc) * math.sqrt(k) if acc else float("nan")
    return {"shots": shots, "accepted": acc, "P_acc": pa, "P_acc_err": pa_err,
            "eps_out_avg": ea, "eps_out_avg_err": ea_err,
            "eps_out_block": eb, "eps_out_block_err": eb_err}
