"""Keep rules for the f = 1 escape at one cultivation input: every mask-kept shot of both basis runs is decoded (Escape.run
flag_cert), so failures on shots where the certified-pattern discards of the two bases disagree are sampled.  Residual
frames of cultivation and N_CLEAN zero frames (the escape's floor) go through the same rules with common noise per slice.
Rules: cert = none (gap only), either (discard on a hit in either basis), z, x (hits of one basis only); cuts = (cut_Z, cut_X).
p_total = (fails_res + w fails_clean) / (kept_res + w kept_clean), w = (accepted - nontrivial) / N_CLEAN;
failure b' != a' (mod d) as in cult_escape_e2e.  One json line per variant -> runs/surface/escape_keep_rules.jsonl.

    PYTHONPATH=.:scratchpad python scratchpad/escape_keep_rules.py p shots n_clean jobs cult_config variant[,variant...]
"""
import json
import sys
import time

import numpy as np
from joblib import Parallel, delayed

from cult_escape_e2e import cult_chunk, qudit_map
from primepower_msd.graft_run import VARIANTS
from primepower_msd.graft_zd import Escape

SLICE = 50_000
CUTS = ((0, 0), (10, 10), (12, 12), (13, 13), (14, 13), (15, 13), (14, 14), (15, 15), (16, 16))
_cache = {}


def run_slice(variant, basis, frames, n, p, seed):
    """frames None: n clean frames."""
    if (variant, basis) not in _cache:
        _cache[(variant, basis)] = Escape(4, basis, **VARIANTS[variant])
    inj = np.zeros((n, 19, 2), dtype=np.int64) if frames is None else frames
    r = _cache[(variant, basis)].run(p, n, seed=seed, decoder="peel-gap", inject=inj, per_shot=True, chunk_seeded=True, flag_cert=True)
    return r["shot_keep"], r["shot_cert"], r["shot_flip"].astype(np.int8), np.nan_to_num(r["shot_gap"], nan=0).astype(np.float32)


def run_all(variant, F, n_clean, p, jobs):
    """(slices of residual frames, then clean) -> arrays by basis; both bases of a slice share a seed."""
    tasks = [(F[i:i + SLICE], len(F[i:i + SLICE])) for i in range(0, len(F), SLICE)]
    tasks += [(None, min(SLICE, n_clean - i)) for i in range(0, n_clean, SLICE)]
    out = Parallel(n_jobs=jobs)(delayed(run_slice)(variant, b, fr, n, p, 7000 + k) for k, (fr, n) in enumerate(tasks) for b in "ZX")
    return [[np.concatenate([out[2 * k + i][c] for k in range(len(tasks))]) for c in range(4)] for i in range(2)]


def analyse(Z, X, nt, acc, shots, n_clean):
    (kz, cz, az, gz), (kx, cx, bx, gx) = Z, X
    bad = (bx.astype(int) - az.astype(int)) % 4 != 0
    mask = kz & kx
    w = (acc - nt) / max(n_clean, 1)
    rows = {}
    for rule, hit in (("none", np.zeros_like(cz)), ("either", cz | cx), ("z", cz), ("x", cx)):
        for c_z, c_x in CUTS:
            k = mask & ~hit & (gz >= c_z) & (gx >= c_x)
            kr, kc = int(k[:nt].sum()), int(k[nt:].sum())
            fr, fc = int((k & bad)[:nt].sum()), int((k & bad)[nt:].sum())
            n = kr + w * kc
            rows[f"{rule}:{c_z},{c_x}"] = dict(kept_res=kr, fails_res=fr, kept_clean=kc, fails_clean=fc, keep_clean=kc / max(n_clean, 1),
                                               p_res=fr / n if n else None, p_clean=w * fc / n if n else None, p_total=(fr + w * fc) / n if n else None,
                                               attempts_per_kept=shots / n if n else None)
    return dict(mask_disagree=float((kz != kx).mean()), cert_z=float((cz & mask).sum() / mask.sum()), cert_x=float((cx & mask).sum() / mask.sum()),
                cert_disagree=float(((cz != cx) & mask).sum() / mask.sum()), rules=rows)


def main():
    p = float(sys.argv[1]); shots = int(sys.argv[2]); n_clean = int(sys.argv[3]); jobs = int(sys.argv[4]); cfg = sys.argv[5]
    t = time.time()
    perm = qudit_map()
    nch = max(jobs, shots // 250_000)
    res = Parallel(n_jobs=jobs)(delayed(cult_chunk)(cfg, p, shots // nch, 1000 + k) for k in range(nch))
    acc = sum(r[0] for r in res)
    Fc = np.concatenate([r[2] for r in res])
    F = np.zeros((len(Fc), 19, 2), dtype=np.int64)
    F[:, perm, :] = Fc
    print("cultivation", cfg, "accepted", acc, "nontrivial", len(F), "s", round(time.time() - t), flush=True)
    for variant in sys.argv[6].split(","):
        t = time.time()
        Z, X = run_all(variant, F, n_clean, p, jobs)
        out = dict(config=cfg, variant=variant, p=p, shots=shots, accepted=acc, nontrivial=len(F), n_clean=n_clean,
                   **analyse(Z, X, len(F), acc, shots, n_clean), seconds=round(time.time() - t, 1))
        with open("runs/surface/escape_keep_rules.jsonl", "a") as f:
            f.write(json.dumps(out) + "\n")
        print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
