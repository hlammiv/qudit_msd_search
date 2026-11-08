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

from cult_escape_e2e import D, cult_chunk, qudit_map
from primepower_msd.graft_run import VARIANTS
from primepower_msd.graft_zd import Escape

SLICE = int(__import__("os").environ.get("EKR_SLICE", 50_000))
SEED_OFF = int(__import__("os").environ.get("EKR_SEED", 0))   # offset of every seed: independent replicas
CUTS = ((0, 0), (10, 10), (12, 12), (13, 13), (14, 13), (15, 13), (14, 14), (15, 15), (16, 16))
_cache = {}
SYN_CUTS = tuple((c, c) for c in range(0, 17)) + ((12, 11), (13, 12), (14, 12), (14, 13), (15, 13), (15, 14), (16, 14))


def RULES(cz, cx, only_none=False):
    rules = (("none", np.zeros_like(cz)), ("either", cz | cx), ("z", cz), ("x", cx))
    return rules[:1] if only_none else rules


def run_slice(variant, basis, frames, n, p, seed):
    """frames None: n clean frames."""
    if (variant, basis) not in _cache:
        _cache[(variant, basis)] = Escape(D, basis, **VARIANTS[variant])
    inj = np.zeros((n, 19, 2), dtype=np.int64) if frames is None else frames
    r = _cache[(variant, basis)].run(p, n, seed=seed, decoder="peel-gap", inject=inj, per_shot=True, chunk_seeded=True, flag_cert=True)
    g2 = r["shot_gap_syn"] if _cache[(variant, basis)].gap_window in ("syndrome", "closed") else r["shot_gap"]
    return (r["shot_keep"], r["shot_cert"], r["shot_flip"].astype(np.int8), np.nan_to_num(r["shot_gap"], nan=0).astype(np.float32),
            np.nan_to_num(g2, nan=0).astype(np.float32))


def run_all(variant, F, n_clean, p, jobs):
    """(slices of residual frames, then clean) -> arrays by basis; both bases of a slice share a seed."""
    tasks = [(F[i:i + SLICE], len(F[i:i + SLICE])) for i in range(0, len(F), SLICE)]
    tasks += [(None, min(SLICE, n_clean - i)) for i in range(0, n_clean, SLICE)]
    out = Parallel(n_jobs=jobs)(delayed(run_slice)(variant, b, fr, n, p, 7000 + SEED_OFF + k) for k, (fr, n) in enumerate(tasks) for b in "ZX")
    return [[np.concatenate([out[2 * k + i][c] for k in range(len(tasks))]) for c in range(5)] for i in range(2)]


def analyse(Z, X, nt, acc, shots, n_clean, cuts=CUTS, only_none=False):
    (kz, cz, az, gz), (kx, cx, bx, gx) = Z, X          # (keep, cert hit, readout, gap)
    bad = (bx.astype(int) - az.astype(int)) % D != 0
    mask = kz & kx
    w = (acc - nt) / max(n_clean, 1)
    rows = {}
    for rule, hit in RULES(cz, cx, only_none):
        for c_z, c_x in cuts:
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
    res = Parallel(n_jobs=jobs)(delayed(cult_chunk)(cfg, p, shots // nch, 1000 + SEED_OFF + k) for k in range(nch))
    acc = sum(r[0] for r in res)
    Fc = np.concatenate([r[2] for r in res])
    F = np.zeros((len(Fc), 19, 2), dtype=np.int64)
    F[:, perm, :] = Fc
    print("cultivation", cfg, "accepted", acc, "nontrivial", len(F), "s", round(time.time() - t), flush=True)
    for variant in sys.argv[6].split(","):
        t = time.time()
        Z, X = run_all(variant, F, n_clean, p, jobs)
        out = dict(config=cfg, variant=variant, p=p, d=D, seed_off=SEED_OFF, shots=shots, accepted=acc, nontrivial=len(F), n_clean=n_clean,
                   **analyse(Z[:4], X[:4], len(F), acc, shots, n_clean), seconds=round(time.time() - t, 1))
        if VARIANTS[variant].get("gap_window") in ("syndrome", "closed"):    # the same shots, kept on the syndrome-only gap
            syn = analyse(Z[:3] + Z[4:], X[:3] + X[4:], len(F), acc, shots, n_clean, SYN_CUTS, only_none=True)
            out["rules_readout_fine"] = analyse(Z[:4], X[:4], len(F), acc, shots, n_clean, SYN_CUTS, only_none=True)["rules"]
            out["rules_syndrome_gap"] = syn["rules"]
            gz, gx, sz, sx = Z[3], X[3], Z[4], X[4]
            ok = Z[0] & X[0]
            out["gap_corr"] = [float(np.corrcoef(gz[ok], sz[ok])[0, 1]), float(np.corrcoef(gx[ok], sx[ok])[0, 1])]
            out["gap_max"] = dict(readout=[float(gz.max()), float(gx.max())], syndrome=[float(sz.max()), float(sx.max())])
        with open("runs/surface/escape_keep_rules.jsonl", "a") as f:
            f.write(json.dumps(out) + "\n")
        print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
