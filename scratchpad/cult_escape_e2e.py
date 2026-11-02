"""End-to-end cultivation -> escape at d = 4: the residual frames of accepted cultivation shots (growth 3 -> 5, r1c2) are
injected into the escape (graft_zd.Escape, recipe one+gap, gap cut 15) at the end of its noiseless colour-state
preparation, each frame into both the Z-basis run (decodes its X part: Z_L readout error a') and the X-basis run (its Z
part: X_L readout error b'), with independent escape noise.  Kept: both runs keep it and both gaps >= 15.  The residual
is harmless iff it is a power of the cultivated state's check X_L Z_L^sigma, i.e. b' = sigma a' (mod d) (calibrated by
injecting X_L, Z_L, X_L Z_L).  Frames are taken without gauge randomisation (bare errors).  Trivial frames see the escape's
own floor (MEAS_DISTILL 5h), so only nontrivial frames are run.  One json line per configuration ->
runs/cultivation/cult_escape_e2e.jsonl.

    PYTHONPATH=. python scratchpad/cult_escape_e2e.py p shots jobs config[,config...]
    config: serial | serial-hooks | par-hooks | par-hooks-f0.1 | par-hooks-f1
"""
import json
import sys
import time

import numpy as np
from joblib import Parallel, delayed

from primepower_msd.cult_growth import GrowCult
from primepower_msd.cultivation import color_code_plaquettes
from primepower_msd.graft_run import VARIANTS
from primepower_msd.graft_zd import Escape, load_graft
from primepower_msd.zdsim import FrameSim, add_tick_idle

CONFIGS = {"serial": ({}, 0.0), "serial-hooks": ({"hooks": True}, 0.0),
           "par-hooks": ({"hooks": True, "parallel": True, "final_block": True}, 0.0),
           "par-hooks-f0.1": ({"hooks": True, "parallel": True, "final_block": True}, 0.1),
           "par-hooks-f1": ({"hooks": True, "parallel": True, "final_block": True}, 1.0)}
CUT = 15


def cult_chunk(cfg, p, shots, seed, batch=20_000):
    kw, f = CONFIGS[cfg]
    g = GrowCult(4, p, 1, 1, rounds3=1, checks3=2, rounds5=1, checks5=2, **kw)
    if f:
        g.c = add_tick_idle(g.c, f * p)
    rng = np.random.default_rng(seed)
    acc, frames, lookup_fail, done = 0, [], 0, 0
    while done < shots:
        b = min(batch, shots - done)
        s = FrameSim(g.c, b, rng, gauge=False).run()       # bare error frames: no stabiliser / logical-stabiliser gauge
        ok = ~np.any(s.detector_values(), axis=1)
        x, z = s.x[ok][:, g.data5], s.z[ok][:, g.data5]
        fl, _ = g.logical_failures(x, z)
        lookup_fail += int(fl.sum()); acc += int(ok.sum())
        nt = np.any(x != 0, axis=1) | np.any(z != 0, axis=1)
        frames.append(np.stack([x[nt], z[nt]], axis=-1))
        done += b
    return acc, lookup_fail, np.concatenate(frames)


def qudit_map():
    """cultivation data index -> position in the escape's colour_q (plaquette-structure isomorphism; tau agrees)."""
    G = load_graft(4)
    col = G["codes"]["colour"]
    cq = sorted(set(q for ch in col for q in ch.support))
    Xs = {frozenset(ch.support) for ch in col if ch.typ == "X"}
    H, _ = color_code_plaquettes(5)
    Pc = [frozenset(np.nonzero(H[r])[0]) for r in range(len(H))]
    inc = [[r for r in range(len(H)) if H[r, i]] for i in range(H.shape[1])]
    order = sorted(range(H.shape[1]), key=lambda i: -len(inc[i]))
    qmap = {}

    def bt(k):
        if k == len(order):
            return all(frozenset(qmap[i] for i in P) in Xs for P in Pc)
        i = order[k]
        for e in cq:
            if e in qmap.values():
                continue
            qmap[i] = e
            if all(any(frozenset(qmap[j] for j in Pc[r] if j in qmap) <= X and len(X) == len(Pc[r]) for X in Xs) for r in inc[i]) and bt(k + 1):
                return True
            del qmap[i]
        return False
    assert bt(0)
    return np.array([cq.index(qmap[i]) for i in range(H.shape[1])])


def escape_injected(basis, frames, p, seed):
    ex = Escape(4, basis, **VARIANTS["one+gap"])
    r = ex.run(p, len(frames), seed=seed, decoder="peel-gap", inject=frames, per_shot=True)
    return r["shot_keep"], r["shot_flip"], r["shot_gap"]


def main():
    p = float(sys.argv[1]); shots = int(sys.argv[2]); jobs = int(sys.argv[3])
    perm = qudit_map()
    for cfg in sys.argv[4].split(","):
        t = time.time()
        nch = max(jobs, shots // 250_000)
        res = Parallel(n_jobs=jobs)(delayed(cult_chunk)(cfg, p, shots // nch, 1000 + k) for k in range(nch))
        acc = sum(r[0] for r in res); lookup = sum(r[1] for r in res)
        F = np.concatenate([r[2] for r in res])
        inj = np.zeros((len(F), 19, 2), dtype=np.int64)
        inj[:, perm, :] = F                                     # cultivation data i -> escape colour position perm[i]
        (kz, az, gz), (kx, bx, gx) = Parallel(n_jobs=2)(delayed(escape_injected)(b, inj, p, 7 + i) for i, b in enumerate("ZX"))
        sigma = 1
        keep = kz & kx
        keep_cut = keep & (np.nan_to_num(gz, nan=0) >= CUT) & (np.nan_to_num(gx, nan=0) >= CUT)
        bad = (bx - sigma * az) % 4 != 0
        out = dict(config=cfg, p=p, shots=shots, accepted=acc, accept=acc / shots, lookup_fails=lookup, nontrivial=len(F),
                   kept=int(keep.sum()), fails=int((keep & bad).sum()), kept_cut=int(keep_cut.sum()),
                   fails_cut=int((keep_cut & bad).sum()), cut=CUT,
                   e2e_fails_per_accepted=int((keep_cut & bad).sum()) / max(acc, 1),
                   e2e_fails_per_accepted_nocut=int((keep & bad).sum()) / max(acc, 1),
                   extra_discard=int((~keep_cut).sum()) / max(acc, 1), seconds=round(time.time() - t, 1))
        with open("runs/cultivation/cult_escape_e2e.jsonl", "a") as f:
            f.write(json.dumps(out) + "\n")
        print(json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
