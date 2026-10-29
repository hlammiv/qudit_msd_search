"""Escape study driver.  python -m primepower_msd.graft_run select basis     -> hyperedge share / discard per colour combo
                        python -m primepower_msd.graft_run run basis cX cZ shots p1,p2,... [jobs] [modes]  -> escape vs final-memory vs D=15
                        python -m primepower_msd.graft_run variants - shots p1,p2,... jobs [names] [cert_jobs] [d] [replicas]  -> variants, both bases, joint postselection"""
from __future__ import annotations

import json
import sys
import time

from joblib import Parallel, delayed


def _select(basis, cX, cZ):
    from .graft_zd import Escape
    t = time.time()
    try:
        ex = Escape(4, basis, post_colours=(cX, cZ))
    except RuntimeError as e:
        return dict(cX=cX, cZ=cZ, error=str(e)[:200])
    r = ex.run(1e-3, 2000, seed=1)
    r.update(cX=cX, cZ=cZ, build_s=round(time.time() - t, 1), n_post=int(len(ex.post_idx)))
    return r


def _point(mode, basis, p, shots, cX, cZ):
    from .graft_zd import Escape
    from .surface_zd import Experiment
    t = time.time()
    if mode == "surface15":
        ex = Experiment(15, 4, basis, rounds=10)
        r = ex.run(p, shots, seed=hash((mode, basis, p)) % 2**31)
        r["mode"] = mode; r["discard"] = 0.0
    else:
        ex = Escape(4, basis, mode=mode, post_colours=(cX, cZ) if mode == "escape" else (None, None))
        r = ex.run(p, shots, seed=hash((mode, basis, p)) % 2**31)
    r["seconds"] = round(time.time() - t, 1)
    return r


VARIANTS = {
    "one":        dict(post_colours=(2, 0)),
    "all":        dict(post_colours=({0, 1, 2}, {0, 1, 2})),
    "one+flags":  dict(post_colours=(2, 0), flag_weight=5),
    "one+pairs":  dict(post_colours=(2, 0)),          # + second-order certification on the colour region
    "all+pairs":  dict(post_colours=({0, 1, 2}, {0, 1, 2})),
    "one+pairs-cap":   dict(post_colours=(2, 0)),     # colour pairs, degenerate first-order mass capped at 5 (0.5% at 1e-3)
    "one+pairs2":      dict(post_colours=(2, 0)),     # colour x everything pairs
    "one+pairs2-cap":  dict(post_colours=(2, 0)),
    "one+hyper":       dict(post_colours=(2, 0)),     # BP+LSD hypergraph decoder (hyperdecoder.py) instead of 2-adic matching
    "all+hyper":       dict(post_colours=({0, 1, 2}, {0, 1, 2})),
    "one+gap":         dict(post_colours=(2, 0)),     # 2-adic matcher with the complementary gap recorded (gap_curve)
    "one+pairs-cap+gap": dict(post_colours=(2, 0)),
}
DECODER = {"one+hyper": "hyper", "all+hyper": "hyper", "one+gap": "peel-gap", "one+pairs-cap+gap": "peel-gap"}
CERT = {  # variant -> (scope, cache tag, degenerate cap)
    "one+pairs": ("colour", "one", None), "all+pairs": ("colour", "all", None),
    "one+pairs-cap": ("colour", "one", 2.0), "one+pairs2": ("colour-all", "one_all", None),
    "one+pairs2-cap": ("colour-all", "one_all", 2.0), "one+pairs-cap+gap": ("colour", "one", 2.0),
}


def _variant(name, basis, p, shots, d=4, replica=0):
    from .graft_zd import Escape
    t = time.time()
    ex = Escape(d, basis, **VARIANTS[name])
    stats = {}
    if name in CERT:
        scope, tag, cap = CERT[name]
        stats = ex.certify_pairs(p, n_jobs=int(sys.argv[7]) if len(sys.argv) > 7 else 1, scope=scope, degenerate_cap=cap,
                                 cache=f"runs/surface/cert_pairs_{tag}_{basis}_{p:g}{'' if d == 4 else f'_d{d}'}.json")
    r = ex.run(p, shots, seed=hash((name, basis, p, replica)) % 2**31, decoder=DECODER.get(name, "peel"))
    r.update(variant=name, seconds=round(time.time() - t, 1), n_post=int(len(ex.mask)), pair_stats=stats, replica=replica)
    return {k: (v.item() if hasattr(v, "item") else v) for k, v in r.items()}


def main():
    cmd = sys.argv[1]; basis = sys.argv[2]
    if cmd == "variants":        # variants shots p1,p2,... jobs [names]
        shots = int(sys.argv[3]); ps = [float(x) for x in sys.argv[4].split(",")]; jobs = int(sys.argv[5])
        names = sys.argv[6].split(",") if len(sys.argv) > 6 else ["one", "all", "one+flags"]
        d = int(sys.argv[8]) if len(sys.argv) > 8 else 4                 # argv[7] = certification jobs
        replicas = int(sys.argv[9]) if len(sys.argv) > 9 else 1         # independent seeds per (variant, basis, p)
        tasks = [(nm, b, p, shots, d, k) for nm in names for b in "XZ" for p in ps for k in range(replicas)]
        out = "runs/surface/graft_variants.jsonl" if d == 4 else f"runs/surface/graft_variants_d{d}.jsonl"
        res = []
        for r in Parallel(n_jobs=jobs, return_as="generator_unordered")(delayed(_variant)(*t) for t in tasks):
            with open(out, "a") as f:                 # one line per finished task: a crash or kill loses only the tasks in flight
                f.write(json.dumps(r) + "\n")
            res.append(r)
        for r in res:
            print(r["variant"], r["basis"], r["p"], "p_fail %.3g" % r["p_fail"], "fails", r["fails"], "kept", r["kept"],
                  "discard %.4f" % r["discard"], "bad_mass %.3f" % r["bad_pattern_mass"], "pairs", r.get("n_pair_patterns"),
                  "pair_mass %.1f" % r.get("pair_mass", 0), "degenerate", r["pair_stats"].get("degenerate"), "kept", r["pair_stats"].get("degenerate_kept"),
                  "qudits", r["n_qudits"], "s", r["seconds"])
            if r.get("gap_curve"):
                print("   gap curve (threshold, total discard, p_fail, fails):",
                      " ".join("%.0f:%.3f:%.2e:%d" % tuple(c) for c in r["gap_curve"]))
                print("   failing shots' gaps:", [round(x, 1) for x in r["fail_gaps"]])
        return
    if cmd == "select":
        combos = [(None, None)] + [(a, b) for a in range(3) for b in range(3)]
        res = Parallel(n_jobs=3)(delayed(_select)(basis, a, b) for a, b in combos)
        for r in res:
            print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items() if k in ("cX", "cZ", "hyper_share", "discard", "p_fail", "fails", "kept", "undetectable_share", "conflicts", "n_post", "error")})
        json.dump(res, open(f"runs/surface/graft_select_{basis}.json", "w"), indent=1)
    else:
        cX, cZ = int(sys.argv[3]), int(sys.argv[4]); shots = int(sys.argv[5]); ps = [float(x) for x in sys.argv[6].split(",")]
        modes = sys.argv[8].split(",") if len(sys.argv) > 8 else ("escape", "final", "surface15")
        tasks = [(m, basis, p, shots, cX, cZ) for p in ps for m in modes]
        res = Parallel(n_jobs=int(sys.argv[7]) if len(sys.argv) > 7 else 6)(delayed(_point)(*t) for t in tasks)
        with open("runs/surface/graft_d4.jsonl", "a") as f:
            for r in res:
                f.write(json.dumps(r) + "\n")
        for r in res:
            print(r["mode"], r["basis"], r["p"], "p_fail", r["p_fail"], "fails", r["fails"], "kept", r.get("kept"), "discard", round(r["discard"], 4),
                  "hyper_discard", r.get("hyper_discard"), "bad_patterns", r.get("n_bad_patterns"), "s", r.get("seconds"))


if __name__ == "__main__":
    main()
