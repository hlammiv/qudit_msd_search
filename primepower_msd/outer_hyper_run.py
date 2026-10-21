"""Driver: order-7 (s=3) outer codes of row weight w by violator-guided switching; full certificate at the end.

    PYTHONPATH=. python -m primepower_msd.outer_hyper_run --n 7000 --w 21 --seed 0 --out runs/meas_distill/outer7_w21_n7000_s0.json
"""
import argparse
import json
import time

import numpy as np

from .outer_hyper import certify, check_bits_from, construct2


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--w", type=int, default=21)
    ap.add_argument("--s", type=int, default=3)
    ap.add_argument("--D", type=int, default=7)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max_switch", type=int, default=3_000_000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    t0 = time.time()
    r = construct2(a.n, a.s, a.w, a.D, seed=a.seed, verbose=True)
    rec = {"n": a.n, "w": a.w, "s": a.s, "D": a.D, "seed": a.seed, "seconds": round(time.time() - t0)}
    if r is None:
        rec["status"] = "failed"
    else:
        bc, m, sw = r
        cb = check_bits_from(bc, m, a.w)
        ok = certify(bc, cb, a.D) < 0
        rec.update({"status": "certified" if ok else "certificate failed", "checks": int(m), "switches": int(sw),
                    "bit_checks": bc.tolist()})
    json.dump(rec, open(a.out, "w"))
    print(rec["status"], rec["seconds"], "s", flush=True)


if __name__ == "__main__":
    main()
