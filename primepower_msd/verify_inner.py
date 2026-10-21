"""Exact re-check of a stored inner-search code: balance, normal (hosting) basis and Z_d distance.

    PYTHONPATH=. python -m primepower_msd.verify_inner runs/meas_distill/inner7/k40_s0.json 45 --d 4 8
"""
import argparse
import json
import time
import numpy as np

from .meas_distill import analyze, t_per_output, gamma_hhpw


def load(path, k):
    e = json.load(open(path))["found"][str(k)]
    S = np.array([[int(c) for c in row] for row in e["rows"]], dtype=np.int64)
    m = len(e["ell"])
    pts = np.array([[(v >> j) & 1 for j in range(m)] for v in e["support"]], dtype=np.int64)
    tau = 1 - 2 * ((pts @ np.array(e["ell"])) % 2)
    return S, tau, e


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("k", type=int)
    ap.add_argument("--d", type=int, nargs="+", default=[4])
    ap.add_argument("--Dmax", type=int, default=6)
    a = ap.parse_args(argv)
    S, tau, e = load(a.path, a.k)
    n = S.shape[1]
    for d in a.d:
        t = time.time()
        r = analyze(S, tau, d, a.Dmax, label=f"{a.path}:{a.k}")
        print(json.dumps({kk: v for kk, v in r.items() if not isinstance(v, (list, np.ndarray))}, default=str),
              f"{time.time() - t:.0f} s", flush=True)
    D = 7
    print("n", n, "k", a.k, "T/out(order 7)", round(t_per_output(n, a.k, D), 3), "gamma", round(gamma_hhpw(n, a.k, D), 4))


if __name__ == "__main__":
    main()
