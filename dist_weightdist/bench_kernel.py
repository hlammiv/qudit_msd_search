"""Single-core throughput probe for the gray-code incremental weight-update inner loop.

This measures the FLOOR rate the distributed plan rests on: a numba-njit scalar emulation of
one gray-code step -- each consecutive message differs in one trit, so the codeword changes by
+/- one generator row, and we rescan ONLY that row's support to update the Hamming weight.
A production trit-packed C/Rust kernel (2 bits/trit, 32 trits/64-bit word, popcount-style
weight delta) does the per-step support scan word-parallel and is ~an order of magnitude
faster; the numba scalar number here is a conservative lower bound the correctness angle cites.

Run:  python -m dist_weightdist.bench_kernel
"""
from __future__ import annotations

import time

import numpy as np
from numba import njit


@njit(cache=True)
def _gray_steps(cw, support, rowvals, weight, nsteps):
    """Apply nsteps incremental gray-code updates over a row's support; return final weight."""
    w = weight
    s = len(support)
    for _ in range(nsteps):
        for t in range(s):
            p = support[t]
            old = cw[p]
            new = old + rowvals[t]
            if new >= 3:
                new -= 3
            if old != 0:
                w -= 1
            if new != 0:
                w += 1
            cw[p] = new
    return w


def measure(n: int, support_frac: float = 0.6, nsteps: int = 5_000_000) -> float:
    rng = np.random.default_rng(0)
    cw = rng.integers(0, 3, size=n).astype(np.int8)
    s = max(1, int(support_frac * n))
    support = np.arange(s, dtype=np.int64)
    rowvals = rng.integers(1, 3, size=s).astype(np.int8)
    w = int(np.count_nonzero(cw))
    _gray_steps(cw.copy(), support, rowvals, w, 10)  # JIT warmup
    t0 = time.time()
    _gray_steps(cw, support, rowvals, w, nsteps)
    return nsteps / (time.time() - t0)


def project(rate_per_core: float, threads: int = 52, K: int = 26) -> None:
    total = 3**K
    secs = total / (rate_per_core * threads)
    print(
        f"  3^{K} = {total:.3e} codewords / ({rate_per_core:,.0f}/s/core x {threads} threads)"
        f" = {secs:,.0f} s = {secs / 3600:.1f} h"
    )


if __name__ == "__main__":
    for n in (667, 1939):
        r = measure(n)
        print(f"n={n}: numba scalar gray-step floor ~ {r:,.0f} codewords/s/core")
    print("Projection for 3^26 on 52 threads at the scalar floor and trit-packed estimates:")
    base = measure(1939)
    for mult, tag in ((1, "scalar floor"), (3, "trit-packed ~3x"), (6, "trit-packed ~6x")):
        print(f" [{tag}]")
        project(base * mult)
