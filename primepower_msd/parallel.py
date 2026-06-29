"""Process-based parallelism for the embarrassingly-parallel ring-code searches.

Each search is a sweep of independent random trials accumulating into a hit set, so it parallelizes
cleanly across cores (same pattern as `qmsd.search`'s joblib `n_jobs`).  Workers are top-level
functions taking (seed, trials); they return picklable summaries that the driver merges.  The GIL is
irrelevant because we use process workers (loky backend).
"""

from __future__ import annotations

import os

# CPU hygiene: keep each worker process single-threaded in BLAS so n_jobs processes do not
# oversubscribe the cores (set before numpy/joblib import their backends; loky inherits the env).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from joblib import Parallel, delayed


def default_jobs() -> int:
    """Conservative core count: respect QMSD_JOBS, else leave ~1/3 of the cores free.

    Avoids pinning the whole machine (the user asked to stay resource-careful).
    """
    env = os.environ.get("QMSD_JOBS")
    if env:
        return max(1, int(env))
    cpu = os.cpu_count() or 4
    return max(1, cpu - max(2, cpu // 3))     # e.g. 20 cores -> 14 workers, 6 free


def run_chunks(worker, n_chunks: int, trials_per_chunk: int, n_jobs: int | None = None,
               base_seed: int = 20_240_629):
    """Run `worker(seed, trials_per_chunk)` over `n_chunks` distinct seeds in parallel; return results.

    Total trials = n_chunks * trials_per_chunk, spread across n_jobs processes.  Seeds are spaced by a
    large prime so the per-worker RNG streams do not overlap.
    """
    if n_jobs is None:
        n_jobs = default_jobs()
    out = Parallel(n_jobs=n_jobs, backend="loky", verbose=0)(
        delayed(worker)(base_seed + i * 1_000_003, trials_per_chunk) for i in range(n_chunks)
    )
    return [r for r in out if r is not None]
