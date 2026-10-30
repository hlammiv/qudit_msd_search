# MITM right-scan throughput bench

Measures probes / s / thread of the meet-in-the-middle right-hand scans in `qmsd/mindist_nb.py`
on production-shaped random instances, old kernels (per-probe `np.searchsorted`) against the
bucket-indexed kernels (`method="index"`, the default since this branch).

```
ulimit -v 4000000
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONPATH=. \
    python scratchpad/bench_mitm/bench.py --old --new          # full table, ~10 min
    python scratchpad/bench_mitm/bench.py --new --inst C       # one instance
    python scratchpad/bench_mitm/bench.py --same-slice --inst A2 --split c6_3+3   # old == new check
```

* `bench.py`       driver: instances, probe counts, timings, markdown table at the end
* `old_kernels.py` frozen copies of the original kernels with an `[i_lo, i_hi)` leading-column
                   bound added (byte-for-byte loop bodies), so the slow path can be timed on a slice
* `RESULTS.md`     the numbers recorded on the development machine (i7-12700H, 1 thread)

Instances are random matrices with the shapes of the real codes (the hash / lookup cost does
not depend on the code structure): A p=7 31x400, B p=23 15x420, C p=7 56x1746 (2+3 slice only),
plus A2 (p=7 31x200) and B2 (p=23 15x60) whose normalised a=3 tables fit the 3 GB budget.
