"""MITM right-hand-scan throughput benchmark (probes / s / thread) on production-shaped instances.

Usage (from the worktree root, NUMBA_NUM_THREADS=1 for per-thread numbers):

    ulimit -v 4000000; OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONPATH=. \
        python scratchpad/bench_mitm/bench.py [--old] [--new] [--inst NAME ...] [--split S ...]

Instances (random seeded matrices; only the SHAPE matters for throughput):
    A  p=7  r=31 n=400    "class block" of the p=7 codes      -> 2+3 split (a=2 table, 2.9e6)
    A2 p=7  r=31 n=200    same rows, n halved so the a=3 normalised table (4.7e7, 1 GB) fits 3 GB
                                                              -> 3+3 and 3+4 splits
    B  p=23 r=15 n=420    the p=23 instance                   -> 2+3 split (a=2 table, 4.3e7, 0.7 GB)
    B2 p=23 r=15 n=60     p=23 with an a=3 normalised table (1.7e7, 0.3 GB) -> 3+3 and 3+4
    C  p=7  r=56 n=1746   the full p=7 code, 2+3 split only (a=2 table 5.5e7, 0.9 GB), scan
                          limited to a slice of leading columns

Every split is timed on a slice [i_lo, i_hi) of leading columns (the kernels prange over the
leading column; the old kernels were copied with that bound added, see old_kernels.py).  The
probe count of a slice is exact: sum_i C(n-1-i, b-1) (p-1)^{b-1 or b}.  Old and new kernels are
run on the same slice and their results compared.
"""
import argparse
import math
import os
import resource
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from qmsd import mindist_nb as mnb                      # noqa: E402
from qmsd.mindist_nb import _powers_u64                 # noqa: E402
import old_kernels as old                                # noqa: E402

INSTANCES = {
    "A":  (7, 31, 400),
    "A2": (7, 31, 200),
    "B":  (23, 15, 420),
    "B2": (23, 15, 60),
    "C":  (7, 56, 1746),
}
# split -> (instances, slice of leading columns)
SPLITS = {
    "w5_2+3":  {"A": (0, 40), "B": (0, 6), "C": (0, 30)},
    "c5_2+3":  {"A": (0, 40), "B": (0, 6), "C": (0, 30)},
    "w6_3+3":  {"A2": (0, 12), "B2": (0, 20)},
    "c6_3+3":  {"A2": (0, 12), "B2": (0, 20)},
    "w7_3+4":  {"A2": (150, 170), "B2": (20, 40)},
    "c7_3+4":  {"A2": (150, 170), "B2": (20, 40)},
}
# old kernels are ~10-100x slower: give them a shorter slice (same per-probe work)
OLD_SLICES = {
    "w5_2+3":  {"A": (0, 4), "B": (0, 1), "C": (1500, 1506)},
    "c5_2+3":  {"A": (0, 4), "B": (0, 1), "C": (1500, 1506)},
    "w6_3+3":  {"A2": (0, 2), "B2": (0, 4)},
    "c6_3+3":  {"A2": (0, 2), "B2": (0, 4)},
    "w7_3+4":  {"A2": (165, 170), "B2": (36, 40)},
    "c7_3+4":  {"A2": (165, 170), "B2": (36, 40)},
}


def nprobes(split, p, n, lo, hi):
    if split.endswith("2+3"):           # right triples, first coeff 1
        return sum(math.comb(n - 1 - i, 2) for i in range(lo, hi)) * (p - 1) ** 2
    if split.endswith("3+3"):           # right triples, all patterns
        return sum(math.comb(n - 1 - i, 2) for i in range(lo, hi)) * (p - 1) ** 3
    return sum(math.comb(n - 1 - i, 3) for i in range(lo, hi)) * (p - 1) ** 4


def instance(name):
    p, r, n = INSTANCES[name]
    rng = np.random.default_rng(12345 + len(name))
    H = rng.integers(0, p, size=(r, n)).astype(np.int64)
    H1 = rng.integers(0, p, size=(3, n)).astype(np.int64)
    return p, np.ascontiguousarray(H), np.ascontiguousarray(H1)


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def table_bytes(*arrs):
    return sum(a.nbytes for a in arrs if a is not None)


def timed(fn, *args):
    t0 = time.perf_counter()
    out = fn(*args)
    return out, time.perf_counter() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", action="store_true")
    ap.add_argument("--new", action="store_true")
    ap.add_argument("--inst", nargs="*", default=list(INSTANCES))
    ap.add_argument("--split", nargs="*", default=list(SPLITS))
    ap.add_argument("--same-slice", action="store_true", help="time old on the new slices too")
    ap.add_argument("--merge", type=int, default=0, metavar="LOG2CAP",
                    help="also time the sort+merge-join 2+3 kernel with a 2^LOG2CAP buffer")
    args = ap.parse_args()
    if not (args.old or args.new):
        args.old = args.new = True
    nthreads = int(os.environ.get("NUMBA_NUM_THREADS", "0") or 0)
    print(f"NUMBA_NUM_THREADS={nthreads}  rss_start={rss_gb():.2f} GB")
    rows = []
    for name in args.inst:
        p, H, H1 = instance(name)
        r, n = H.shape
        powers = _powers_u64(p, r)
        want = [s for s in args.split if name in SPLITS[s]]
        if not want:
            continue
        print(f"\n=== instance {name}: p={p} r={r} n={n}")
        tables = {}
        if any(s.endswith("2+3") for s in want):
            (lc, ls), tb = timed(mnb._build_left, H, p, 2, powers)
            print(f"  a=2 table: N={lc.shape[0]:.3g} {table_bytes(lc, ls)/1e9:.2f} GB  build {tb:.1f}s")
            tables["a2"] = (lc, ls, None, tb)
            del lc, ls
            if "c5_2+3" in want:
                (lc, ls, lf), tb = timed(mnb._build_left2_coef, H, p, powers)
                tables["a2c"] = (lc, ls, lf, tb)
                del lc, ls, lf
        if any(s.endswith("3+3") or s.endswith("3+4") for s in want):
            (lc, ls, lf), tb = timed(mnb._build_left3_norm_coef, H, p, powers)
            print(f"  a=3 norm table: N={lc.shape[0]:.3g} {table_bytes(lc, ls, lf)/1e9:.2f} GB  build {tb:.1f}s")
            tables["a3"] = (lc, ls, lf, tb)
            del lc, ls, lf
        idx = {}
        if args.new:
            for key, (lc, ls, lf, tb) in tables.items():
                (bstart, shift), ti = timed(mnb._left_index, lc)
                idx[key] = (bstart, shift)
                print(f"  index[{key}]: K={64-shift} {bstart.nbytes/1e6:.0f} MB  build {ti:.2f}s "
                      f"(= {100*ti/max(tb,1e-9):.1f}% of table build)")
        for split in want:
            lo, hi = SPLITS[split][name]
            olo, ohi = SPLITS[split][name] if args.same_slice else OLD_SLICES[split][name]
            res = {}
            if args.old:
                np_old = nprobes(split, p, n, olo, ohi)
                res["old"] = run(split, "old", H, H1, p, powers, tables, idx, olo, ohi)
                res["old"] = (res["old"][0], res["old"][1], np_old, (olo, ohi))
            if args.new:
                np_new = nprobes(split, p, n, lo, hi)
                res["new"] = run(split, "new", H, H1, p, powers, tables, idx, lo, hi)
                res["new"] = (res["new"][0], res["new"][1], np_new, (lo, hi))
                if args.old and args.same_slice:
                    check_same(split, res["old"][0], res["new"][0])
            if args.merge and split == "w5_2+3":
                import merge_kernels
                lc, ls, _, _ = tables["a2"]
                out, dt = timed(merge_kernels.scan_b3_merge, H, p, powers, lc, ls, lo, hi,
                                1 << args.merge)
                res["merge"] = (out, dt, nprobes(split, p, n, lo, hi), (lo, hi))
                if args.new:
                    assert out == res["new"][0], ("merge != index", out, res["new"][0])
            for which, (out, dt, npb, sl) in res.items():
                rate = npb / dt / max(nthreads, 1)
                rows.append((name, split, which, sl, npb, dt, rate))
                print(f"  {split:8s} {which}: cols {sl} probes={npb:.3g} time={dt:.2f}s "
                      f"-> {rate:.3g} probes/s/thread   out={summ(out)}")
        del tables, idx
        print(f"  rss_max so far {rss_gb():.2f} GB")
    print("\n| instance | split | path | cols | probes | s | probes/s/thread |")
    print("|---|---|---|---|---|---|---|")
    for name, split, which, sl, npb, dt, rate in rows:
        print(f"| {name} | {split} | {which} | {sl} | {npb:.3g} | {dt:.2f} | {rate:.3g} |")


def summ(out):
    if isinstance(out, tuple):
        return tuple(int(np.asarray(o).sum()) for o in out)
    return out


def check_same(split, a, b):
    if isinstance(a, tuple):
        for x, y in zip(a, b):
            assert np.array_equal(np.asarray(x), np.asarray(y)), split
    else:
        assert a == b, (split, a, b)
    print(f"  {split}: old == new OK")


def run(split, which, H, H1, p, powers, tables, idx, lo, hi):
    n = H.shape[1]
    if split == "w5_2+3":
        lc, ls, _, _ = tables["a2"]
        if which == "old":
            return timed(old.scan_b3_old, H, p, powers, lc, ls, lo, hi)
        bs, sh = idx["a2"]
        return timed(mnb._scan_right_b3_idx, H, p, powers, lc, ls, bs, sh, lo, hi)
    if split == "c5_2+3":
        lc, ls, lf, _ = tables["a2c"]
        if which == "old":
            return timed(old.census5_old, H, H1, p, powers, lc, ls, lf, lo, hi)
        bs, sh = idx["a2c"]
        return timed(mnb._census_weight5_idx, H, H1, p, powers, lc, ls, lf, bs, sh, lo, hi)
    lc, ls, lf, _ = tables["a3"]
    if split == "w6_3+3":
        found = np.zeros(1, dtype=np.int64); wit = np.full(6, -1, dtype=np.int32)
        if which == "old":
            return timed(old.scan3_full_old, H, p, powers, lc, ls, found, wit, lo, hi)
        bs, sh = idx["a3"]
        return timed(mnb._scan_right3_idx, H, p, powers, lc, ls, bs, sh, found, wit, lo, hi)
    if split == "c6_3+3":
        if which == "old":
            return timed(old.census6_old, H, H1, p, powers, lc, ls, lf, lo, hi)
        bs, sh = idx["a3"]
        return timed(mnb._census_weight6_idx, H, H1, p, powers, lc, ls, lf, bs, sh, lo, hi)
    if split == "w7_3+4":
        found = np.zeros(1, dtype=np.int64); wit = np.full(7, -1, dtype=np.int32)
        if which == "old":
            return timed(mnb._scan_right4_full, H, p, powers, lc, ls, lo, hi, found, wit)
        bs, sh = idx["a3"]
        return timed(mnb._scan_right4_idx, H, p, powers, lc, ls, bs, sh, found, wit, lo, hi)
    if split == "c7_3+4":
        if which == "old":
            return timed(mnb._census_weight7_3x4, H, H1, p, powers, lc, ls, lf, lo, hi)
        bs, sh = idx["a3"]
        return timed(mnb._census_weight7_idx, H, H1, p, powers, lc, ls, lf, bs, sh, lo, hi)
    raise ValueError(split)


if __name__ == "__main__":
    main()
