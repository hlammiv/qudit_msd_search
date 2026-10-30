# MITM right-scan throughput: before / after

Machine: i7-12700H (P-cores 4.7 GHz, L2 1.25 MB/core, L3 24 MB), DDR5, shared box.
`NUMBA_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `ulimit -v 4000000`; numbers are probes / s for
the one thread (a probe = one right-hand (subset, coefficient pattern) hashed and looked up).
Old = the original kernels (per-probe `np.searchsorted` + two `% p` per row); new =
`method="index"` (bucket index + staged lookups + division-free hash), commit on this branch.
Each pair of rows is timed on a slice of leading columns of the SAME instance and left table
(old on a shorter slice, same per-probe work; `--same-slice` runs confirmed old == new output).

| instance | split | left table (entries, GB) | old probes/s | new probes/s | speed-up |
|---|---|---|---|---|---|
| A  p=7 31x400   | 2+3 witness (w5)   | 2.9e6, 0.05 | 1.21e6 | 8.04e6 | 6.6x |
| A  p=7 31x400   | 2+3 census (c5)    | 2.9e6, 0.05 | 1.16e6 | 6.41e6 | 5.5x |
| B  p=23 15x420  | 2+3 witness (w5)   | 4.3e7, 0.68 | 0.97e6 | 6.45e6 | 6.6x |
| B  p=23 15x420  | 2+3 census (c5)    | 4.3e7, 0.68 | 0.98e6 | 6.22e6 | 6.3x |
| C  p=7 56x1746  | 2+3 witness (w5)   | 5.5e7, 0.88 | 0.71e6 | 5.28e6 | 7.4x |
| A2 p=7 31x200   | 3+3 witness (w6)   | 4.7e7, 1.04 | 0.81e6 | 5.68e6 | 7.0x |
| A2 p=7 31x200   | 3+3 census (c6)    | 4.7e7, 1.04 | 0.86e6 | 5.28e6 | 6.2x |
| A2 p=7 31x200   | 3+4 witness (w7)   | 4.7e7, 1.04 | 0.85e6 | 5.75e6 | 6.7x |
| A2 p=7 31x200   | 3+4 census (c7)    | 4.7e7, 1.04 | 0.87e6 | 5.67e6 | 6.5x |
| B2 p=23 15x60   | 3+3 witness (w6)   | 1.7e7, 0.36 | 0.88e6 | 6.09e6 | 7.0x |
| B2 p=23 15x60   | 3+3 census (c6)    | 1.7e7, 0.36 | 0.86e6 | 5.76e6 | 6.7x |

Old A2 numbers are the solo re-run (a first pass taken while another single-thread bench ran on
the box gave 0.56-0.64e6, i.e. the old path degrades first under contention).  A2 / B2 are
the 3+3 / 3+4 shapes scaled so the normalised a=3 table fits 3 GB (the real 31x400 block needs
7.6 GB; the per-probe cost only grows with table size, so the speed-up there is, if anything,
larger -- see the diagnostic below).

## Where the time goes (2+3, p=7, r=31; `diag.py`)

| left table | old-style cost | index v1 (no staging) | index v2 (staged) |
|---|---|---|---|
| 1 MB (cache resident)  | --     | 68 ns/probe | 57 ns/probe |
| 2.8 MB (L3)            | --     | 76 ns/probe | 55 ns/probe |
| 43 MB (DRAM)           | 830 ns | 166 ns/probe | 99 ns/probe |

* The old kernel is ~800 ns/probe even when the table is nearly cache resident (instance A):
  the two 64-bit `% p` per row (62 divisions per probe at r = 31) cost as much as the binary
  search.  The division-free hash (pre-reduced scaled columns + `red`/`neg` tables, partial sums
  reused over the innermost column) brings the compute floor to ~55 ns/probe.
* The bucket index (top-K hash bits, 2^K <= N/8, i.e. ~8-16 codes per bucket) replaces the
  ~log2(N) dependent misses by one index read plus one or two adjacent table lines.  Without
  staging the misses of consecutive probes were still serialised by the data-dependent bucket
  walk (166 ns at 43 MB); staging each innermost column into a batch of <= 1024 targets --
  pass 1 hashes and reads the two index entries (independent loads), pass 2 walks the buckets
  with early exit -- lets the core overlap them (99 ns at 43 MB, ~170 ns at 1 GB tables).
* Sort + merge-join alternative (`merge_kernels.py`, argsort of 2^20 targets per thread, then a
  galloping merge against the sorted table): 4.0e6 probes/s on A, 3.5e6 on B, 3.0e6 on C, i.e.
  about half the index path (still 3-4x the old kernel, as it shares the division-free hash).  The sort itself costs ~100 ns/target and the galloping merge still touches
  ~log2(N/M) new lines per target; the index path reads strictly fewer lines and needs no
  per-thread buffer, so it is the one wired in.

## Memory and setup

| table | entries | table bytes | index K | index bytes | index build | table build (argsort) |
|---|---|---|---|---|---|---|
| A  a=2      | 2.9e6 | 46 MB  | 19 | 4 MB  | 0.01 s | 1.3 s |
| B  a=2      | 4.3e7 | 680 MB | 23 | 67 MB | 0.18 s | 18 s |
| C  a=2      | 5.5e7 | 880 MB | 23 | 67 MB | 0.21 s | 30 s |
| A2 a=3 norm | 4.7e7 | 1.04 GB | 23 | 67 MB | 0.23 s | 28 s |
| B2 a=3 norm | 1.7e7 | 360 MB | 21 | 17 MB | 0.07 s | 7 s |

Index = 8 bytes x 2^K with 2^K <= N/8: at most N bytes, 3-6% of the table (which is 16-22 B per
entry); build is a single pass over the sorted codes, ~1% of the argsort that builds the table.
Per thread the staged kernels add the partial-sum tables ((p-1)^2 r and (p-1)^3 r int64 for the
3+4 / 2+4 streams: 1.3 MB at p = 23, r = 15; 11 MB at p = 31, r = 50) and three 1024-entry
batch arrays; peak RSS of every bench run was unchanged to within 20 MB (2.10 vs 2.08 GB on A2).
Per-call overhead on tiny inputs: 0.4 ms vs 0.16 ms (5x8 matrix, 2+4 path), negligible against
the test-suite and search-loop granularity.

## Raw logs

`/tmp/claude-1000/bench_{old,new}_{A,A2,B,B2,C}.log` on the dev box (not committed);
`bench_new_*_v1.log` hold the unstaged index-v1 numbers (A 4.6e6, A2 3.0-3.2e6, B2 3.8e6).
