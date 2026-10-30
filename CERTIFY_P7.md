# Certification queue for the p = 7 sweep candidates (and the p >= 17 ladder)

*Started 2026-10-29. Status is updated in place; raw outputs under `runs/certify/` (p = 7) and
`runs/m2_line_sweep/` (p >= 17). All jobs run at nice 19 on lenore and yield to the owner's
circuit-knitting runs; wall times below are for the thread counts stated.*

## What is being certified and why

The multi-seed sweep (`WEIGHTED_SUPPORTS.md`, "Multi-seed sweep") produced 79 candidates at k = 265
and 51 at k = 270 on the sphere-charge code [[2016 - k, k, d]]_7: each has plane distance >= 7 and no
class-block word of weight <= 5, and 129 of the 130 carry an explicit planar weight-7 word (so d <= 7
is already proven). What is missing is the lower bound d >= 7, which needs three scans the cheap
filters cannot do:

| stage | what it excludes | method | cost (8 threads) |
|---|---|---|---|
| block6 | weight-6 words inside the class block | 3+3 MITM on the block (`weight6_witness_3x3`, ~6 GB table) | ~40 min |
| full5 | weight <= 5 words anywhere | serial MITM d <= 4, thread-parallel 2+3 weight-5 scan on the full code (`weight5_exists_parallel`) | ~10 h |
| mixed6 | weight-6 words touching a non-class column | kernel rows force the non-class part into RM(20,4) (d = 5); its full-support weight-5/6 words are line-supported (lemma, `scratchpad/mixed_w6_check.py`); enumerate those lines inside the non-class points and test extension by <= 1 class column | ~10 min |
| mixed7 | weight-7 words touching a non-class column (needed only for a d = 8 claim) | the weight-7 words of RM(20,4) are line-supported too (`scratchpad/rm20_w7_classify.py`, WEIGHTED_SUPPORTS.md "Weight-7 words of RM(20,4)"); line words of weight 5/6/7 inside the non-class points, extension by <= 2 class columns via a direction table (`scratchpad/mixed_w7_check.py`) | 5 s |
| witness7 | (upper bound) | the stored planar witness, verified against G0; for seed 407 the 3+4 block scan | seconds, or hours for 407 |

Driver: `scratchpad/certify_candidate.py K SEED OUT.json` (resumable; each stage appends to the JSON).
A verdict of `d = 7 exact` makes the code a certified entry: k = 265 gives gamma = 0.960, k = 270 gives
0.959, against the current record 0.982.

## Queue (lenore, `~/repro_aanha/runs/certify/`)

| # | candidate | why | status | result |
|---|---|---|---|---|
| 1 | k = 270 seed 407 | the only draw with plane distance 8; d = 7 (gamma 0.959) or d = 8 (gamma 0.898) | **d >= 7 certified** 2026-10-29 (block6 73 min, full5 10.0 h, mixed6 3.5 min: 1 197 504 line-word candidates, 0 extensions); **mixed7 done** 2026-10-30 (`scratchpad/mixed_w7_check.py`, 233 376 line-word candidates, 0 extensions, 5 s): every weight-7 word lies in the class block; witness scan deferred | [[1746,270,>=7]]_7, gamma <= 0.959; d = 8 iff the block 3+4 scan is empty |
| 2 | k = 270 seed 19 | first k = 270 survivor, planar witness in hand | running: block6 passed (81 min), in full5 | |
| 3 | k = 265 seed 8 | first k = 265 survivor, planar witness in hand | queued after 2 | |

**Throttling (2026-10-29/30).** Lenore runs at load 31-39 with the owner's circuit-knitting jobs, so the
nice-19 jobs get 40-65% of their threads and the measured probe rate is ~1-3e5/s/thread (the idle figure
was ~1e6). Revised ETAs: 407's full weight-5 scan ~16 h after its start; each further candidate ~1 day;
the p = 29 / p = 31 certificates ~0.5 / ~1.5 days. Seed 407's witness stage (a 3+4 block scan with no
early exit if d = 8, ~1000 thread-hours) is **skipped by a watcher** (`lenore:~/repro_aanha/skip407.sh`)
so the chain proceeds to the two candidates with stored witnesses; 407 keeps the saved verdict
"d >= 7" and the scan is deferred to an idle lenore. The weight-7 census (~4000 thread-hours at this
rate) likewise needs an idle box: it stays queued but will not finish at the current throttle.

Seed 407's witness stage is a full 3+4 scan of the class block with early exit. Since 2026-10-30 an
empty scan *proves* d = 8: the weight-7 analogue of the line lemma (every word of weight <= 7 of
RM(20,4) is a line word; `scratchpad/rm20_w7_classify.py`) plus the mixed7 lookup above exclude every
weight-7 word with a non-class column for this draw (0 of 233 376 candidates extend).

## The p >= 17 ladder (m = 2 line-spread codes)

Method: `scratchpad/m2_line_sweep.py` (greedy line-spread sampler with at most L outputs per line;
exact per-line minimum weight by null-space enumeration; full-code certificate = serial MITM d <= 4 +
parallel weight-5 scan; d = 6 exact when a weight-6 line word exists and no weight <= 5 word does).
Line screens run on loranne (4 cores, 2 GB); full certificates on lenore (4 threads each) because the
MITM tables at n ~ 430-770 exceed loranne's memory and the local box reaps them under memory pressure.

| p | code | gamma | status | file |
|---|---|---|---|---|
| 19 | [[296,65,6]]_19 | 0.846 | certified (5/6 draws at k = 65; all 6 at k = 68 fail) | `runs/m2_line_sweep/p19_k65.jsonl` |
| 23 | [[424,105,6]]_23 | 0.779 | certified (15 residual checks = accidental-word limit) | `runs/m2_line_sweep/local/p23_k105_certified.jsonl` |
| 29 | [[666,175,6]]_29 ? | 0.746 if d = 6 | certificate running on lenore (4 threads) | `lenore:~/repro_aanha/runs_m2/p29_k175_certified.jsonl` |
| 31 | [[766,195,6]]_31 ? | 0.763 if d = 6 | certificate running on lenore (4 threads) | `lenore:~/repro_aanha/runs_m2/p31_k195_certified.jsonl` |
| 37 | | ~0.745 vs RS 0.800 | not started: marginal gain, 20 h certificate | |

The route stops paying at p ~ 37, where the RS family's gamma falls below what d = 6 at the
accidental-word limit can give.

## Follow-up work (not started)

- Weight-7 census on the class block of [[1756,260,7]]_7 (cost claim): kernel `count_weight7_per_output`
  (3+4, 35-fold multiplicity asserted, tested against null-space enumeration) written 2026-10-29;
  `scratchpad/census7_1756.py` is QUEUED on lenore behind the certification chain (8 threads, ~50 h;
  output `runs/certify/census7_1756.json`). It counts block words only, which since 2026-10-30 is
  all of them (next item).
- ~~Weight-7 mixed-word lemma~~ done 2026-10-30 (lines only; `scratchpad/mixed_w7_check.py`, seed 407:
  0 mixed words of weight <= 7; same for the [[1756,260,7]]_7 draw,
  `runs/weighted/p7_m4_sphere_1756_260_mixed_w7_check.json`). The queued weight-7 census of the
  [[1756,260,7]]_7 class block therefore counts ALL weight-7 words of that code, not a lower bound.
- Catalog import + performance records for every certified candidate (exact per-output enumeration
  is 7^{dim A} with dim A ~ 56-61: out of reach; use the leading-order surrogate from the census).
