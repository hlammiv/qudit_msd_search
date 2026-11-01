# RESULTS: consolidated final results with provenance

*2026-10-31. One table per result, each number copied from the lab-notebook document or data file
named next to it. The 2026-06 adversarial-verification ledger that previously lived at this path is now
`VERIFICATION_LEDGER_2026-06.md` (its p = 17/19 and A_d items are carried over below where still
current). Conventions: gamma = ln(n/k)/ln d; eps is the per-input error convention of arXiv:2610.08012
(the repo's delta = p/(p-1) eps); "exact" means the certificate's lower and upper bounds meet
(`CATALOG.md`). Where two source documents disagree the discrepancy is stated in a bracketed remark
(`[disagreement: ...]`) rather than resolved.*

Source documents: `COVERAGE.md`, `WEIGHTED_SUPPORTS.md`, `CERTIFY_P7.md`, `PERFORMANCE.md`,
`SMALL_CODES.md`, `CATALOG.md`, `NEXT_STEPS.md` (section H), `D_NEGATIVE_RESULT.md`. Data:
`runs/weighted/`, `runs/certify/`, `runs/coverage/`, `qmsd/data/catalog/`, `qmsd/data/performance/`.

## 1. Headline results

- **First proved gamma < 1 code at p = 11: [[13332,1309,22]]_11, gamma = 0.751** (d = 22 exact: the theorem
  plus an explicit weight-22 word). Output set = the full zero set of a cubic; distance follows from a one-line
  theorem (g*e lies in RM_11(29,4), distance 22), no scan needed (`LEVEL_SETS.md`, 2026-10-31;
  catalog `imported-p11-n13332-k1309-c7aae3253e15`). Built at p = 17: [[78455,5066,>=68]]_17,
  gamma <= 0.649 (proved); at p = 19 a quartic gives [[123284,7037,>=57]]_19, gamma <= 0.708 (rank check pending).

- **First gamma < 1 magic-state-distillation codes at p = 7**: [[1745,271,7]]_7, gamma 0.957 (four independent
  draws, all d = 7 exact, 2026-11-01); [[1746,270,7]]_7, gamma 0.959 (two
  independent draws, plane-sweep seeds 19 and 407, both d = 7 exact, 2026-10-30), [[1751,265,7]]_7
  (0.970, exact) and the first one found, [[1756,260,7]]_7 (0.982, exact, 2026-10-27). All are weighted
  sphere-charge supports on F_7^4 (`WEIGHTED_SUPPORTS.md`, `CERTIFY_P7.md`, `runs/certify/`). No
  referenced paper has a gamma < 1 code at p = 7 (`WEIGHTED_SUPPORTS.md`, "Result").
- **The m = 2 line-spread ladder**: [[237,52,6]]_17 gamma 0.847 (330x smaller than the paper's analytic
  [[77540,5981,15]]_17), [[296,65,6]]_19 0.846, [[424,105,6]]_23 0.779, [[666,175,6]]_29 0.746,
  [[766,195,>=6]]_31 <= 0.763 (0.703 if d = 7) (`COVERAGE.md`, `CERTIFY_P7.md`).
- **Reed-Solomon family p >= 23, closed form** d = r_max - k + 2, A_d = C(n,d)(p-1): min-gamma member
  from [[18,5,4]]_23 (0.924) to [[73,24,9]]_97 (0.506), n <= 73; 329 (p, k) members at primes 5..97 with
  exact thresholds (`runs/coverage/rs_family_all.json`, `SMALL_CODES.md`).
- **No small gamma < 1 code at p <= 13**: forced-grid theorem (evaluation codes) + flat cap, p = 11 misses by
  one distance unit (d = 4, need 5), p = 13 by one (`D_NEGATIVE_RESULT.md`); weighted supports do not
  reopen p = 3, 11, 13 at any certifiable size (`WEIGHTED_SUPPORTS.md`, "Do weighted supports reopen").
- **Rigorous cost**: the gamma < 1 code [[1756,260,7]]_7 keeps a *proved* saving of 4.14x over the p = 7
  standard library at 1e-2 -> 1e-30, where the d = 5, 6 codes keep only 1.27x / 1.33x of their
  surrogate 4.11x / 4.32x; rigorous one-round thresholds 3.4e-3 / 1.4e-3 / 2.5e-3
  (`PERFORMANCE.md`, "Rigorous bounds behind the surrogate"; `runs/weighted/rigorous_savings.json`).
- **Weight-<= 7 line lemma for RM_7(20,4)** (every word of weight <= 7 is supported on an affine line;
  A_5 = 17 287 200, A_6 = 11 524 800, A_7 = 18 110 400), which makes the mixed stage of every p = 7
  certificate a 5-second lookup and shows the full weight-5 scan is logically redundant
  (`WEIGHTED_SUPPORTS.md`, "Weight-7 words of RM(20,4)"; `CERTIFY_P7.md`).
- **The p = 7 window edge is a law**: the accidental-word count on the 402-column class block with 291 - k
  checks gives P(no weight-6 block word) = 0.92 at k = 270, 0.58 at 271, 0.02 at 272, 2e-12 at 273; all
  20 symmetric (orbit-union) survivors at k = 272 / 275 have weight-6 block words (a sampled census
  shows symmetric blocks carry ~30x the law's weight-6 count). The law matches exact counts at p = 11. The direction is
  closed at p = 7 for the sphere charge (`WEIGHTED_SUPPORTS.md`, "Orbit-union output sets").

## 2. Best certified code per prime, p = 5 .. 97

From `COVERAGE.md` (2026-10-28, updated through 2026-10-31). `artifact_id` is the catalog record
(`PYTHONPATH=. python3 -m qmsd catalog show <id>`); built-in records (paper codes, `new-*`) are generated
in code (`CATALOG.md`). p = 3 is included for completeness (gamma > 1, the paper's code).

| p | code | gamma | d certificate | family | artifact_id | data |
|---|---|---|---|---|---|---|
| 3 | [[206,37,4]]_3 (paper's puncture set; A_d 572 ours vs 880) | 1.239 | exact, MITM | punctured RM(3,5), cap-set A_d optimisation | `optimized-qutrit-ad572` (paper: `paper-table3-p3-m5-n206-k37`) | `qmsd/data/performance/`, `PERFORMANCE.md` |
| 5 | [[519,106,5]]_5 (paper) | 0.987 | exact, MITM (10 s) | punctured RM(5,4) | `paper-table3-p5-m4-n519-k106` (A_d 2180; ours `optimized-p5-flagship-ad1904`, A_d 1904) | `PERFORMANCE.md` |
| 7 | [[1746,270,7]]_7 (seed 19; also seed 407) | 0.959 | exact: block6 + mixed6 (+ full5 cross-check), planar weight-7 witness | weighted RM(7,4), sphere charge | `imported-p7-n1746-k270-c93a88f0735a` (seed 19), `imported-p7-n1746-k270-c833577e4463` (seed 407) | `runs/certify/p7_k270_seed19.json`, `p7_k270_seed407.json`, `p7_certified_2026-10-30_bundle.json` |
| 11 | [[100,21,4]]_11 | 1.126 | exact, MITM | punctured RM(6,2), m = 2 flat-spread search (loranne) | `imported-p11-n100-k21-20afd92fd535` | `runs/coverage/p11_m2_k21.json` |
| 13 | [[145,24,5]]_13 (600 trials at k = 24) | 1.118 | exact, MITM | punctured RM(7,2), m = 2 search | `imported-p13-n145-k24-74cc1dfb13cd` | `runs/coverage/small/p13_m2_k24_t600.json` |
| 17 | [[237,52,6]]_17 | 0.847 | exact, MITM two-way; A_6 = 147 856 | punctured RM(10,2) | `new-p17-d6` | `PERFORMANCE.md`, `NEXT_STEPS.md` B |
| 19 | [[296,65,6]]_19 (was [[293,68,5]]_19, 0.908) | 0.846 | no weight <= 5 (MITM) + explicit weight-6 line word | punctured RM, m = 2 line-spread sweep | `imported-p19-n296-k65-211e1fbb4b1a` (older: `new-p19-d5`) | `runs/m2_line_sweep/p19_k65.jsonl` |
| 23 | [[424,105,6]]_23 (also [[429,100,6]]_23, 0.813); RS [[18,5,4]]_23 0.924 | 0.779 | no weight <= 5 (MITM, 1 h) + explicit weight-6 line word | m = 2 line-spread sweep | `imported-p23-n424-k105-9d5204d5ccf4` | `runs/m2_line_sweep/local/p23_k105_certified.jsonl` |
| 29 | [[666,175,6]]_29; RS [[22,7,4]]_29 0.826 | 0.746 | no weight <= 5 (MITM, 86 142 s at 4 throttled threads) + explicit weight-6 line word | m = 2 line-spread sweep | `imported-p29-n666-k175-4f8ad426d19b` | `runs/certify/p29_k175_certified.jsonl`, `p29_k175_bundle.json` |
| 31 | **[[766,195,>=6]]_31** (PARTIAL: d = 6 or 7); RS [[25,6,5]]_31 0.887 | **<= 0.763 (0.703 if d = 7)** | lower bound 6: no weight <= 5 (MITM, 27 920 s, indexed kernels); upper bound 7: no line word below weight 7, explicit weight-7 line word; weight-6 exclusion (~4e14 probes) out of reach | m = 2 line-spread sweep | `imported-p31-n766-k195-8cc6fbf9e4da` (status `partial`) | `runs/certify/p31_k195_certified.jsonl`, `p31_k195_bundle.json` |
| 37 | [[29,8,5]]_37 | 0.800 | exact, MITM; eps* = 0.383 | RS (m = 1) | `imported-p37-n29-k8-98fc2b500b10` | `runs/coverage/rs_family_p23_97.json` |
| 41 | [[31,10,5]]_41 | 0.703 | exact, MITM; eps* = 0.371 | RS | `imported-p41-n31-k10-3defb2a31a1f` | same |
| 43 | [[34,9,6]]_43 | 0.742 | closed-form MDS (MITM tables scale as (p-1)^3: infeasible from here) | RS | `imported-p43-n34-k9-3ca915a0fe45` | same |
| 47 | [[36,11,6]]_47 | 0.662 | closed-form MDS | RS | `imported-p47-n36-k11-20b8354e8790` | same |
| 53 | [[40,13,6]]_53 | 0.627 | closed-form MDS | RS | `imported-p53-n40-k13-8a65b6c9152e` | same |
| 59 | [[45,14,7]]_59 | 0.600 | closed-form MDS | RS | `imported-p59-n45-k14-ef0f30cfe20a` | same |
| 61 | [[47,14,7]]_61 | 0.622 | closed-form MDS | RS | `imported-p61-n47-k14-bc9d88fce81f` | same |
| 67 | [[51,16,7]]_67 | 0.596 | closed-form MDS | RS | `imported-p67-n51-k16-523bb7ac9d24` | same |
| 71 | [[53,18,7]]_71 | 0.555 | closed-form MDS | RS | `imported-p71-n53-k18-4cfa4c492a64` | same |
| 73 | [[56,17,8]]_73 | 0.573 | closed-form MDS | RS | `imported-p73-n56-k17-2bf6c612ca56` | same |
| 79 | [[60,19,8]]_79 | 0.553 | closed-form MDS | RS | `imported-p79-n60-k19-aca65b2fc440` | same |
| 83 | [[62,21,8]]_83 | 0.521 | closed-form MDS | RS | `imported-p83-n62-k21-fedef75ce505` | same |
| 89 | [[67,22,9]]_89 | 0.507 | closed-form MDS | RS | `imported-p89-n67-k22-048f4920965d` | same |
| 97 | [[73,24,9]]_97 | 0.506 | closed-form MDS | RS | `imported-p97-n73-k24-6ec269ebdce9` | same |

Remarks (all from `COVERAGE.md` unless stated).

- p = 5 is the only p <= 13 crossing; it is Singleton-infeasible at m <= 3 (`D_NEGATIVE_RESULT.md`).
- p = 11 equals the 2D-cap prediction of the no-go. p = 13: the prediction is [[143,26,5]]_13 (gamma 1.059);
  600 trials at k = 25-27 and 3000 flat-spread trials at k = 26 all give d = 4, so it is unrealised.
  (D_NEGATIVE_RESULT.md section 3 lists [[143,26,5]] with d = 5 as the *2D-cap prediction*, needed 6; the
  realised distance at k = 26 is 4 in every trial.)
- p = 19: 6 of 6 draws at k = 62 and 5 of 6 at k = 65 (<= 7 outputs per line) give d = 6; all 6 at k = 68
  fail, so k = 65 is the limit. p = 23: all 180 line-spread draws at k = 90-110 have line distance 6;
  k = 105 leaves 15 residual checks (the accidental-word limit). p = 29: k = 175 leaves 15 checks.
- p = 37: not started (~0.745 vs RS 0.800, 20 h certificate); the ladder stops paying at p ~ 37 (`CERTIFY_P7.md`).
- Acceptance probabilities (Monte Carlo, 2-4 x 10^5 shots, within 1% of (1-eps)^n; `runs/certify/pacc_mc.json`,
  `scratchpad/pacc_mc.py`): P_acc(1e-3) = 0.79 ([[237,52,6]]_17), 0.74 ([[296,65,6]]_19), 0.66 ([[424,105,6]]_23),
  0.17 ([[1756,260,7]]_7), 0.18 ([[1746,270,7]]_7); at 1e-2 the p = 7 codes accept < 1e-5, so they are
  second- or third-round codes fed at eps ~ 1e-4.

## 3. The p = 7 weighted-support family

Construction (`WEIGHTED_SUPPORTS.md`, "Construction"; `qmsd/weighted.py`): charge g = x1^2+x2^2+x3^2+x4^2 on
F_7^4, |X| = 2016, rank 326 of RM(7,4) on X; outputs are k random points of the cube class g in {1, 6}
(672 points). The 35 polynomials (g-1)(g-6) RM(3,4) vanish on every class point, so the class rank is
291 (the ceiling for k) and the 672 - k non-output class points ("class block") see only 291 - k checks.
(The "Beyond the sphere" table in WEIGHTED_SUPPORTS.md lists the sphere at r = 6, where the ranks are
210 / 195 / 15; at r = 7 they are 326 / 291 / 35.)

| code | gamma | d | certificate stages (machine, wall) | census (A_d, L_d per output, c_eps) | artifact_id / data |
|---|---|---|---|---|---|
| [[1746,270,7]]_7 seed 19 | 0.959 | 7 exact | block6 4899 s (82 min), full5 35 187 s (9.8 h), mixed6 186 s (0 of 1 197 504 line words extend), planar weight-7 witness G0 cols [141,473,717,857,921,1125,1354] coeffs (1,1,1,5,1,1,1); mixed7: no mixed word of weight <= 7 (6 hits, min weight 9). lenore, 8 threads nice 19, 2026-10-29/30 | weight-7 census RUNNING (`scratchpad/census7_draw.py`, lenore; ~30 weight-7 classes found so far vs 27 predicted) | `imported-p7-n1746-k270-c93a88f0735a`; `runs/certify/p7_k270_seed19.json`, `p7_k270_seed19_mixed7.json` |
| [[1746,270,7]]_7 seed 407 | 0.959 | 7 exact | block6 4375 s (73 min), full5 36 049 s (10.0 h), mixed6 209 s, mixed7 0 of 233 376 candidates extend, weight-7 word from the 3+4 block scan (820 s, 14 min on the indexed kernels) cols [177,337,492,14,20,651,1160] coeffs (1,1,3,3,3,4,1). d_plane = 8 (the only such draw; d = 8, gamma 0.898, is closed) | none (same n, k as [[1746,270,6]]_7, one more unit of distance) | `imported-p7-n1746-k270-c833577e4463`; `runs/certify/p7_k270_seed407.json`, `_mixed7.json`, `_bundle.json` |
| [[1751,265,7]]_7 seed 8 | 0.970 | 7 exact | block6 5275 s (88 min), full5 skipped (implied by block6 + mixed6), mixed6 235 s, planar weight-7 witness cols [317,373,438,607,667,735,1630] coeffs (1,5,1,1,5,1,1); 2026-10-30 | none | `imported-p7-n1751-k265-0ecb29c584c7`; `runs/certify/p7_k265_seed8.json` |
| [[1756,260,7]]_7 | 0.982 | 7 exact | (i) no weight <= 5 in the full code, 2+3 scan 12 243 s on lenore; (ii) no weight-6 word on the 412 class columns (3+3); (iii) mixed6: 76 242 point sets of span 2-4 enumerated, 1 197 504 line words tested, 428 s; witness 3+4 scan 5.8 h, cols [523,1045,1117,19,490,1260,1556] coeffs (3,3,1,1,1,3,1); mixed7 0 of 233 376 | A_7 = 24 (4 scalar classes, 140 raw hits / 35), L_7 = 1.29, c_eps = 4.6e-6; complete; 119 533 s (33 h) on lenore's indexed 3+4 kernel, 2026-10-31 | `imported-p7-n1756-k260-9d263f2240e4`; `runs/weighted/p7_m4_sphere_1756_260_cols.json`, `..._7_census.json`, `..._mixed_w6_check.json`, `..._mixed_w7_check.json`, `runs/certify/census7_1756.json` |
| [[1746,270,6]]_7 | 1.042 | 6 exact | no weight <= 4 (numba MITM), no weight 5 (full 2+3 scan, 9919 s), weight-6 witness cols [1657,1669,1716,1584,1624,1727] coeffs (3,3,3,3,1,1) from the class block (3+3, 510 s); mixed_w7_check (2026-10-31): no mixed word of weight <= 7, so A_6 is complete | A_6 = 24 (4 classes), L_6 = 1.33, c_eps = 2.9e-5; 1806 s | `imported-p7-n1746-k270-271d5366cc68`; `runs/weighted/p7_m4_sphere_1746_270_cols.json`, `..._6_census.json`, `..._6_bundle.json` |
| [[2093,308,5]]_7 (unit space, free-form, m = 4) | 1.19 | 5 exact | MITM | A_5 = 1560 (260 classes), each hitting 86% of the outputs, L_d = 1338, c_eps = 0.17; 11 030 s (3.1 h) | `new-p7-m4-n2093-k308`; `runs/weighted/p7_m4_unit_2093_308_5_census.json` |
| [[263,24,5]]_7 (ellipsoid charge x1^2+x2^2+x3^2+1, m = 3, outputs in the 2-cube class) | 1.49 | 5 exact | MITM (flat_spread cap search, 10-12 trials) | L_d = 49, c_eps = 0.0063 (surrogate) | not in the catalog; `runs/weighted/p7_m3_scores_2026-10-26.json`, `NEXT_STEPS.md` H |

Also certified: [[1751,265,6]]_7, gamma 1.054 (bracket seed 7, k = 265 sample 1; no weight <= 5 in the full
code, 9604 s; `runs/weighted/p7_m4_sphere_bracket_2026-10-27.json`). (The "Result" table of WEIGHTED_SUPPORTS.md, the 0.960 slip for k = 265 in CERTIFY_P7.md /
WEIGHTED_SUPPORTS.md, and COVERAGE.md's "seed 19 only" remark were corrected on 2026-10-31 to match the
certificate files: gamma 0.970 at k = 265, full5 run for seeds 19 and 407.)

**The three-part certificate** (`WEIGHTED_SUPPORTS.md`, "class-block lemma"; `CERTIFY_P7.md`). A word of
weight <= 6 either (a) lies entirely in the class block, excluded by a 3+3 meet-in-the-middle on the
<= 412 class columns (`mindist_nb.weight6_witness_3x3`, stage block6), or (b) touches a non-class
column, in which case the 35 kernel rows force its non-class part u = (g-1)(g-6) e_N into
RM(3,4)^perp = RM_7(20,4) (minimum distance 5), so u is a weight-5/6 word of RM(20,4) with at most one
class column; every such word is line-supported (lemma), and all 1 197 504 line words inside the
non-class points are tested for extension to ker G0 (stage mixed6, `scratchpad/mixed_w6_check.py`).
Hence block6 + mixed6 prove d >= 7; the thread-parallel full 2+3 weight-5 scan (stage full5, ~10 h) is an
independent cross-check that is logically redundant (noted 2026-10-30). An explicit weight-7 word,
verified against G0 (rank 6, every 6-subset of full rank), makes d = 7 exact. The weight-7 analogue
(stage mixed7, 5 s, 233 376 candidates via a 484 009-key direction table) is what a d = 8 claim needs.

**Sweep statistics** (`scratchpad/plane_sweep.py`, `runs/weighted/plane_sweep/summary.json`; a survivor has
plane distance d_plane >= 7 and no class-block word of weight <= 5): k = 250: 4/15 (27%); 255: 4/22 (18%);
260: 4/23 (17%); 265: 79/600 (13%); 270: 51/600 (8%); 275: 0/63 (0%). 129 of the 130 survivors carry an
explicit planar weight-7 word (`flat_witness_k*.jsonl`); d_plane was exact on 8 of 9 draws with known d
(`runs/weighted/p7_m4_sphere_plane_distance_check.json`).

**Window-edge law** (`WEIGHTED_SUPPORTS.md`, "Orbit-union output sets"; model
E[projective weight-w words] = C(402, w)(p-1)^(w-1)/p^(291-k)):

| k | checks | E[weight-5] | E[weight-6] | P(no weight-6 block word) | E[weight-7] |
|---|---|---|---|---|---|
| 269 | 22 | 3e-5 | 0.011 | 0.99 | 3.8 |
| 270 | 21 | 2e-4 | 0.079 | 0.92 | 27 |
| 271 | 20 | 1.4e-3 | 0.55 | 0.58 | 187 |
| 272 | 19 | 1e-2 | 3.9 | 0.02 | 1300 |
| 273 | 18 | 0.07 | 27 | 2e-12 | 9200 |
| 275 | 16 | 3.3 | 1300 | 0 | 4.5e5 |

The window closes at k = 271-272: at k <= 270 a draw that passes the plane filter almost always has d = 7
(both certified k = 270 survivors did), at k = 271 about half do (gamma 0.957), at k = 272 one in fifty
(gamma 0.955), from k = 273 none. Shift-symmetric orbit-union sets improve the *planar* population (30%
survivors at k = 272, 8 of 60 with d_plane 8, against 8% random at k = 270) but all 18 k = 272 survivors
and both k = 275 survivors have a non-planar weight-6 class-block word (`certify_candidate.py --file`,
6 min each on lenore); groups containing -I are rank-capped at k = 260. The direction is therefore
closed at p = 7 for the sphere charge (`runs/weighted/orbit_sweep/`, `scratchpad/orbit_sweep.py`,
`scratchpad/orbit_rank_probe.py`, `scratchpad/window_law.py`). The confirming plane sweep at k = 271
(`runs/weighted/plane_sweep/plane_sweep_k271.jsonl`) holds 40 draws with d_plane 6/5/4/3 on 26/11/2/1 of
them and none with d_plane >= 7, so the 58% figure is the model's, not a measured pass rate.

## 4. Cost: surrogate versus proved savings

Baseline: the p = 7 standard library (RS [[6,1]], RM [[48,1]], cyclic [[7,1,3]] / [[28,1,5]] / [[29,1,6]] /
[[49,1,7]], doubled [[20,1,4]]_7); cost = expected raw states per output; saving = cost(library) /
cost(library + code) (`PERFORMANCE.md`, `NEXT_STEPS.md` H). Surrogate = leading-order per-output curve from
the complete weight-d census. Proved = protocol optimiser run on the pessimistic rigorous curve
(`qmsd/rigorous_curve.py`: weight-d term exact, every heavier accepted error counted as a failure with the
Singleton count N_w <= C(n,w) p^(w-d+1)); data `runs/weighted/rigorous_savings.json`
(`scratchpad/rigorous_saving.py`).

| code | rigorous threshold | 1e-2 -> 1e-30: proved / optimistic (surrogate) | 1e-3 -> 1e-30: proved / optimistic (surrogate) | proved protocol at 1e-2 -> 1e-30 |
|---|---|---|---|---|
| [[2093,308,5]]_7 | 1.4e-3 | 1.27x / 4.11x (4.11x) | 2.94x / 2.94x (2.94x) | [[20,1,4]]_7 dbl -> code (166.2 vs 210.3) |
| [[1746,270,6]]_7 | 2.5e-3 | 1.33x / 4.32x (4.32x) | 3.09x / 3.60x (3.60x) | [[20,1,4]]_7 dbl -> code (158.1 vs 210.3) |
| [[1756,260,7]]_7 | 3.4e-3 | 4.14x / 4.14x (4.14x) | 2.96x / 3.44x (3.44x) | [[7,1,3]]_7 cyc-B -> code (50.8 vs 210.3) |

Surrogate table for the whole weighted family (`PERFORMANCE.md`, "Weighted-support caps at p = 7"):
[[263,24,5]]_7 2.55x (210 -> 82) / 2.02x (28.8 -> 14.3) / 2.13x (141 -> 66) at 1e-2 -> 1e-30 / 1e-3 -> 1e-15 /
1e-3 -> 1e-30; [[270,24,5]]_7 2.49 / 1.95 / 2.07; [[266,21,5]]_7 2.21 / 1.74 / 1.84; [[319,24,5]]_7 (unit)
2.11 / 1.57 / 1.50; [[235,10,5]]_7 1.19 / none / none; [[2093,308,5]]_7 4.11 / none / 2.94; [[1746,270,6]]_7
4.32 (210 -> 48.7) / none / 3.60 (141 -> 39.2); [[1756,260,7]]_7 4.14 (210 -> 50.8) / none / 3.44 (141 -> 40.9).
(PERFORMANCE.md's table cell for [[2093,308,5]]_7 at 1e-2 -> 1e-15 corrected from "none" to 1.03x on
2026-10-31, matching rigorous_savings.json; the sphere codes give 1.03x / 1.08x there.)

The reversal: the Singleton tail C(n,d+1) p^2 (eps/(p-1))^(d+1) at eps = 1e-6 sits 1e9 (d = 5), 6e13 (d = 6)
and 1e16 (d = 7) above the respective surrogates, but its absolute level falls with d (f_hi(1e-6) = 1e-13,
2e-27, 7e-32), so after one [[7,1,3]]_7 round only the d = 7 code provably reaches 1e-30 and keeps its full
4.14x, while the d = 5, 6 codes must enter from a [[20,1,4]]_7 round and keep 1.3x. The surrogate prefers
rate ([[1746,270,6]]_7 at 48.7 beats 50.8) because it assumes the random-model-expected absence of
weight-8 words; the proof favours distance, and the rigorous thresholds (one round provably improves
below 0.14-0.34%) are the first rigorous threshold statements for these codes (`PERFORMANCE.md`).

Exact-engine context (`PERFORMANCE.md`): the exact per-output coefficient is 2.6-5x A_d/k on small codes and
~85x on [[519,106,5]]_5 (1.49 eps^5, eps* = 0.04314); A_d optimisation lowers the leading coefficient
(1.16x for [[206,37,4]]_3: 17.6 -> 15.2 eps^4) but leaves eps* = 0.06160 unchanged; against the corrected
p = 3 baseline (Prakash-Saha d = 2 codes included) the qutrit search codes save 1.2-1.6x, not the earlier
20x; [[519,106,5]]_5 saves 6.2x at 1e-2 -> 1e-30; [[237,52,6]]_17 >= 3.5x and [[293,68,5]]_19 >= 4.2x
(rigorous brackets). `[disagreement: PERFORMANCE.md gives the Lemma-17 bound for [[237,52,6]]_17 as
eps* >= 0.0462 in the results table and >= 0.034 in the catalog-run table; [[230,13,6]]_3 >= 0.0499 vs
>= 0.045 likewise; for [[2093,308,5]]_7 the Lemma-17 bound >= 0.0018 and the rigorous-curve threshold
1.4e-3 are two different lower bounds.]`

## 5. Structural and negative results

- **No gamma < 1 at p <= 13 in the RM framework** (`D_NEGATIVE_RESULT.md`). Flat cap: at the gamma < 1 puncture
  density every codimension-1 flat over-punctures, so a 2D (m = 2), plane (m = 3) or hyperplane (m = 4) word
  caps d. p = 11: RM_11(6,2), k = 21 gives d = 4, need 5 (gamma 1.126); p = 13: need 6 at k = 26, best gamma
  1.059 predicted, 1.118 realised (section 2). Forced-grid theorem (rigorous): a column subset S of F_p^m is
  triorthogonal iff 1_S in RM(t, m), t = m(p-1) - 3r - 1; at every p <= 13 target 2t < p forces S in
  {empty, full grid}, so the only evaluation-code construction is ordinary RM puncturing. Boundary engine:
  the normal rational curve is triorthogonal iff 3(r-1) < p-1, so the triorthogonal-MDS distance ceiling is
  floor((p+1)/3) + 1 = 3, 3, 5, 5, 7, 7 at p = 5, 7, 11, 13, 17, 19; d = 6 first becomes available at p = 17.
  Non-evaluation routes (cyclic/BCH, one-point AG, Artin-Schreier, direct moment-search caps) are each
  closed (cyclic collapses to d <= 2, moment caps to d <= 3). Open: a high-dimensional general-position tail
  (estimated < 5%) and the unweighted p = 7 mid-k m = 4 RM window k in [110, 312]. Scope: unweighted
  evaluation codes only; the weighted p = 7 codes escape the forced-grid argument.
- **Weighted census closure at p = 3, 11, 13** (`WEIGHTED_SUPPORTS.md`, "Do weighted supports reopen";
  `scratchpad/weighted_census_general.py`, `scratchpad/flat_cap_check.py`, `runs/weighted/census_p_le_13/`):

  | p, m | r | best model window | needed d | exact flat words at that k | verdict |
  |---|---|---|---|---|---|
  | 3, 6 | 3 | none (rho = 78, \|X\| >= 432) | - | - | no window |
  | 5, 4 | 4 | none (rho = 70, \|X\| >= 475) | - | - | no window; unweighted [[519,106,5]]_5 stands |
  | 13, 2 | 7 | parabola x1^2 + x2, k = 19, gamma 0.95 | 8 | lines: weight 4-5 (5/5 draws) | dead |
  | 11, 3 | 8 | cubic sum x_i^3 + c, k = 138, gamma 0.93 | 9 | lines: weight 5 (2/2 draws) | dead |
  | 13, 3 | 11 | quadric x1^2+x2^2+x3^2+1, k = 221, gamma 0.96 | 9 | lines: weight 5 (1 draw) | dead |
  | 3, 7 | 4 | x7 + c, k = 196, gamma 0.957 at d = 7 | 7 | full-code MITM: d = 5, 4 (2/2 draws) | dead |

  Reason: RM(r, m) restricted to a line is the full function space only when r >= p-1; at p = 7, m = 4 we
  have r = 7 >= 6 (no line words), at p = 11, 13 the condition fails at m = 3 and holds only at m = 4
  (n > 10^4). At p = 3, 5, 11 every unit is a cube, so there is one class and no class-block structure.
- **Weight-<= 7 line lemma for RM_7(20,4)** (`WEIGHTED_SUPPORTS.md`; `scratchpad/rm20_w7_classify.py`, 7 min,
  `qmsd/pointsets.py`, `tests/test_pointsets.py`). A set of <= 7 points of F_7^4 fails to impose independent
  conditions on the cubics iff it contains 5 collinear points, and carries a full-support annihilator iff it
  is 5, 6 or 7 points of one line (frame enumeration: 163 185 + 6 435 689 + 2 869 210 sets of span 2/3/4;
  dependent 1 515 / 60 / 0; full-support 0 / 0 / 0; cross-check 6 x 10^5 random 7-sets, 0 violations). The
  lemma stops at weight 8 (the 8 points of an ellipse carry a full-support weight-8 word).
- **full5 redundancy** (`CERTIFY_P7.md`): block6 + mixed6 already prove d >= 7; the 10-hour full weight-5 scan
  is a cross-check (seeds 407 and 19 ran it, consistent; seed 8 and later run `--skip-full5`, cutting a
  certificate from ~12 h to ~1.5 h at 8 throttled threads).
- **Kernel speedup** (`CERTIFY_P7.md`, `scratchpad/bench_mitm/RESULTS.md`): indexed MITM right-scans give
  5.5-7.4x more probes per second per thread on every split at identical output.
- **Degree-5 charges at r = 6** (`WEIGHTED_SUPPORTS.md`): the five-hyperplane charge (|X| = 1110, no class
  block) gives d = 5 at k = 180 (gamma 1.020) and k = 185 (1.000, 3/3 draws) and collapses above; a near miss.
- **Earlier structural items** (git `c90af2c`, `NEXT_STEPS.md` A/B): the line-spread "d + 1" rung is refuted
  for [[237,52,6]]_17 and [[293,68,5]]_19 by 2D full-span words (they are optimal at their (n, k)); the
  qutrit m = 7 cap code [[1968,219,>=10]]_3 is refuted (true d = 1, 83 zero columns; catalog `refuted-qutrit-cap`).

## 6. The small-code program (`SMALL_CODES.md`)

- **Reed-Solomon family [[p-1-(k-1), k, r-k+2]]_p**, exact closed-form MDS performance engine (`qmsd.rs_family`,
  `runs/coverage/rs_family_all.json`, `scratchpad/rs_family_table.py`): 329 members at the 23 primes 5..97, every
  k admissible, distance independent of the punctured columns. k = 1 members [[p-1, 1, r+1]]_p are the best
  single-output codes at every prime >= 11: eps* 0.532 ([[10,1,4]]_11), 0.602 ([[16,1,6]]_17), 0.643
  ([[22,1,8]]_23, c = 1.28e-4), 0.671 ([[28,1,10]]_29), 0.708 ([[40,1,14]]_41), 0.742 ([[58,1,20]]_59),
  0.775 ([[96,1,32]]_97, c = 1.05e-36); P_acc(1e-2) falls from 0.90 (p = 11) to 0.38 (p = 97). Min-gamma
  members: [[9,2,3]]_11 1.369, [[11,2,3]]_13 1.552, [[13,4,3]]_17 1.073, [[15,4,3]]_19 1.203, then the
  gamma < 1 ladder of section 2 (eps* 0.396 at [[18,5,4]]_23 .. 0.391 at [[73,24,9]]_97). Below p = 11 the RS
  family is trivial (d = 2).
- **Where m = 2 search codes beat RS (p = 11, 13):** RS stops at d = 4, while the m = 2 front has d = 7 at
  k = 1-3 and d = 6 up to k = 10 at p = 11 ([[120,1,7]], [[118,3,6]], [[115,6,6]], [[113,8,6]], [[111,10,6]])
  and d = 7 at k = 3, d = 6 at k = 6 at p = 13 ([[166,3,7]], [[163,6,6]]); MITM-exact distances, exact curves
  not enumerable. Searches on loranne, `runs/coverage/small/`.
- **Where RS wins (p >= 17):** [[16,1,6]]_17 is 16 qudits; the m = 2 draws at k = 6-22 are 260-280 qudits with
  only d >= 6 certified. m = 2 pays only once its rate beats RS's n/k = p-1, i.e. from k of order p
  ([[237,52,6]]_17, [[296,65,6]]_19, [[424,105,6]]_23, ...).
- **p <= 7:** the cyclic / doubling / paper codes ([[7,1,3]]_7 eps* 0.446, [[44,5,3]]_7 0.173, ...).

## 7. Comparison with the literature (`COVERAGE.md`, figure `runs/coverage/gamma_vs_p.png`, `scratchpad/gamma_vs_p_plot.py`)

The figure plots the paper's asymptotic optimum gamma_0(p) (arXiv:2510.10852 Fig. 6 / Table 1; `qmsd asymptotic`
reproduces every tabulated value to 5 digits and extends it to p = 97), the paper's finite codes (Table 2
Manhattan family, Table 3 search codes), the best k = 1 codes of arXiv:2610.08012 (cyclic) and arXiv:2610.29457
(doubling) ([[37,1,5]]_3 2.24, [[16,1,4]]_5 2.00, [[7,1,3]]_7 1.77), and our best code per prime; bottom panel:
block sizes.

| p | paper's best finite code | gamma | n | ours | gamma | n |
|---|---|---|---|---|---|---|
| 3 | [[206,37,4]]_3 (Table 3) / Manhattan | 1.239 / 0.990 | 206 / 1.9e15 | same code (A_d 572 vs 880) | 1.239 | 206 |
| 5 | [[519,106,5]]_5 (Table 3) / Manhattan | 0.987 / 0.980 | 519 / 1.5e11 | same code | 0.987 | 519 |
| 7 | Manhattan | 0.990 | 9.6e10 | [[1746,270,7]]_7 | 0.959 | 1746 |
| 11 | Manhattan | 0.964 | 1.9e7 | [[100,21,4]]_11 | 1.126 | 100 |
| 13 | Manhattan | 0.975 | 6.1e7 | [[145,24,5]]_13 | 1.118 | 145 |
| 17 | Manhattan | 0.946 | 77540 | [[237,52,6]]_17 | 0.847 | 237 |
| 19 | Manhattan | 0.967 | 121470 | [[296,65,6]]_19 | 0.846 | 296 |
| 23 | [[17,6,3]]_23 (RS) | 0.948 | 17 | [[424,105,6]]_23 (m = 2); RS [[18,5,4]]_23 0.924 | 0.779 | 424 |
| 29 | (RS family cited) | | | [[666,175,6]]_29; RS 0.826 | 0.746 | 666 |
| 31 | (RS family cited) | | | [[766,195,>=6]]_31; RS 0.887 | <= 0.763 | 766 |
| 37-97 | (RS family cited, not tabulated) | | | RS min-gamma member | 0.800 -> 0.506 | 29 -> 73 |

(COVERAGE.md's comparison table was brought in line with its per-prime table on 2026-10-31.)
Other p = 7 literature (`WEIGHTED_SUPPORTS.md`): Saha-Prakash k = 1 codes gamma 2.6-3.3 and analytic
[[2331,70,9]]_7 (1.595); cyclic codes of arXiv:2610.08012 gamma 1.77-2.1; doubled [[20,1,4]]_7 2.16;
Krishna-Tillich reach gamma < 1 only from p = 23; gamma_0(7) = 0.508 is an n -> infinity limit.
Reading (COVERAGE.md): every finite code sits far above gamma_0(p) (0.98 vs 0.51 at p = 7); at p = 7, 17, 19 we
beat the paper's Manhattan gamma at block sizes 10^3-10^8 times smaller; at p = 11, 13 nobody crosses 1 below
n ~ 10^7; from p = 23 the RS family is the whole story.

## 8. Reproducibility

- **Catalog v2** (`CATALOG.md`): `qmsd/data/catalog/*.json`, schema `qmsd.catalog-record.v2`; every record
  rebuilds from its `construction` block (`punctured_rm` / `reed_solomon` / `weighted_rm` / `manhattan`) and
  claims exactly what its `certificate` says (`lower_bound`, `upper_bound`, `exact`, methods, witness, stages,
  machine, date); status derived (`confirmed` / `partial` / `candidate` / `refuted`); `performance` and
  `geometry` blocks are pointers. 245 v1 files migrated 2026-10-29; 733 files + built-in records = 765 loaded.
  (CATALOG.md's partial example updated 2026-10-31 to [[766,195,>=6]]_31.)
- **CLI** (`README.md`, `python -m qmsd --help`; run with `ulimit -v 4000000`, `OMP_NUM_THREADS=1`,
  `PYTHONPATH=.`): `search`, `catalog import|migrate|show`, `reconstruct --label`, `asymptotic --p`,
  `perf <id> | --all`, `protocol --p P --eps0 E --target T [--bounded <id>]`, `assess <id> [--all --store]`,
  `value <ids> [--surrogate]`, `best`.
- **Tests**: `tests/` (232 test functions), including `test_catalog_v2.py` (rebuild + certificate +
  witness of every record), `test_mindist.py`, `test_objective.py`, `test_pointsets.py`,
  `test_rigorous_curve.py`, `test_rs_family.py`, `test_performance.py` (paper Table 2 to every digit),
  `test_protocol.py`, `test_weighted.py`.
- **Machines**: local box (15 GB; heavy jobs capped at 4 GB, single-thread BLAS); lenore (remote, 32 cores,
  jobs at nice 19 under the owner's load: ~1-3e5 probes/s/thread throttled vs ~1e6 idle; all p = 7
  certificates, censuses, the p = 29/31 weight-5 scans); loranne (4 cores, 2 GB: plane sweeps, line screens,
  the p = 11/13 m = 2 searches). `CERTIFY_P7.md`, `COVERAGE.md`, `SMALL_CODES.md`.
- **Scripts behind each table** (`scratchpad/`): section 2: `m2_line_sweep.py`, `line_witness.py`,
  `rs_family.py`, `pacc_mc.py`; section 3: `plane_sweep.py`, `plane_sweep_summary.py`, `flat_distance4.py`,
  `flat_witness.py`, `block_test.py`, `certify_candidate.py`, `mixed_w6_check.py`, `mixed_w7_check.py`,
  `rm20_w7_classify.py`, `rm20_w7_sample_check.py`, `census7_1756.py`, `census7_draw.py`, `orbit_rank_probe.py`,
  `orbit_sweep.py`, `window_law.py`, `census_r6.py`; section 4: `rigorous_saving.py` (+ `qmsd/rigorous_curve.py`,
  `qmsd/objective.py`, `qmsd/protocol.py`); section 5: `weighted_census_general.py`, `flat_cap_check.py`,
  `p3m7_probe.py`, `verify_p2.py`, `verify_proof.py`; section 6: `rs_family_table.py`, `small_codes_build.py`;
  section 7: `gamma_vs_p_plot.py`; catalog: `catalog_enrich_v2.py`. Engines: `qmsd/mindist_nb.py`
  (`count_weight5_per_output` 2+3, `weight6_witness_3x3` / `count_weight6_per_output` 3+3,
  `weight7_witness_3x4` / `count_weight7_per_output` 3+4, `method="index"`), `qmsd/weighted.py`,
  `qmsd/pointsets.py`, `qmsd/performance.py`, `qmsd/rs_family.py`.

## 9. Open items (`NEXT_STEPS.md` H, `CERTIFY_P7.md`, `COVERAGE.md`, `PERFORMANCE.md`)

- **Weight-7 census of [[1746,270,7]]_7 seed 19** (`scratchpad/census7_draw.py`, lenore): running; ~30 weight-7
  classes found so far against 27 predicted; makes the seed-19 cost exact (at leading order the k = 270 d = 7
  codes can only be cheaper than [[1746,270,6]]_7 unless L_7 exceeded ~1e5). Seed 407's census (~33 h) not started.
- **Aanha's Monte Carlo tasks** (`docs/aanha_mc_guide/main.tex`, `docs/aanha_mc_guide2/main.tex`, figures and
  `numbers.json` from `make_figures.py`): guide I: the sampler `qmsd/mc_distill.py`, a point and a curve, the two
  numbers the paper quotes, threshold by bisection; guide II: Task A thresholds of the large codes by sampling
  (bisection on a cheap code first), Task B within-block independence and the two-round experiment, Task C
  biased noise. Deliverables not yet returned.
- **Exact d of [[766,195,>=6]]_31**: needs a weight-6 exclusion of ~4e14 probes (2+4), out of reach; gamma is
  0.763 or 0.703.
- **p = 13 [[143,26,5]]_13 prediction** (gamma 1.059): unrealised; 600 + 3000 trials at k = 26 give d = 4 and the
  arc-climb / cap-set / near-cap samplers stall at k = 15 at p = 13.
- **Weight-8 words of RM_7(20,4)** (theory; needed only for d = 9 certificates, which also need a 4+4 kernel);
  the other degree-2 charges at m = 4; a charge with a vanishing-stabiliser block at r = 6.
- **Cost model gaps** (`PERFORMANCE.md`): exact per-output curves for the large codes are out of reach (7^56-7^66,
  19^10 words); costs count magic states only (no Clifford circuit cost, footprint or latency); the weighted
  codes consume three species of U^mu magic state, counted as one raw state each (`WEIGHTED_SUPPORTS.md`,
  "Caveats"); under a single-species accounting a conversion cost would be needed.
- **The unweighted p = 7 mid-k m = 4 RM window and the < 5% high-dimensional tail** of the p <= 13 no-go
  (`D_NEGATIVE_RESULT.md` section 7).
- **Manuscript**: not started (`NEXT_STEPS.md` E: "Consolidate / write-up").
