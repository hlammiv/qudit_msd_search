# Longer-term todos (qudit-MSD code discovery)

Last updated 2026-06-29. Status tags: ✅ done · 🔄 running · ⏸ deferred · 🔭 research.

## A. Push the new γ<1 codes lower — ✅ TESTED + REFUTED (2026-06-30)
The line-spread d+1 rung was tested for both flagships and **refuted by the 2D cap**: a lower-max_line set
reaches the higher *line* bound (p=17 `d_lines=7` at max_line=5; p=19 `d_lines=7` at max_line=6), but a 2D
full-span codeword punctures to the *current* distance (p=17 weight-6 ⇒ d=6, via `mindist_balanced(d_max=6)`,
57 min; p=19 weight-5 ⇒ d=5, via `min_dependent_columns(d_max=5)`, 27 s). So **`[[237,52,6]]₁₇` γ=0.847 and
`[[293,68,5]]₁₉` γ=0.9076 are at their 2D ceilings — optimal at their (n,k)**, and p=19's d is now pinned exact.
Caveat: one candidate set tested each; the 2D cap is structural (matches p=11/13), so very likely general.

## B. Fully characterize the headline codes (A_d)
- ✅ **EXACT A₆ for `[[237,52,6]]₁₇` = 147,856 — DONE 2026-06-30** via the rebalanced counter
  `qmsd/weightcount_balanced.py` (8.97e9 stream, ~5 h lenore). [Former "threshold δ_th=2.74>1" was a
  leading-order extrapolation — REFUTED 2026-10-20; rigorous ε\*≥0.0462 per output, see `PERFORMANCE.md`.] The flagship is now fully characterized. *(history below:)* The blocker:
  `count_weight_d` splits d=6 as d1=3 ⇒ left table `C(237,3)·16³ ≈ 9e9` entries (~72 GB) ⇒ **OOM**; and
  `structured_ad(jmax=1)` gives a *wrong* 0 at p=17 m=2 (contradicts the verified line-supported d=6 witness —
  not validated there). FIX = the count-analog of `mindist_balanced`: fix the left block's leading coefficient
  (left → `16²·C(237,3) ≈ 5.6e8`, ~4.5 GB) and STREAM the 9e9 right, with care for the counting multiplicity
  (each weight-6 codeword counted `C(6,3)=20`× and per-scalar). Then A₅ for `[[293,68]]₁₉` (same wall, smaller).
  Compute ~hours on lenore. **DEFERRED** — confirmation-grade only: the threshold is already settled (see below),
  so this gives the exact δ_out curve + paper completeness, not a blocker.
- RS family A_d: ✅ already exact via the MDS closed form `A_d = C(n,d)(p−1)`.

## C. Deeper / broader A_d optimization
- 🔄 Flagship `[[519,106,5]]₅` A₅ sweep on lenore (push below 1904).
- Deeper sweeps on qutrit `[[206,37,4]]₃` (A_d=572) and `[[112,13,3]]₅` (A_d=396); apply the validated m=4
  `structured_ad` + `weightcount` to other m=4 codes.

## D. 🔭 Non-RM triorthogonal codes (the only γ<1 route for p≤13)
The flat-cap closed *punctured RM* for p≤13 (m=2 2D-capped, m=3 plane-capped). A different triorthogonal family —
algebraic-geometry/Hermitian codes, or a direct triorthogonal-matrix search — is the only way past it. Big,
high-risk; possibly its own paper.

## E. Consolidate / write-up
The session is a paper's worth: new γ<1 codes (p=17/19), the optimal RS family (p≥23), the **flat-cap no-go
boundary**, the A_d optimizations (qutrit/ququint), the distillation analysis, and the tooling. Synthesize.

## F. ✅ (2026-06-30): p=7 m=4 hyperplane-capped — high-k closed; ⚠ mid-k open (corrected by adversarial verify)
Probed p=7 m=4 (N=2401, D=326) at the γ<1 density (high k): true d = 2–3 (via `min_dependent_columns(d_max=3)`,
e.g. k=314 ⇒ d=3, k=320 ⇒ d=2), vs d>6.6 needed for γ<1. The codim-1 flats (hyperplanes, 7³ pts) over-puncture
(~k/7 each) — same mechanism as p=11/13 m=3 (plane) and p≤13 m=2 (2D). **CORRECTION (adversarial verification,
RESULTS.md):** this is rigorous only for **k ≥ 313** (Hamming/MDS on redundancy R=326−k); the **mid-k window
k ∈ [110, 312] is NOT certified** (`min_dependent_columns` overflows at R>21 rows there, and the "~k/7" heuristic
is loose — binding weight-21 2-flat words cost only ~k/114 per hyperplane). So "airtight for p=7" was an
**overstatement**: p=7 is closed only for k ≥ 313; the mid-k m=4 window is the one open loose end (and is D's best
non-RM target — see D_PLAN.md §2.3). **p=11, 13 closures are airtight** (fully validated). (p=5 is the exception:
it DOES cross at m=4 — the flagship `[[519,106,5]]` γ=0.987, A₅=1904 — just no *smaller* code.)

## G. 🔄 Protocol-cost search objective (Option C; started 2026-10-22, project for Aanha)
After arXiv:2610.29457 the qutrit standard library includes the Prakash–Saha [[6m+2,3m−2,2]]₃ codes. Against
it, our codes help only as **middle/late rounds after d=2 first rounds**: [[206,37,4]]₃ saves 1.2–1.5×,
[[215,28,5]]₃ 1.1× and [[204,39,3]]₃ 1.09×, while 7 of 10 qutrit catalog codes save nothing. So γ is the
wrong search target. **Objective:** `qmsd/objective.py` `marginal_saving(candidate, standard_library(p))`.
Candidates use `LeadingOrderCurve` from `structured_ad` (per-output L_d = Σ_w (w−d)·hist[w]/k), which matches
exact curves to ≲3% on the 10 qutrit catalog codes (`validate`). CLI: `python -m qmsd value <ids> [--surrogate]`.
**Next:** plug `marginal_saving` into the search ranking (cheap screen: d, k/n, then the surrogate on the top
few), and target d=3–4 at k/n ≈ 0.15–0.2, where [[206,37,4]] sits.

## H. 🔭 Weighted RM supports (p=7 lever; started 2026-10-22)
Charge picture (arXiv:2610.29457 §4): any codeword μ of RM(3r,m)^⊥ = RM(m(p−1)−1−3r, m) gives a weighted
triorthogonal support X = supp μ, which is shorter than F_p^m whenever 3r < m(p−1)−1. At p=7, m=3 (r=5) the budget is
degree 2: μ = x₁ (294 pts) or x₁(x₁−1) (245 pts, vs 343), full rank 55–56 (`qmsd/weighted.py`,
`tests/test_weighted.py`). Outputs must sit at cube-charge points (a third of them at p=7). Inputs at non-cube
charges consume U^μ magic states (counted as one input each, as for the paper's [[20,1,4]]₇).
**Charge census (2026-10-25, `qmsd/weighted.py: quadric, class_points`).** Up to affine equivalence
every degree-2 charge is a diagonal quadric plus a linear term or constant (39 forms at p = 7, m = 3).
The three charge classes (cubes, 2·cubes, 4·cubes) are interchangeable by rescaling g, so outputs can
sit in any one class. x₁(x₁−1) puts its 98 cube points on 2 planes (49 per plane); the parabolic
quadric x₁²+x₂²+x₃ (|X| = 294, 98 points per class) and x₁²+x₂²+x₃²+1 (|X| = 287, 84 points in the
2·cube class) spread them to ≤ 16 per plane and ≤ 4 per line.

**Cap search on them (flat_spread, 10–12 trials, d certified by MITM):**

| k | unit F₇³ (343) | x₁(x₁−1) (245) | x₁²+x₂²+x₃ (294) | x₁²+x₂²+x₃²+1 (287) |
|---|---|---|---|---|
| 10 | [[333,10,5]] γ 2.18 | [[235,10,5]] γ 1.96 | [[284,10,5]] γ 2.08 | [[277,10,5]] γ 2.06 |
| 18 | [[325,18,5]] γ 1.80 | none | [[276,18,5]] γ 1.70 | [[269,18,5]] γ 1.68 |
| 21 | [[322,21,5]] γ 1.70 | none | [[273,21,5]] γ 1.59 | [[266,21,5]] γ 1.58 |
| 24 | [[319,24,5]] γ 1.61 | — | [[270,24,5]] γ 1.50 | **[[263,24,5]] γ 1.49** |
| 25–28 | d = 3–4 | d ≤ 2 | d = 3–4 | d = 3–4 |

The spread charges remove the stall: d = 5 now reaches the same k = 24 ceiling as the unit space
(the cap bound of the sampler, not the weighting), with ~56 fewer qudits. [[263,24,5]]₇ γ 1.49 is
the best p = 7 d = 5 code of this size (paper: [[322,21,5]]-type unit codes, γ ≥ 1.70).
**Scored against the p = 7 standard library (2026-10-26, lenore; `objective.surrogate_from_matrices`,
streamed MITM, 90–170 s per code).** Library: RS [[6,1]], RM [[48,1]], cyclic [[7,1,3]] / [[28,1,5]] /
[[29,1,6]] / [[49,1,7]], doubled [[20,1,4]]₇. Saving = cost(library) / cost(library + code):

| code | L_d (per output) | c_ε | 1e-2→1e-30 | 1e-3→1e-15 | 1e-3→1e-30 | protocol it enters |
|---|---|---|---|---|---|---|
| **[[263,24,5]]₇** ellipsoid (2·cube) | 49 | 0.0063 | **2.55×** (210 → 82) | **2.02×** (28.8 → 14.3) | **2.13×** (141 → 66) | [[7,1,3]] → code; code alone; RS → code |
| [[270,24,5]]₇ parabolic | 75 | 0.0096 | 2.49× | 1.95× | 2.07× | same |
| [[266,21,5]]₇ ellipsoid | 48 | 0.0062 | 2.21× | 1.74× | 1.84× | same |
| [[319,24,5]]₇ unit | 138 | 0.0177 | 2.11× | 1.57× | 1.50× | same |
| [[235,10,5]]₇ x₁(x₁−1) | 24 | 0.0031 | 1.19× | none | none | [[7,1,3]] → code |

No code enters at 1e-2 → 1e-15 (two rounds of [[7,1,3]]₇, cost 52.6, are already cheap). This is the
first p = 7 search code that beats the cyclic + doubled library, and by ~2×, not the 1.1–1.6× of the
qutrit codes. The weighting itself is worth ~1.2× over the unit-space code of the same k and d
(2.13× vs 1.50× at 1e-3 → 1e-30), on top of the 56 fewer qudits.

Caveats: (i) the surrogate is leading-order only. In every protocol the code is applied at input
noise ≤ 1e-3 (n·ε ≤ 0.26), where on the exact [[44,5,3]]₇ curves the leading order is within 4% of
the exact ε_out and P_acc matches (1−ε)ⁿ to four digits; the margins (6e-18 vs 1e-15, 6e-31/6e-33 vs
1e-30) are ≥ 150×. An exact per-output curve needs 7^{dim G0 ≈ 32} words, out of reach; a tighter
Lemma-17 bound than the dim-5 one that left Aanha's [[322,21,5]]₇ undetermined would be the way to
make this rigorous. (ii) The leading-order "threshold" 0.72 is meaningless. (iii) The ellipsoid
code's outputs sit in the 2·cube class; rescaling g by 4 (4·2 = 1 mod 7) makes them cube-class
(ordinary magic outputs), after which the inputs at non-cube charges are U^μ states, counted as
one raw state each as for [[20,1,4]]₇.
**m = 4 (2026-10-26).** Engine first: `mindist_nb.count_weight5_per_output` (thread-parallel 2+3
census with per-output classification; exact vs the streamed MITM on three codes) makes the
surrogate affordable at n ≈ 2000. Scored the free-form unit-space **[[2093,308,5]]₇** (γ 1.19,
`p7_2093_308_5.json`): A_5 = 1560 only (260 classes) but each weight-5 error hits 86% of the
outputs, L_d = 1338, c_ε = 0.17 → saving **4.11× at 1e-2 → 1e-30** (210 → 51, via [[7,1,3]]₇ →
code at ε ≈ 1e-6) and **2.94× at 1e-3 → 1e-30**; none at 1e-3 → 1e-15 (P_acc(1e-3) = 0.12).
Data: `runs/weighted/p7_m4_unit_2093_308_5_census.json`. So the high-k m = 4 regime beats the
m = 3 weighted caps (2.5×) on cost, as the rate n/k = 6.8 vs 11 suggests.
Weighted m = 4 (r = 7, degree-2 charges; census `scratchpad/charge_census_m4.py`): the sphere
x₁²+x₂²+x₃²+x₄² has |X| = 2016 and 672 class points with ≤ 21 on a plane. But outputs confined to
one class lie on the degree-4 surface (g−1)(g−6) = 0, so the 35 polynomials (g−1)(g−6)·RM(3,4)
vanish on every class point: the restricted rank is 291, and near it the non-output class points
are barely checked (k = 290: d = 1, 56 class columns with zero residual check; k = 280: d = 4).
**[[1746,270,6]]₇, γ = 1.042 (2026-10-26)** — the best p = 7 γ so far (free-form unit space:
1.19). Certificate: no weight ≤ 4 (numba MITM, lenore and local), no weight-5 (thread-parallel
full 2+3 scan, 9919 s on lenore), weight-6 witness on columns [1657, 1669, 1716, 1584, 1624, 1727]
of G0 with coefficients (3, 3, 3, 3, 1, 1), all six on class points. Puncture set, charge and
witness: `runs/weighted/p7_m4_sphere_1746_270_cols.json`. The witness came from the **class
block**: the 402 non-output class columns see only 21 independent checks (33 nonzero rows of
G0, rank 21), so low-weight words live there; a 3+3 MITM on that block
(`mindist_nb.weight6_witness_3x3`, 510 s) found it, where the 2+4 scan on the full code would
have taken weeks (measured probe rate ~1e6/s/thread under load). Random-code heuristic for the
block ([402, 381]₇): ~0.5 weight-6 words and ~160 weight-7 words expected, so d = 6 at k = 270
is "typical"; words mixing in non-class columns face 35 extra checks (expected ~1e-9).
**Cost (census 2026-10-26):** A_6 = 24 (4 classes) on the class block, L_6 = 1.33, c_ε = 2.9e-5 →
saving **4.32× at 1e-2 → 1e-30** (210 → 48.7, [[7,1,3]]_7 cyc-B → [[1746,270,6]]_7) and **3.60× at
1e-3 → 1e-30** (141 → 39.2, [[6,1,?]]_7 RS → [[1746,270,6]]_7); 1.08× at 1e-2 → 1e-15; none at
1e-3 → 1e-15. Slightly better than [[2093,308,5]]₇ (4.11 / 2.94) with 347 fewer qudits.
Bracketing (`wm4_bracket.py`, class-block filter): k = 275 samples give d_block 4, 5; k = 265
samples give 3 and 6; k = 260: 4 — d depends on the sample as much as on k, so the search at
fixed k needs many cheap block-filtered draws (or a geometric constraint on the outputs), not one
draw. The k = 265 block-6 sample was fully certified: **[[1751,265,6]]₇, γ = 1.054** (no weight
≤ 5 in the full code, 9604 s), confirming that the block distance was the true distance there
(`runs/weighted/p7_m4_sphere_bracket_2026-10-27.json`, seed 7).
**[[1756,260,7]]₇, γ = 0.982 — the first γ < 1 search code at p = 7 (2026-10-27).** Sphere
charge, outputs = 260 random class points (`wm4_bracket.py`, seed 7, k = 260 sample 1; set and
certificate in `runs/weighted/p7_m4_sphere_1756_260_cols.json`). Three-part distance proof:
(i) no weight ≤ 5 in the full G0 (serial MITM for ≤ 4, thread-parallel full 2+3 scan for 5,
12243 s); (ii) no weight-6 word on the 412 class columns (3+3 MITM); (iii) no weight-6 word
touching a non-class column: the 35 kernel rows (g−1)(g−6)·RM(3,4) force the non-class part
u = (g−1)(g−6)·e_N into RM(3,4)^⊥ = RM(20,4) (d = 5), so u has weight 5 or 6 and the class
part at most one column; every full-support weight-5/6 word of RM(20,4) is line-supported
(verified by enumerating point sets up to affine equivalence, 76242 sets of span 2–4, none
carries a full-support word), and all 1,197,504 line words inside N fail to extend to ker G0
(`scratchpad/mixed_w6_check.py`, `runs/weighted/p7_m4_sphere_1756_260_mixed_w6_check.json`).
A weight-7 witness (3+4 on the class block, 5.8 h; G0 columns [523, 1045, 1117, 19, 490, 1260,
1556], coefficients (3, 3, 1, 1, 1, 3, 1)) makes d = 7 exact. Weight-7 words are rare: the first
block of 8 leading columns held none, so the random-model ~160 is an overestimate. The cost needs
the full weight-7 census (3+4 on the block, ~1.5e12 probes: days at the current probe rate).
**Todo (queued 2026-10-27, not started):**
- [~] Cost of [[1756,260,7]]₇ — weight-7 counting kernel in progress 2026-10-29 (see `CERTIFY_P7.md`). (a) rigorous bound now — count every accepted error of weight ≥ 7
  as a failure (ε_out ≤ Σ_{w≥7} C(n,w) ε^w / P_acc; ≤ 3e-24 behind one [[7,1,3]]₇ round) and
  run the protocol optimizer on that curve; (b) exact leading order — full weight-7 census on
  the class block (3+4, ~1.5e12 probes, ~12 h on an idle lenore) → `marginal_saving`.
- [x] Multi-seed sweep, cheap stages (plane filter + block weight ≤ 5), run on loranne
  2026-10-27/28: 1323 draws; survivor rate 13% at k = 265, 8% at k = 270, 0 at 275 (accidental
  non-planar words take over; verified). 79 + 51 candidates in `runs/weighted/plane_sweep/`
  (`summary.json`); k = 270 seed 407 has d_plane = 8. See WEIGHTED_SUPPORTS.md.
- [~] Certify the candidates on lenore — STARTED 2026-10-29 at 8 threads, nice 19 (chain 407 → k270 s19 → k265 s8; driver
  `scratchpad/certify_candidate.py`, status in `CERTIFY_P7.md`). Original plan: order k = 270 seed 407 first
  (3+3 weight-6 block scan ~10 min, full weight ≤ 5 scan ~2.75 h, then weight-7 witness; no
  witness after the full 3+4 scan would prove d = 8, γ 0.898), then other k = 270 survivors
  (γ 0.959 at d = 7), then k = 265 (γ 0.960). Stop after the first success per k.
- [ ] The other degree-2 charges at m = 4 (non-square coefficient, hyperbolic form): recompute
  the class-block rank and rerun the cheap sweep; not started.
- [ ] **Symmetric output sets (orbit sampler), queued 2026-10-30.** Draw the outputs as a union of
  orbits of a subgroup G ≤ O(4, F₇) acting on the 672 class points of the sphere, instead of a uniform
  random subset. If d ≥ 7 survives at k ≈ 265–270, ker G0 is G-invariant and every distance scan
  (weight-5 full scan, weight-6 block scan, weight-7 census / witness) shrinks by ~|G| by enumerating
  one orbit representative per word. Cheap first test (half a day, loranne-sized): enumerate small
  subgroups (|G| 20–60), build orbit unions of size 260–275, run the plane filter and the random-count
  model; only then write the orbit-aware wrappers around the existing MITM kernels (the kernels do not
  change). Prompted by the binary fold-and-lift certificates of Zenodo 23466024 (not applicable
  directly: characteristic 2 and cyclic structure). Risk: symmetry can force collinear
  configurations and planar words; the census decides.
**Next:** finish the weight-7 witness / census; bracket further (k = 255–265) using the class block
as the cheap filter before the 2.75 h full weight-5 proof; then plug `marginal_saving` into the
search as the ranking. r = 6 (degree-5 charges) was probed 2026-10-27: five-hyperplane charge
(|X| = 1110, no class block) gives d = 5 at k = 180/185 (γ 1.020 / 1.000) and collapses above —
a near miss, see WEIGHTED_SUPPORTS.md.

## I. ⏸ Multi-output doubling (arXiv:2610.29457 open problem) — deferred by request 2026-10-22.

---
## Settled this session (context for the above)
- Distillation/threshold for `[[237,52,6]]₁₇`: [the "δ_th>1" leading-order argument is REFUTED 2026-10-20 —
  exact thresholds are finite; rigorous ε\*≥0.0462 per output, see `PERFORMANCE.md`]. Beats `[[14,3,4]]₁₇` (same C≈4.6, higher d ⇒ γ<1 vs >1)
  and the paper's analytic `[[77540,5981,15]]` (smaller + lower γ). Only real size cost is circuit-level
  (noisy-Clifford), outside the standard model.
