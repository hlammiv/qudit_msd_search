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
**First data (flat_spread sampler, 10 trials each, d certified):**

| k | unit F₇³ | weighted μ = x₁(x₁−1) |
|---|---|---|
| 4 | [[339,4,5]], γ = 2.76 | [[241,4,5]], γ = 2.55 |
| 10 | [[333,10,5]], γ = 2.18 | **[[235,10,5]], γ = 1.96** |
| 18 / 21 | [[325,18,5]] / [[322,21,5]], γ = 1.80 / 1.70 | sampler finds no cap |

Same distance, ~100 fewer qudits. It stalls at k ≥ 18 because the cube-charge points (μ ∈ {1, 6}) lie on only
2 parallel planes (98 points), and caps inside two planes are small.
**Next:** choose charge polynomials whose cube-valued points are spread across F₇³, e.g. irreducible quadrics,
or allow non-cube outputs and track their Clifford class. Then derive a distance law for weighted supports
(flats ∩ X), and try m=4, where the degree budget is larger.

## I. ⏸ Multi-output doubling (arXiv:2610.29457 open problem) — deferred by request 2026-10-22.

---
## Settled this session (context for the above)
- Distillation/threshold for `[[237,52,6]]₁₇`: [the "δ_th>1" leading-order argument is REFUTED 2026-10-20 —
  exact thresholds are finite; rigorous ε\*≥0.0462 per output, see `PERFORMANCE.md`]. Beats `[[14,3,4]]₁₇` (same C≈4.6, higher d ⇒ γ<1 vs >1)
  and the paper's analytic `[[77540,5981,15]]` (smaller + lower γ). Only real size cost is circuit-level
  (noisy-Clifford), outside the standard model.
