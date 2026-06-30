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
- ⏸ **EXACT A₆ for `[[237,52,6]]₁₇` — needs a REBALANCED weight-counter (added 2026-06-29).** The blocker:
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

## F. Complete the closure
A quick **p=7 m=4** probe — the one feasible case we haven't explicitly run — to make the p≤13 no-go airtight.

---
## Settled this session (context for the above)
- Distillation/threshold for `[[237,52,6]]₁₇`: protocol threshold is generous (δ_th>1 — suppresses any input
  error — unless A₆>2.3e7; RM-sparsity ⇒ A₆~10³–10⁴). Beats `[[14,3,4]]₁₇` (same C≈4.6, higher d ⇒ γ<1 vs >1)
  and the paper's analytic `[[77540,5981,15]]` (smaller + lower γ). Only real size cost is circuit-level
  (noisy-Clifford), outside the standard model.
