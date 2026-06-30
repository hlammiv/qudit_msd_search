# Longer-term todos (qudit-MSD code discovery)

Last updated 2026-06-29. Status tags: ✅ done · 🔄 running · ⏸ deferred · 🔭 research.

## A. Push the new γ<1 codes lower (line-spread → higher d)
The line-spread that pinned p=17 d=6 has a next rung: lower `max_ℓ|S∩ℓ|` ⇒ higher d.
- **p=17:** `[[237,52,6]]` γ=0.847 → **d≥7 ⇒ γ≈0.78** (needs `max_line≤5`; avg ≈3, plausible).
- **p=19:** `[[293,68,≥5]]` γ=0.908 → **d≥6 ⇒ γ≈0.82** (d_RM=13 ⇒ `max_line≤7`); also pin p=19's exact d.
- Cost: the d+1 certification needs the `d_max=6` MITM (the ~9e9 wall; memory-safe via `mindist_balanced`,
  compute-heavy). At p=17 the 2D codewords stayed ≥6, so d=7 is a *plausible probe*, not a sure thing.

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
