# Exact distillation performance: leading order and threshold

*2026-10-20. Engine: `qmsd/performance.py`. Store: `qmsd/data/performance/<artifact_id>.json`
(schema `qmsd.performance.v1`). CLI: `python -m qmsd perf --all`.*

This applies the single-round analysis of arXiv:2610.08012 (Prakash, *Cyclic triorthogonal
codes in prime dimension*, §8) to every code in our catalog.

## Model and definitions

- **Noise.** Each input is (1−ε)|M⟩⟨M| + ε/(p−1)·Σ_{j≠0} Z^j|M⟩⟨M|Z^{−j}. This is twirled
  depolarizing noise. **ε is the paper's convention.** The repo's δ is δ = p/(p−1)·ε.
- **Acceptance and failure.** Post-selection accepts Z(j) iff j ∈ A⊥, where
  A = rowspan(G0). Output i is wrong iff in addition j·h_i ≠ 0, where h_i is the canonical
  logical row (δ_ij on the punctures).
- **Two output-error notions** (identical for k = 1):
  - **per output**: (1/k) Σ_i P(output i wrong | accept). This is what a recursion consumes.
  - **block**: P(any output wrong | accept). This is the paper's formula.
- **Leading order.** ε_out ≈ c·ε^d with c = L_d/(p−1)^d, where L_d counts the weight-d
  accepted logical errors.
- **Threshold.** ε\* is the smallest root of ε_out(ε) = ε in (0, (p−1)/p).
  ε = (p−1)/p, the fully depolarized input, is always a trivial fixed point.

**Exactness.** Everything is computed from primal weight enumerators, via
P(j ∈ C⊥) = |C|⁻¹ Σ_x C_x (1−η)^x with η = pε/(p−1). The logical part is combined into one
integer polynomial before it is evaluated, which avoids cancellation.

## Validation

- **Paper Table 2.** All 11 codes we can build are reproduced to every printed digit: c, d, ε\*
  and P_acc(0.01) (`tests/test_performance.py`). The cyclic codes come from the spectral
  supports via `qmsd/cyclic.py`; [[19,1,3]]₃ is also built from Example 14's explicit H.
- **Catalog A_d.** The block L_d reproduces the catalog A_d for every exactly computed
  Table-3 code.
- **Monte Carlo.** `qmsd/mc_distill.py` agrees with the exact P_acc and ε_out (both kinds).
- **C kernel.** `dist_weightdist/core_c/dwq_coset` (any q ≤ 127, A plus cosets A + h_i) is
  bit-identical to numpy on p = 3, 5, 7 codes, and passes a brute-force self-check for q = 2…13.

## Results (ε convention)

| code | per-output ε_out | ε\* per output | ε\* block | how |
|---|---|---|---|---|
| [[80,1,5]]₃ | 65 ε⁵ | 0.1622 | 0.1622 | exact |
| [[79,2,4]]₃ | 8.13 ε⁴ | 0.1512 | 0.1463 | exact |
| [[72,9,3]]₃ | 31.0 ε³ | 0.0622 | 0.0539 | exact |
| [[206,37,4]]₃ paper (A_d=880) | 17.6 ε⁴ | 0.06160 | — | exact |
| [[206,37,4]]₃ ours (A_d=572) | 15.2 ε⁴ | 0.06160 | — | exact |
| [[200,43,3]]₃ | 127 ε³ | 0.0284 | — | exact |
| [[667,62,4]]₃ | 132 ε⁴ | 0.0210 | — | exact (lenore C kernel, 3 s) |
| [[215,28,5]]₃ | 6.81 ε⁵ | 0.1039 | — | exact (lenore C kernel, ~25 min) |
| [[24,1,3]]₅ | 1.5 ε³ | 0.2633 | 0.2633 | exact |
| [[20,5,2]]₅ | 38.3 ε² | 0.0230 | 0.0189 | exact |
| [[112,13,3]]₅ paper (A_d=512) | 3.43 ε³ | 0.07738 | — | exact |
| [[112,13,3]]₅ ours (A_d=396) | 1.82 ε³ | 0.07765 | — | exact |
| **[[519,106,5]]₅ ours (A_d=1904), γ<1** | **1.49 ε⁵** | **0.04314** | — | exact on 16 of 106 outputs (spread 1e-9) |
| [[44,5,3]]₇ (10 codes, A_d=600/606) | 1.11–1.14 ε³ | 0.1732 | 0.1672 | exact |
| [[237,52,6]]₁₇, γ<1 | — | **≥ 0.0462** (subcode dim 5 of 14) | — | rigorous bound (Lemma 17) |
| [[519,106,5]]₅ paper (A_d=2180) | — | ≥ 0.0173 | — | rigorous bound |

## Findings

1. **A_d/k under-counts the per-output error, and the gap grows with k.** One low-weight
   logical error flips many outputs at once, because the logical rows h_i are dense. The exact
   per-output coefficient is 2.6–5× the A_d/k value on the small codes. On the flagship
   [[519,106,5]]₅ it is **~85×**: a typical weight-5 error corrupts ~85 of the 106 outputs.
   Every δ_out and cost figure built on `distillation.delta_out_avg` (eq 39) is optimistic
   per output. The STATUS.md "δ_out≈8e-18" example is one such figure. γ, d and A_d are
   unaffected.
2. **A_d optimization does not move the threshold.** Our A_d-optimized codes lower the
   leading coefficient (per output: 1.16× for [[206,37,4]]₃, 1.89× for [[112,13,3]]₅), but ε\*
   agrees to 3–4 digits. ε\* is set by the high-weight part of the weight distribution. The
   "1.54× / 8.6× / 8500×" compounding in QUTRIT_PARETO.md and NEW_CODES.md used the total
   (block) count. The exact per-output recursion from δ_in = 0.01 gives 1.12× / 1.85× / 13×
   (rounds 1/2/3), with absolute errors ~14× above eq 39; see the corrections there.
3. **The p = 17 "δ_th = 2.74 > 1" claim is refuted.** That figure extrapolated the leading
   order (RESULTS.md §6, NEXT_STEPS.md and NEW_CODES.md are corrected;
   `p17_d6_code.json` still carries the old note). Exact thresholds are finite. The γ<1
   flagship has ε\* = 0.043.
4. **Per-output thresholds are identical across outputs** to ~1e-9 on the flagship, so
   sampling a few outputs gives the exact per-output threshold. The leading coefficients
   vary by about ±3%.
5. **The γ<1 codes have low thresholds (ε\* ≈ 0.03–0.05),** and so do the high-rate qutrit
   codes ([[200,43,3]]₃: 0.028; [[667,62,4]]₃: 0.021). By comparison, [[215,28,5]]₃ reaches
   0.104, the small k = 1 codes 0.16–0.26, and the paper's cyclic codes up to 0.45. The
   low-threshold codes need pre-distillation by a small high-threshold code before they pay off.

## Engines / tiers (`performance.assign_tier`)

| tier | when | engine |
|---|---|---|
| 1 | p^{dim A}·k ≤ 3e8 | local numpy (`weightdist.coset_weight_enumerators`) |
| 2 | ≤ 1e13 | `dwq_coset` C/OpenMP kernel (`dist_weightdist/coset_run.py`), lenore |
| 3 | larger | rigorous Lemma-17 subcode lower bound (`threshold_lower_bound_subcode`) |
| 4 | analytic Manhattan codes | none (no stored puncture set) |

The block error needs p^{dim A + k} words, so it is feasible only for small k.

## Not yet done

- [[293,68,5]]₁₉ (19¹⁰ ≈ 6e12 words per output; ~20 h on lenore even with 4 outputs).
  [[230,13,6]]₃, [[690,39,5]]₃ and [[2093,308,5]]₇ are beyond exact reach; they get bounds
  only.
- Exact per-output leading order for the large-p codes (needs the explicit min-weight codewords).
- A multi-round cost model that uses the per-output curves and accounts for correlated output
  errors.
