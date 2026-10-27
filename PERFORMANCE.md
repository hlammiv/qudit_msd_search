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
| [[519,106,5]]₅ paper (A_d=2180) | 1.62 ε⁵ | 0.04314 | — | exact on 16 of 106 outputs (B₅ = 2180 reproduced) |
| [[230,13,6]]₃ | — | ≥ 0.0499 | — | rigorous bound (subcode dim 13 of 38) |
| [[690,39,5]]₃ | — | ≥ 0.0134 | — | rigorous bound (dim 12 of 39) |
| [[293,68,5]]₁₉, γ<1 | — | ≥ 0.0279 | — | rigorous bound (dim 4 of 10) |
| [[2093,308,5]]₇ | — | ≥ 0.0018 | — | rigorous bound (dim 5 of 18; weak) |

## Findings

1. **A_d/k under-counts the per-output error, and the gap grows with k.** One low-weight
   logical error flips many outputs at once, because the logical rows h_i are dense. The exact
   per-output coefficient is 2.6–5× the A_d/k value on the small codes. On the flagship
   [[519,106,5]]₅ it is **~85×**: a typical weight-5 error corrupts ~85 of the 106 outputs.
   Every δ_out and cost figure built on `distillation.delta_out_avg` (eq 39) is optimistic
   per output. The STATUS.md "δ_out≈8e-18" example is one such figure. γ, d and A_d are
   unaffected.
2. **A_d optimization does not move the threshold.** Our A_d-optimized codes lower the
   leading coefficient (per output: 1.16× for [[206,37,4]]₃, 1.89× for [[112,13,3]]₅, 1.09× for
   the [[519,106,5]]₅ flagship), but ε\* agrees to 3–4 digits (flagship: 0.043139 for both). ε\* is set by the high-weight part of the weight distribution. The
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

## Multi-round protocols: what the codes actually cost

`qmsd/protocol.py` / `python -m qmsd protocol --p P --eps0 E --target T...` searches every code
sequence of depth ≤ 6 drawn from a library of codes with exact per-output curves:
- **Catalog codes.** Our tier-1/2 codes.
- **The paper's codes.** Its cyclic codes (App. C), plus the shortest RS and RM codes.
- **arXiv:2610.29457 codes (added 2026-10-22).** The doubled codes of its Table 3, and for p = 3
  the Prakash–Saha [[6m+2, 3m−2, 2]]₃ codes, m ≤ 12 (`qmsd/doubling.py`).
- **Excluded:** codes with a non-cube λ, e.g. [[12,1,3]]₇, whose output is in another Clifford class.

**Cost model.** The cost is the expected number of raw states per final output,
Π nᵢ/(kᵢ·P_acc,ᵢ). Each round takes its inputs from distinct blocks of the previous round, so
the per-output marginal is the correct input.

| p | ε₀ → target | best protocol (cost) | best with k=1 codes only | recursive smallest code |
|---|---|---|---|---|
| 5 | 1e-2 → 1e-15 | [[4,1,2]] → [[519,106,5]] (**22**) | 106 ([[11,1,3]] → [[112,13,3]], no 519) | 267 ([[4,1,2]]⁴) |
| 5 | 1e-2 → 1e-30 | [[4,1,2]]² → [[519,106,5]] (**82**) | 314 (no 519) | 1,066 |
| 5 | 1e-3 → 1e-10 | [[519,106,5]] (**8.2**) | 16 | 16 |
| 5 | 1e-3 → 1e-30 | [[519,106,5]]² (**40**) | 138 (no 519) | 257 |
| 7 | 1e-2 → 1e-10 | [[28,1,5]] cyc (**37**) | same | 230 ([[6,1,2]]³) |
| 7 | 1e-2 → 1e-30 | [[7,1,3]] cyc → [[28,1,5]] cyc (**210**) | same | 8,270 |
| 7 | 1e-3 → 1e-10 | [[20,1,4]] doubled (**20.4**) | same | |
| 7 | 1e-3 → 1e-30 | [[7,1,3]] cyc → [[20,1,4]] doubled (**141**; 175 before the doubled code) | same | 1,304 |

**p = 3 after adding the Prakash–Saha codes (2026-10-22).** These d = 2 codes have rate up to
~1/2, so they replace both the recursive and the k = 1 baselines. The earlier p = 3 rows
(e.g. 1,818 → 234 at 1e-2 → 1e-30) compared against a baseline without them and overstated the
gain about 10×.

| ε₀ → target | best protocol (cost) | standard codes only | before PS codes were in the library |
|---|---|---|---|
| 1e-2 → 1e-10 | [[20,7,2]] PS → [[206,37,4]] (**26.3**) | 40.3 | |
| 1e-2 → 1e-15 | [[20,7,2]] → [[44,19,2]] → [[206,37,4]] (**48.7**) | 58.9 | 70 |
| 1e-2 → 1e-30 | [[14,4,2]] → [[74,34,2]] → [[206,37,4]] → [[74,34,2]] (**113.9**) | 155.8 | 234 |
| 1e-3 → 1e-10 | [[206,37,4]] (**6.8**) | 11.2 | 6.8 |
| 1e-3 → 1e-15 | [[74,34,2]] → [[206,37,4]] (**13.2**) | 16.6 | |
| 1e-3 → 1e-30 | [[56,25,2]] → [[206,37,4]] → [[74,34,2]] (**29.0**) | 43.6 | 32 |

For p = 5, the "k=1 only" column is the best protocol without the flagship; it may use k>1
codes such as [[112,13,3]]₅.

**Findings:**
- **High-rate search codes win, even though their thresholds are low.** One round of a small
  high-threshold code ([[8,1,2]]₃, [[4,1,2]]₅, [[11,1,3]]₅) brings ε from 1e-2 to about 1e-4.
  A high-rate code then finishes the job.
  - At p = 3 the Prakash–Saha d = 2 codes do most of the work. [[206,37,4]]₃ is the only
    search code in the optimal protocols, and it lowers cost 1.2–1.6× below the standard codes.
  - At p = 5 the γ<1 flagship [[519,106,5]]₅ gives the cheapest protocol at every target:
    3–13× below recursive RS. This holds with its exact per-output error, which is ~85×
    worse than A_d/k.
- **γ is not a good proxy for cost at these targets.** The p = 3 winners have γ>1. What matters
  is n/k together with how far one round pushes ε past the target.
- **At p = 7 the paper's cyclic codes and the doubled [[20,1,4]]₇ dominate.** [[20,1,4]]₇ has
  non-cube coefficients, so its inputs include U^η states, which (as in arXiv:2610.29457) we
  count as one raw state each. Our [[44,5,3]]₇ codes never appear in an optimal protocol. The
  weighted-support cap codes below ([[263,24,5]]₇) are the first p = 7 search codes that beat this
  library, by ~2× in the leading-order surrogate.
- **Pre-distillation is forced by the thresholds.** The ε\* ≈ 0.02–0.06 codes never appear
  as the first round at ε₀ = 1e-2.

**Caveats:**
- The cost counts magic states only. It ignores Clifford/encoder circuit cost, footprint and
  latency. It also ignores batch granularity: a 519-qudit block yielding 106 outputs, fed
  from distinct blocks, implies large factories.
- The flagship curve is averaged over 16 sampled outputs. Its threshold is exact to 1e-9; its
  leading coefficient is good to about ±3%.
- Below ε = 1e-7 the leading order c·ε^d is used, with relative error about n·ε.
- p = 17/19 are not in the study, because their exact curves are out of reach (bounds only).

## Large p: rigorous cost bounds without exact curves

When a code's exact curve is out of reach (p^{dim A} too large), `protocol.BoundedCurve` gives
rigorous two-sided bounds on its per-output map. The inputs are a small enumerable subcode
A′ ⊆ A, the exact d, and A_d if it is known:

| quantity | lower bound | upper bound |
|---|---|---|
| ε_out(ε) | (A_d/k)·q^d(1−ε)^{n−d} / P_acc(A′) | min(Lemma-17 curve of A′, [A_d q^d(1−ε)^{n−d} + P(A′ accepts wt ≥ d+1)] / P_acc,lo) |
| P_acc(ε) | (1−ε)^n + A_d q^d(1−ε)^{n−d} | P_acc(A′) |

Notation: q = ε/(p−1). If A_d is unknown, p−1 replaces it in the lower bound and the
truncation term of the upper bound is dropped.

**From curves to costs.** Running the protocol search on the pessimistic curves gives an upper
bound on the optimal cost; running it on the optimistic curves gives a lower bound
(`protocol_bounds`, or `python -m qmsd protocol --p P --bounded <artifact_id>`). Two facts
make this valid:
- P_acc decreasing in ε is a theorem.
- ε_out increasing in ε is assumed. It holds on every exact curve we checked.

**Validation.** On [[72,9,3]]₃, [[112,13,3]]₅ and [[200,43,3]]₃, deliberately small subcodes
gave bounds that bracket the exact curves everywhere tested. The cost brackets contain the
exact optimum, and are tight to 0.01% when the code enters at ε ≲ 1e-3
(`tests/test_protocol.py`).

**Limitation.** The threshold gets only a lower bound. The optimistic curve keeps just the
weight-d term, so it never crosses ε.

**Results.** All bounds are rigorous. The libraries are RS plus the [[p,1,2⌊p/6⌋+1]] cyclic
family at p = 17 and 19 (the paper searched cyclic codes only up to p = 7), and the full p = 3
library.

| code | ε₀ → target | best without it | optimal cost with it | proven saving |
|---|---|---|---|---|
| [[237,52,6]]₁₇ | 1e-2 → 1e-30 | 301 ([[16,1,6]]²) | **[85.6, 85.6]** ([[16,1,6]] → [[237,52,6]]) | 3.5× |
| [[237,52,6]]₁₇ | 1e-3 → 1e-30 | 260 | [26.3, 74.1] | ≥ 3.5× |
| [[293,68,5]]₁₉ (A_d unknown) | 1e-2 → 1e-30 | 388 | [92.9, 99.1] ([[19,1,7]] → [[293,68,5]]) | ≥ 3.9× |
| [[293,68,5]]₁₉ | 1e-3 → 1e-30 | 330 | [24.9, 79.0] | ≥ 4.2× |
| [[230,13,6]]₃ | four cases | 26–234 | lower bound = best without it, or at most 16% below | none proven: not needed at p = 3 |

At 1e-2/1e-3 → 1e-15 the p = 17/19 codes are not needed: a single RS or cyclic round
overshoots the target anyway.

## Assessing one code: `python -m qmsd assess <artifact_id> [...]`

For each catalog code this reports:
- the per-output leading order and the threshold;
- P_acc at ε = 1e-3;
- the cheapest protocol that **uses** the code, against the cheapest protocol **without** it.
  The other rounds come from `library(p)`. The baseline keeps the code's twins, e.g. the paper
  version of an A_d-optimized code.

**Exact or bounds.** The figures are exact when the curve is computable: the histograms are
stored, or the code is tier 1. Otherwise every figure is a rigorous bound or bracket, from
`BoundedCurve`. A code needs a stored puncture set and an exact distance; A_d is optional.

**Saving column.** For exact codes it is a number. For bounded codes it reads "≥" when the whole
bracket beats the baseline, and "undetermined" when the bracket straddles it.

`--all --store` writes every assessable catalog code to `qmsd/data/performance/assess/`.
The explorer shows each code's best saving over the standard codes in a column, and the
full table in its Distillation tab.

**Two baselines.** The tables use two:
- "without": all other codes, twins included. This answers whether this exact code is needed.
- "standard": RS / RM / cyclic codes, plus (since 2026-10-22) the doubled codes and, for p = 3,
  the Prakash–Saha [[6m+2, 3m−2, 2]]₃ codes. This answers how much the code beats the published
  constructions.

**Catalog run (re-assessed 2026-10-22 against the corrected baseline).** 84 codes: the 31
catalog codes and Aanha's 53. The analytic Manhattan family, the refuted cap code and the
candidates are skipped. The table gives the best saving vs standard codes, over ε₀ ∈
{1e-2, 1e-3} × target ∈ {1e-10, 1e-15, 1e-30}, depth ≤ 5; "≥" marks a rigorous bound.

| code | ε\* per output | best saving vs standard | before PS codes (2026-10-20) |
|---|---|---|---|
| [[206,37,4]]₃ (both A_d versions) | 0.062 | 1.63× (1e-3 → 1e-10) | 20× |
| [[200,43,3]]₃ | 0.028 | 1.37× (1e-3 → 1e-30) | 20× |
| [[215,28,5]]₃ | 0.104 | 1.17× | 14.6× |
| [[72,9,3]]₃, [[667,62,4]]₃, [[79,2,4]]₃, [[80,1,5]]₃ | 0.062 / 0.021 / 0.15 / 0.16 | none | 11.6× / 8.7× / 3.2× / 1.8× |
| [[230,13,6]]₃ / [[690,39,5]]₃ | ≥ 0.045 / ≥ 0.012 | none | ≥ 6.2× / ≥ 5.3× |
| [[519,106,5]]₅ (both) | 0.043 | 6.2× (1e-2 → 1e-30) | same |
| [[20,5,2]]₅ | 0.023 | 6.2× | same |
| [[293,68,5]]₁₉ | ≥ 0.028 | ≥ 4.2× | same |
| [[112,13,3]]₅ (both) | 0.077 | 3.7× | same |
| [[237,52,6]]₁₇ | ≥ 0.034 | ≥ 3.5× | same |
| [[24,1,3]]₅ | 0.26 | 1.5× | same |
| [[2093,308,5]]₇ | ≥ 0.0018 | **4.1× (1e-2 → 1e-30), surrogate** (weight-5 census below) | same |
| [[44,5,3]]₇ (all 10), [[124,1,4]]₅ | 0.173 / ≥ 0.11 | none | same |

Reading:
- At p = 3 the Prakash–Saha codes remove almost all of the gain. Only [[206,37,4]]₃,
  [[200,43,3]]₃ and [[215,28,5]]₃ still lower cost, by 1.2–1.6×, as middle or late rounds after
  d = 2 front ends. The earlier 3–20× figures came from a baseline without these codes.
- p = 5 and the large-p codes are unchanged. The only new p = 5 standard code, the doubled
  [[33,1,5]]₅, never enters an optimal protocol.
- At p = 7 the paper's cyclic codes and the doubled [[20,1,4]]₇ are better than our catalog
  codes; the weighted-support caps (next section but one) change this.
- Twins, such as the A_d-optimized versions, give the same saving, so A_d optimization does not
  show up in cost.
- The p = 23 Manhattan [[17,6,3]]₂₃ shows 24.7×. This is an artifact of a thin baseline: the
  p = 23 standard library has only [[23,1,7]] and no multi-output RS codes.

**Aanha's imported codes (2026-10-21).** The 53 records merged from `aanha/results` were
independently re-certified:
- Every code was rebuilt from its puncture set, or from its weight-≤w point set for the small
  Manhattan codes, which `performance.puncture_columns_of` now rebuilds explicitly. n, k, full
  rank, γ and d all match. d comes from the certified MITM search, or from Theorem 4 for the
  Manhattan codes with d ≥ 9.
- Exact A_d (`structured_ad`) was written into the 24 d ≤ 4 records that lacked it. It agrees
  with all 10 values she had stored.
- 11 of her 12 seeded searches rerun to identical output: 8 locally, and m6_flatspread,
  m6_seed1 and `arc_climb` (~29 h) on lenore. The 1000-trial m = 6 run needs ~22 h serially
  and timed out at 12 h; a 4-process rerun (6 h, different random stream, so not bit-identical)
  gives codes of the same quality: 6 of its top 10 (n, k, d) match hers, best γ 1.721
  ([[686,43,5]]₃) against her 1.727 ([[668,61,4]]₃). [[686,43,5]]₃ is now in the catalog
  (`imported-p3-n686-k43-8881d380242d`), independently re-certified; it is the best-γ qutrit
  d = 5 code we have (vs the paper's [[690,39,5]]₃, γ 1.785), but like the other m = 6 codes
  it saves nothing against the standard library (ε\* ≥ 0.012, P_acc(1e-3) ≈ 0.50; cheapest
  protocol using it ≥ 69 at 1e-3 → 1e-30 vs 43.6 standard). Exact A_5 is out of reach (> 100 GB).

| her codes | ε\* per output | best saving vs standard (2026-10-22) | before PS codes |
|---|---|---|---|
| [[202,41,3]]₃ … [[204,39,3]]₃ (5 codes, exact) | 0.040 – 0.051 | 1.22 – 1.29× (1e-3 → 1e-30) | 17.8 – 18.9× |
| [[214,29,4]]₃ … [[217,26,4]]₃ (14, incl. 10 × [[215,28,4]]₃) | ≥ 0.048 | ≥ 1.07 – 1.22× (1e-3 → 1e-10) | 13.2 – 15.0× |
| [[223,20,5]]₃ | ≥ 0.046 | none | ≥ 8.4× |
| [[668,61,4]]₃ … [[679,50,4]]₃ (12) | ≥ 0.012 | none | ≥ 7.6 – 9.4× |
| [[689,40,5]]₃ … [[696,33,5]]₃ (6) | ≥ 0.012 | none | ≥ 4.4 – 5.4× |
| [[555,70,3]]₅ / [[590,35,6]]₅ (Manhattan) | ≥ 0.013 / ≥ 0.015 | ≥ 2.9× / ≥ 1.4× | same |
| [[610,15,9]]₅, [[620,5,12]]₅, [[624,1,14]]₅ (Manhattan) | ≥ 0.014 – 0.020 | none | undetermined / none |
| p = 7: [[322,21,5]]₇ (2), [[325,18,5]]₇, [[323,20,3]]₇, [[333,10,4]]₇, [[339,4,5]]₇ (5) | ≥ 0.030 – 0.032 | none | undetermined or none |

Reading:
- Her m = 5 qutrit codes save 1.1–1.3×, in the same range as the paper's [[206,37,4]]₃ and
  [[200,43,3]]₃. Her m = 6 qutrit codes save nothing against the corrected baseline.
- None of them beats the paper's [[206,37,4]]₃. For example, requiring [[215,28,4]]₃ costs
  39.7 at 1e-3 → 1e-30 against 29.0 without it, so they are not needed in the cheapest catalog
  protocol.
- At p = 7 her codes do not enter any optimal protocol.

**Search objective (Option C, `qmsd/objective.py`).** Because γ does not predict these savings,
`marginal_saving(candidate, standard_library(p))` scores a code by cost(standard) /
cost(standard + candidate). Candidates use a leading-order surrogate built from `structured_ad`
(per-output L_d = Σ_w (w − d)·hist[w] / k). On the 10 qutrit catalog codes the surrogate matches
the exact saving to ≤ 3%: [[206,37,4]]₃ 1.503 / 1.503, [[215,28,5]]₃ 1.108 / 1.108,
[[204,39,3]]₃ 1.091 / 1.119, [[203,40,3]]₃ 1.000 / 1.025, and the other six 1.000 / 1.000.
The figures here are lower than the catalog table because the library holds only the
standard codes plus the candidate. CLI: `python -m qmsd value <ids> [--surrogate]`.
`surrogate_from_matrices(G0, H1, p)` gives the same surrogate for any code from its matrices
(weight-d accepted errors enumerated once by the streamed MITM in `structured_ad`, then classified
per output); it reproduces `structured_ad` on the catalog codes and the slow streaming count on
[[235,10,5]]₇ (L_d = 25.2) in 89 s instead of 66 min.

**Weighted-support caps at p = 7 (2026-10-26, surrogate; see NEXT_STEPS.md §H).** Caps on the
support of a degree-2 charge polynomial (arXiv:2610.29457 §4) shorten F₇³ from 343 to 287 points
and keep d = 5 up to k = 24. Scored with the surrogate against the p = 7 standard library:

| code | charge | c_ε | 1e-2 → 1e-30 | 1e-3 → 1e-15 | 1e-3 → 1e-30 |
|---|---|---|---|---|---|
| **[[263,24,5]]₇** | x₁²+x₂²+x₃²+1, outputs in the 2·cube class | 0.0063 | **2.55×** (210 → 82) | **2.02×** (28.8 → 14.3) | **2.13×** (141 → 66) |
| [[270,24,5]]₇ | x₁²+x₂²+x₃ | 0.0096 | 2.49× | 1.95× | 2.07× |
| [[266,21,5]]₇ | x₁²+x₂²+x₃²+1 | 0.0062 | 2.21× | 1.74× | 1.84× |
| [[319,24,5]]₇ | unit (unweighted) | 0.0177 | 2.11× | 1.57× | 1.50× |
| [[235,10,5]]₇ | x₁(x₁−1) | 0.0031 | 1.19× | none | none |
| **[[2093,308,5]]₇** (m = 4, unit, free-form) | — | 0.172 | **4.11×** (210 → 51) | none | **2.94×** (141 → 48) |
| **[[1746,270,6]]₇** (m = 4, sphere charge) | x₁²+x₂²+x₃²+x₄² | 2.9e-5 | **4.32×** (210 → 48.7) | none | **3.60×** (141 → 39.2) |

The last row is the m = 4 free-form code from the γ hunt (γ 1.19), scored with the thread-parallel
weight-5 census (`mindist_nb.count_weight5_per_output`, 3.1 h on lenore under load): it has only
A_5 = 1560 accepted weight-5 errors (260 scalar classes), but each hits 86% of the 308 outputs,
so L_d = 1338 per output and c_ε = 0.17. Its protocols are [[7,1,3]]₇ → [[2093,308,5]]₇ from both
1e-2 and 1e-3, i.e. the code is applied at ε ≈ 1e-6 to 1e-9 (n·ε ≤ 0.002, leading order
essentially exact); at 1e-2 → 1e-30 the output is 1.8e-31, a 5.6× margin. It is not used at
1e-3 → 1e-15 because P_acc(1e-3) = 0.12 makes a direct round cost 55 against 28.8.

**[[1756,260,≥7]]₇ (m = 4, sphere charge, 2026-10-27): γ ≤ 0.982, the first γ < 1 search code
at p = 7** (three-part certificate in NEXT_STEPS.md §H). Cost pending its weight-7 census.

**[[1746,270,6]]₇ (m = 4, sphere charge, 2026-10-26)**: γ = 1.042, the best p = 7 γ before that
(certificate in NEXT_STEPS.md §H and `runs/weighted/p7_m4_sphere_1746_270_cols.json`). Weight-6
census on its class block (`count_weight6_per_output`, 3+3, 30 min): only A_6 = 24 words
(4 scalar classes), each hitting ~6% of the outputs, L_6 = 1.33, c_ε = 2.9e-5. Protocols:
[[7,1,3]]_7 cyc-B → [[1746,270,6]]_7 (1e-2 → 1e-30, cost 48.7) and [[6,1,?]]_7 RS → [[1746,270,6]]_7
(1e-3 → 1e-30, cost 39.2); the code is applied at ε ≈ 1e-6, output ≈ 3.0e-41. Words that
mix in non-class columns are assumed absent (35 extra checks; random-model expectation ~1e-9),
so this is a surrogate on a lower bound of L_6, and at this n·ε the next order is negligible.
Data: `runs/weighted/p7_m4_sphere_1746_270_6_census.json`.

Typical protocols: [[7,1,3]]₇ → [[263,24,5]]₇ (1e-2 → 1e-30) and [[263,24,5]]₇ alone
(1e-3 → 1e-15). The code is always applied at input noise ≤ 1e-3, i.e. n·ε ≤ 0.26, where the
leading order is within 4% of the exact ε_out on the [[44,5,3]]₇ codes and the margins to the
targets are ≥ 150×. This is a surrogate result, not an exact or rigorously bounded one: the exact
per-output curve needs 7^{~32} words, and the dim-5 Lemma-17 bound used for Aanha's
[[322,21,5]]₇ was too loose to decide. Still, it is the first p = 7 search code that enters an
optimal protocol, and the weighting is worth ~1.2–1.4× over the unit-space code of the same k and d.

## Engines / tiers (`performance.assign_tier`)

| tier | when | engine |
|---|---|---|
| 1 | p^{dim A}·k ≤ 3e8 | local numpy (`weightdist.coset_weight_enumerators`) |
| 2 | ≤ 1e13 | `dwq_coset` C/OpenMP kernel (`dist_weightdist/coset_run.py`), lenore |
| 3 | larger | rigorous Lemma-17 subcode lower bound (`threshold_lower_bound_subcode`) |
| 4 | analytic Manhattan codes | none (no stored puncture set) |

The block error needs p^{dim A + k} words, so it is feasible only for small k.

## Not yet done

- [[293,68,5]]₁₉ exactly (19¹⁰ ≈ 6e12 words per output; ~20 h on lenore with 4 outputs).
  [[230,13,6]]₃, [[690,39,5]]₃ and [[2093,308,5]]₇ are beyond exact reach (bounds above). Two
  shortcuts were tested and fail: uniform codeword sampling misses the rare low-weight
  elements that carry the logical signal, and truncating the error-weight series converges
  only once weights ~d+14 are included.
- Exact per-output leading order for the large-p codes (needs the explicit min-weight codewords).
- Protocol costs that include Clifford circuit cost and factory footprint (the current model
  counts magic states only).
