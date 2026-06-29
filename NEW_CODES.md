# New γ<1 codes found by search (this toolkit)

Codes discovered by the `qmsd` search that are **not** in arXiv:2510.10852 (which searched only
p=3,5 numerically and gave analytic codes for higher p). Distances certified by the
adversarially-validated meet-in-the-middle routine (`qmsd.mindist.min_dependent_columns`).

## p = 17

### `[[234, 55, 5]]₁₇`  —  γ = 0.8997  (found 2026-06-27)

- Construction: puncture `RM₁₇(10, 2)` (m=2, r = r_max = 10; triorthogonal since 3·10 = 30 < 2·16 = 32).
- Parameters: n = 234, k = 55, d = 5; full rank; γ = log(234/55)/log(5) = 0.8997.
- Significance: ~330× smaller than the paper's analytic p=17 code `[[77540,5981,15]]`; below the
  Reed–Solomon threshold (p=17 < 23) so not reachable by m=1 RS. First small certified γ<1 code at p=17.
- A_d: not yet computed (pending the MacWilliams-from-small-dual engine).

Puncture columns (1-indexed, per the Appendix-C convention):

```
[1, 3, 6, 9, 12, 16, 27, 28, 29, 45, 48, 57, 58, 70, 74, 79, 80, 85, 90, 94, 96, 104,
 107, 110, 111, 119, 120, 122, 131, 133, 140, 142, 175, 186, 194, 216, 218, 225, 227,
 237, 240, 251, 252, 257, 258, 259, 262, 263, 264, 265, 276, 280, 283, 287, 288]
```

Reproduce / verify:

```python
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.mindist import min_dependent_columns
punc = [1,3,6,9,12,16,27,28,29,45,48,57,58,70,74,79,80,85,90,94,96,104,107,110,111,119,
        120,122,131,133,140,142,175,186,194,216,218,225,227,237,240,251,252,257,258,259,
        262,263,264,265,276,280,283,287,288]
b = build_triorthogonal_code(17, 2, 10, punc)
assert b["full_rank"]
assert min_dependent_columns(b["X_stab"], 17, d_max=7) == 5   # certified distance
# n = 17**2 - len(punc) = 234, k = len(punc) = 55  ->  gamma = 0.8997
```

Note: many other γ<1 codes exist at p=17, m=2 in the high-k window (k≈53–59, d=4–5); this is the
lowest-γ one found so far. p=19 (more headroom) and pushing for d=6 are expected to go lower.

## p = 3 (qutrit), m = 7 — cap puncturing  ⟶  REFUTED

### ~~`[[1968, 219, ≥10]]₃` (γ = 0.9536)~~ — **REFUTED: actual d = 1 ⇒ γ ≥ 1** (2026-06-27)

**This code does NOT achieve γ<1.** Validation (`CAP_VALIDATION.md`) established its true distance is
**1**: the shortened generator `G0` (55×1968) has **83 all-zero columns**, each an exact weight-1
*logical* codeword of `G0^perp` (an `RM₃(9,7)` codeword of weight 77–145 whose support sits almost
entirely inside the cap, leaving 1 point outside). The `≥10` was only a *minimum-weight-class* bound;
a high-weight codeword concentrating on the cap collapses the real distance to 1. Verified exactly and
re-checkably (83 zero columns, syndromes = 0) — not probabilistic.

The full-rank **206-cap** alternative `[[1981,206]]₃` fails **identically** (140 zero columns → d=1).
**Conclusion: cap-puncturing for m=7 qutrit γ<1 is ruled out** — any large (~200-pt) algebraic cap in
AG(7,3) makes pervasively many points' `RM₃(4,7)` evaluations dependent on the cap, producing
degenerate (zero) coordinates = weight-1 logicals. (Earlier ad-hoc "0 zero columns" checks were an
off-by-one indexing artifact; the validated `cap_validate.trivial_low_weight` is authoritative.)

Historical construction (for the record): puncture `RM₃(4,7)` at a 219-pt cap (maximal full-rank
sub-cap of a rank-deficient 236-cap); the min-weight (weight-18 = 2 parallel 2-flats) bound
`max|cap ∩ 2-flat| = 4 ⇒ ≥10` is real but is NOT the code distance. Files: `CAP_QUTRIT_RESULT.md`
(original candidate), `CAP_VALIDATION.md` (refutation), `cap_validate/` (the validated finder).

## p = 3 (qutrit), m = 5 — A_d-optimized code (strictly better than the paper)

### `[[206, 37, 4]]₃` — A_d = 572 (paper: 880)   (2026-06-27, verified)

- Same parameters and cost as the paper's `[[206,37,4]]₃` (n=206, k=37, d=4, C=22.1 at δ_in=0.01),
  but **A_d = 572 vs 880** ⇒ ~1.54× lower output error per round, **compounding** to 8.6× (2 rounds)
  / ~8500× (3 rounds) through the d=4 recursion. Strictly better: same n, k, d, C; lower δ_out.
- Found by deep cap-set search at m=5, k=37: a cap (no 3 collinear punctures) minimizes the affine
  lines that spawn low-weight codewords, so cap-set sampling lowers A_d while keeping d maximal —
  **all 3 d=4 codes found beat 880**.
- A_d = 572 **certified two independent ways**: the MacWilliams engine (`qmsd.weightdist`) and the
  logical count `qmsd.distance.A_d_logical_Z` (= 572, with the Gp stabilizer filter).
- Puncture columns (1-indexed), saved in `qutrit_Ad572.json`:
  [13,22,31,34,35,40,53,61,70,78,80,81,90,91,95,96,109,112,118,121,122,131,136,156,157,180,185,186,189,199,203,212,216,219,228,234,242]

See `QUTRIT_PARETO.md` for the full (γ, C, A_d) study.

## p ≥ 23 (Reed–Solomon, m = 1) — the optimal small γ<1 code FAMILY   (2026-06-27, verified)

The paper lists only `[[17,6,3]]₂₃` as "the smallest known γ<1 code." But at m=1 the quantum code is a
*punctured MDS* code, so the distance is the closed form **d = r_max − k + 2** (`r_max = ⌊(p−2)/3⌋`),
**independent of which columns are punctured** — verified against qmsd (d *and* A_d match exactly for the
p=23 cases). So the whole γ<1 family is exactly enumerable, with **no search and no distance bottleneck**:

- Every prime p≥23 has a FAMILY of γ<1 codes `[[p−k, k, r_max−k+2]]` over the γ<1 range of k — a clean
  (γ, n, C) Pareto trade per prime.
- **At p=23, `[[18,5,4]]₂₃` (γ=0.924) beats the paper's `[[17,6,3]]` (γ=0.948) on γ** — the paper optimized
  smallest block n, not γ. (The full p=23 family: `[[19,4,5]]` 0.97, `[[18,5,4]]` 0.92, `[[17,6,3]]` 0.95.)

Min-γ code per prime (all certified by the MDS closed form; cross-checked in qmsd):

| p | code | γ | p | code | γ | p | code | γ |
|---|---|---|---|---|---|---|---|---|
| 23 | `[[18,5,4]]` | 0.924 | 47 | `[[36,11,6]]` | 0.662 | 73 | `[[56,17,8]]` | 0.573 |
| 29 | `[[22,7,4]]` | 0.826 | 53 | `[[40,13,6]]` | 0.627 | 79 | `[[60,19,8]]` | 0.553 |
| 31 | `[[25,6,5]]` | 0.887 | 59 | `[[45,14,7]]` | 0.600 | 83 | `[[62,21,8]]` | 0.521 |
| 37 | `[[29,8,5]]` | 0.800 | 61 | `[[47,14,7]]` | 0.622 | 89 | `[[67,22,9]]` | 0.507 |
| 41 | `[[31,10,5]]` | 0.703 | 67 | `[[51,16,7]]` | 0.596 | 97 | `[[73,24,9]]` | 0.506 |
| 43 | `[[34,9,6]]` | 0.742 | 71 | `[[53,18,7]]` | 0.555 | | | |

So from the paper's single γ=0.95 code to a certified family down to **γ≈0.51** (block ≤ 73), tracking
γ₀(p)~1/ln p. **`A_d` note (corrected 2026-06-28):** MDS makes `A_d` large and structurally fixed
(`A_d = C(n,d)(p−1)`, e.g. 67320 for `[[18,5,4]]₂₃`), but this is **NOT a suppression downside** — verified
that per-round δ_out *improves* across the ladder (5.2e-10 @p=23 → 5.6e-25 @p=97) because the larger d
(δ_in^d) and p^(d−1) dominate the `A_d` growth. The genuine MDS limitation is the **lost `A_d`-optimization
freedom** (puncture-invariant ⇒ can't cut `A_d` the way the qutrit cap-set did 880→572), which only matters
in the small-p / within-prime regime. Reproduce: `rs_family.py`; cross-prime analysis in `ARCHITECTURE_DIMENSION.md`.
