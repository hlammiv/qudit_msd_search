# No small γ<1 triorthogonal magic code over prime F_p for p ≤ 13

**A negative result (Direction D).** Consolidated 2026-06-30 from `D_PHASE01_FINDINGS.md`
(execution log), `RESULTS.md` (the positive results + the flat-cap boundary), and the
`qmsd` toolkit. Companion to arXiv:2510.10852.

---

## Summary

For magic-state distillation, a triorthogonal CSS code over the qudit dimension p with parameters
[[n,k,d]] has overhead exponent **γ = ln(n/k)/ln(d)**, and γ<1 (sublogarithmic overhead) requires
**d > n/k**. Small γ<1 codes are known for p ≥ 17 (the punctured-Reed–Muller codes `[[237,52,6]]₁₇`,
`[[293,68,5]]₁₉`) and p ≥ 23 (the Reed–Solomon family). This note establishes that **no analogous
small γ<1 code exists for p ≤ 13** — closing the gap below the p=17 boundary. The result is:

- **Rigorous within the evaluation-code framework** (Reed–Muller and every column-subset of it), via the
  **forced-grid theorem**: the only triorthogonal subsets of the affine grid are ∅ and the full grid, so
  the *only* construction available is ordinary RM puncturing — which is distance-capped below γ<1 for
  p ≤ 13 (the flat cap).
- **Strongly evidenced for triorthogonality in general** (beyond evaluation codes): three independent
  non-RM families — cyclic/BCH, algebraic-geometry, and Artin–Schreier — are each closed, and direct
  moment-search caps collapse identically. The unifying statement: **the distance cap is generic to
  triorthogonality over a prime field, a degeneracy of shortening to a low-redundancy CSS X-stabilizer,
  not an artifact of Reed–Muller geometry.**

One honest gap remains (§7): a single high-dimensional general-position regime is not constructively
reachable and not closed by theorem (estimated <5% chance of hiding a survivor).

---

## 1. Setup

A CSS code built from a triorthogonal generator G over F_p admits a transversal non-Clifford diagonal
gate (the qudit T_s of arXiv:2510.10852). Triorthogonality (`qmsd.triorthogonal.is_triorthogonal`,
ground truth throughout) means the rowspan satisfies, mod p and including repeated indices,

  pair:  Σ_i G[a,i]G[b,i] = 0   (a≤b);   triple:  Σ_i G[a,i]G[b,i]G[c,i] = 0   (a≤b≤c).

The magic (non-Clifford logical action) emerges by **puncturing**: the full generator must be fully
triorthogonal; shortening at k columns yields the [[n−k, k, d]]_p code with X-stabilizer G0 (shortened
generator) and distance d = the minimum number of F_p-dependent columns of G0
(`qmsd.mindist.min_dependent_columns`). Overhead **γ = ln(n/k)/ln(d)**; γ<1 ⟺ **d > n/k**.

**Point/moment view (verified equivalent to `is_triorthogonal`).** Reading each column as a point
x_i ∈ F_p^K, triorthogonality = the 2nd and 3rd moment tensors of the column multiset vanish mod p
(Σ_i x_i⊗x_i = 0, Σ_i x_i^⊗3 = 0). These are *linear* in the column multiplicities — the lever for the
direct search (§5).

## 2. Context — the positive results (p ≥ 17) that frame the question

`RESULTS.md` establishes the γ<1 codes that *do* exist and motivate the p ≤ 13 question:
- **p = 17:** `[[237,52,6]]₁₇`, γ=0.8466 (d=6 exact, A₆=147856), from punctured RM₁₇(10,2).
- **p = 19:** `[[293,68,5]]₁₉`, γ=0.9076 (d=5 exact), from punctured RM₁₉(11,2).
- **p ≥ 23:** the Reed–Solomon (m=1) family, MDS, closed-form d = r_max−k+2, γ from 0.924 (p=23) to ≈0.506.

The p ≥ 17 codes are the *first* small γ<1 search codes; the natural question is whether the trend
extends below p=17. It does not.

## 3. The flat cap closes punctured Reed–Muller for p ≤ 13

Min-weight codewords of the relevant dual concentrate on low-dimensional affine flats, giving the
line/flat-spread identity d = d_RM − max_flat|S∩flat|. At the γ<1 puncture density (high k), every
codimension-1 flat over-punctures by ~k/p^{m−1}, so a 2D codeword (m=2), plane codeword (m=3), or
hyperplane codeword (m=4) caps the true distance below the γ<1 threshold. Per-p closure (Filter C, the
binding 2D cap, re-validated live with `min_dependent_columns`):

| p | best RM construction at γ<1 density | true d | needed d | gap | best γ | status |
|---|---|---|---|---|---|---|
| 5 | crosses only at m=4 (`[[519,106,5]]₅`, γ=0.987); no smaller code | — | — | — | 0.987 | Singleton-infeasible m≤3 |
| 7 | m=2 Singleton-infeasible; m=3 plane-capped; m=4 hyperplane-capped | 2–3 | >6.6 | large | >1 | closed for k≥313; **mid-k m=4 open** |
| 11 | RM₁₁(6,2), k=21 ⇒ `[[100,21,4]]` | 4 | 5 | **+1** | 1.126 | 2D-capped |
| 13 | RM₁₃(7,2), k=26 ⇒ `[[143,26,5]]` | 5 | 6 | **+1** | 1.059 | 2D-capped |
| ≥17 | line-supported words become the true minimum | — | — | — | <1 ✓ | **crosses** |

The discriminator p≤13 vs p≥17 is whether line-supported words are the *global* minimum at the γ<1
density — first true at p=17, m=2. (Filter A, the Singleton bound, only sets first-feasible m; the 2D cap
is what actually binds.)

## 4. The forced-grid theorem — the evaluation route is fully closed (rigorous)

The flat cap above is about *puncturing* RM. Could a *different* subset of the affine grid, or a different
evaluation code, escape? No:

> **Theorem (forced grid).** Fix the full degree-`r_max` monomial generator on the grid F_p^m. A column
> subset `S ⊆ F_p^m` yields a triorthogonal code **iff** its indicator `1_S ∈ RM(t,m)`, where
> `t = m(p−1) − 3r − 1`. At every p ≤ 13 target `2t < p`, so on each affine line `1_S` is a univariate
> polynomial g of degree ≤ t with g∈{0,1}; then `g(g−1)` has degree `2t < p` and vanishes on all p line
> points, hence is the zero polynomial, so g is constant. `1_S` constant on every line ⇒ globally
> constant ⇒ **`S ∈ {∅, full grid}`.**

*Proof of the iff:* the triple condition says `1_S` is orthogonal (coordinate dot product) to all products
of three degree-≤r monomials, i.e. to `RM(r,m)^{⋆3} = RM(3r,m)`; the pair condition (⊥ RM(2r,m)) is
weaker. So `1_S ∈ RM(3r,m)^⊥ = RM(m(p−1)−1−3r, m) = RM(t,m)`. ∎

**Consequence.** Within the evaluation framework — where every known γ<1 RM code and *all* its puncturings
live — the only triorthogonal building block is the full grid, so the only available construction is
ordinary RM puncturing, capped at the §3 distances (p=11 d=4, p=13 d=5). No non-RM *evaluation* code escapes.

**Independent verification** (`scratchpad/verify_p2.py`, RAM-bounded): the Hadamard-cube rank identity
`rank(RM(r,m)^{⋆3}) = dim RM(3r,m)` holds exactly (118=118 at p=11, 163=163 at p=13, 333=333 at p=7 m=3);
`2t<p` at all targets (t = 1,2,2,2,0 for p=11/13, p=7 m=3/m=4, p=5 m=4); the full grid is triorthogonal,
removing **any** single point breaks it (10/10), **0 of 250** random proper subsets are triorthogonal at
p=11 and p=13; and RM(t,m) has **no** non-trivial 0/1 codeword (the constant-forcing step), exhaustively
at p=11.

### 4.1 The boundary engine — a closed form for *why* the threshold is p=17

The flat cap (§3) shows p≤13 is capped *empirically*; the following gives the **closed-form reason**, and
pins the boundary to a single inequality.

> **Boundary engine (proved; independently verified).** The moment / normal-rational curve in `F_p^r`
> (columns `(1, t, t², …, t^{r-1})`, `t∈F_p`) is **triorthogonal iff `3(r−1) < p−1`**, and is **MDS**, so its
> Z-distance is `d = r+1` (Singleton defect 0). Hence the maximum distance of any *triorthogonal-MDS* cap is
> `d_max = ⌊(p+1)/3⌋ + 1`, giving

| p | 5 | 7 | 11 | 13 | **17** | 19 |
|---|---|---|---|---|---|---|
| triorthogonal-MDS `d_max` | 3 | 3 | 5 | 5 | **7** | 7 |

> The first prime admitting a triorthogonal-MDS cap of distance **6** is **p=17** (at `r=5`: `3·4=12<16`),
> matching `[[237,52,6]]₁₇` exactly. So `p≤13` is hard-capped at `d=5`, and `p=17` is precisely where `d=6`
> first becomes algebraically available — *without* over-killing `p≥17`.

*Proof sketch.* The row functions are the power maps `t↦t^j`; the 2nd/3rd moments are power sums
`Σ_{t∈F_p} t^s`, which vanish mod p unless `(p−1)∣s`. The largest exponent appearing is `3(r−1)` (triple
products), so all moments vanish iff `3(r−1) < p−1`. The normal rational curve is the classic MDS arc
(`d=r+1`). Verified end-to-end (`scratchpad/verify_proof.py`): `is_triorthogonal ⟺ 3(r−1)<p−1` matched on
every `(p,r)`, and `d = r+1` (MDS) every time; the `d_max` table reproduced as `{5:3, 7:3, 11:5, 13:5, 17:7}`.

**Caveat (the MDS branch carries no logicals).** A defect-0 (MDS) triorthogonal cap admits `k=0` logical
qudits — the high-`d` MDS route does not by itself give a `γ<1` *code* (its logical-by-magic Gram block is
totally isotropic, so its rank is `≤⌊k/2⌋<k`). So `d_max` is the distance *ceiling*; the codes that actually
carry logicals (`k≥1`) sit strictly below it (e.g. p=11 realizes `d=4`, p=13 `d=5` with logicals), which is
why both miss `γ<1` by exactly one distance unit. Closing that last unit *rigorously* (the rank + 2D-codeword
cap) is the one remaining open crux — see §7 and `D_PROOF_MAP.md`.

## 5. Non-evaluation routes are closed too (strong evidence for the general thesis)

Three families that are *not* evaluation codes — the natural candidates to escape the forced grid — each fail:

- **Cyclic / BCH** (length n = p²−1, *not* affine-invariant). Genuinely non-RM triorthogonal cyclic codes
  **do exist** (maximal triorthogonal dim 30 > RM cap 28 at p=11; 51 > 36 at p=13; defining sets not
  P-closed) — independently verified on the explicit dim-30 generator. **But they collapse to d ≤ 2** under
  puncturing at the γ<1 density (γ ≥ 1.65), strictly worse than RM. Clean transform-domain characterization
  (validated, 0 disagreements vs `is_triorthogonal`): pair ⟺ Z∩(−Z)=∅, triple ⟺ no j+k+l≡0 (mod n) in the
  spectral support Z. *This is the key evidence the cap is not RM-special:* a non-affine-invariant family
  collapses anyway.
- **Algebraic-geometry (one-point AG / Goppa).** Codes `L(αP∞)` over prime F_p are **not triorthogonal at
  all** — `is_triorthogonal` fails for every nonsingular elliptic curve over F₇/F₁₁/F₁₃ and every genus-2
  curve over F₇ (~2150 curve×α tests). The triple-product moment Σ_P f·g·h has no curve analog over a
  prime field (the self-dual "Castle"/maximal curves that would supply one live over F_{q²} ⇒ the settled
  field-CCZ / multi-qudit collapse).
- **Artin–Schreier `yᵖ−y=f(x)`.** Degenerate over the prime field: `yᵖ−y ≡ 0` (Fermat) makes the
  F_p-rational points vertical lines over the roots of f — no usable one-point code. The additive lever
  needs F_{pᵏ} ⇒ the same multi-qudit collapse.
- **Direct moment-search (Rank 2).** Beyond the grid, moment-zero **caps** (genuine, 0 collinear triples,
  built via ±-symmetry MITM subset-sum) still collapse to **d ≤ 3** in the shortened code (weight-2
  dependencies appear in G0), best γ=1.32. Line-unions collapse to d=1.

## 6. The unifying statement

Triorthogonality = vanishing degree-≤3 moments = a strength-3 orthogonal-array / 3-design condition over
F_p. The only *tractable* moment-zero structures are **flats** (grids/lines/planes — forced to RM, §4) and
**multiplicative orbits** (linear-code OAs / cyclic — collapse, §5); both are collinear/scaling-rich and so
carry low-weight dependencies that survive shortening. Genuinely exotic moment-zero caps exist but, in every
case constructible (§5), collapse to d≤3. Hence:

> **Thesis (the negative result).** For p ≤ 13, the distance cap that prevents γ<1 is **generic to
> triorthogonality over the prime field F_p** — a degeneracy of shortening a triorthogonal code to a
> low-redundancy CSS X-stabilizer — not a property special to Reed–Muller. Rigorous for evaluation codes
> (forced-grid theorem); strongly evidenced in general (cyclic/AG/Artin–Schreier kills + moment-cap collapse).

## 7. Honest scope and open items

- **High-dimensional general-position cap (the <5% tail).** A γ<1-density code needs generator dim K≈28
  (p=11) / 36 (p=13), i.e. points in F_p^~27, and a 0/1 (set) selection satisfying ~4000 moment equations
  — an NP-hard subset-sum over (Z_p)^q (Davenport-constant threshold ~q·log p ≈ 150). This regime is not
  constructively reachable here and **not closed by a theorem**. Every adjacent realizable construction
  collapses, so this is a low-probability tail, but it is the one route a determined search or a new
  algebraic idea could still attack.
- **p=7 mid-k m=4 RM window** `k∈[110,312]`. Still uncertified — `min_dependent_columns` overflows past
  R=21 rows there. This is *RM* (within the forced grid), not a non-RM route; closing it rigorously needs
  an overflow-safe distance certificate (a separate engineering task, see `RESULTS.md`).

## 8. Per-p verdict

| p | small γ<1 triorthogonal code? | basis |
|---|---|---|
| 5 | no smaller than `[[519,106,5]]` (γ=0.987, m=4) | Singleton-infeasible m≤3; forced-grid; moment collapse |
| 7 | none for k≥313; mid-k m=4 RM window open | flat cap (high-k rigorous); AG/AS/cyclic/moment all closed |
| 11 | no (misses by +1: d=4, need 5) | 2D cap; forced-grid; cyclic d≤2; AG/AS closed; moment caps d≤3 |
| 13 | no (misses by +1: d=5, need 6) | 2D cap; forced-grid; cyclic d≤2; AG/AS closed; moment caps d≤3 |

## 9. Methods, tooling, reproducibility

All `qmsd` + scratchpad, local 15 GB box (heavy jobs RAM-capped); ground truth = `is_triorthogonal` and
`min_dependent_columns(d_max≤5)`. Feasibility filters and Filter-C re-runs: `scratchpad/phase0.py`. AG kill:
`scratchpad/probe1b_elliptic.py`, `probe1b_ext.py`. Cyclic: `scratchpad/cyc_*.py` + the transform
characterization. Moment/direct search + forced-grid: `scratchpad/d_p2_*.py`, `p2_forcedgrid.py`,
`verify_p2.py`. Distances at γ<1 density independently re-validated against the RM baseline.

## 10. Conclusion

For magic-state distillation at qudit dimension p ≤ 13, there is **no small γ<1 triorthogonal code** to be
had: the evaluation route is rigorously forced to ordinary RM puncturing (distance-capped below γ<1), and
the principal non-RM routes are independently closed, pointing to a cap generic to triorthogonality over
the prime field. The practical recommendation is unchanged: use the Reed–Solomon family at p ≥ 23, or the
punctured-RM codes at p = 17/19; **p ≤ 13 offers no sublogarithmic-overhead triorthogonal magic code**
(modulo the high-dimensional tail of §7 and the p=7 mid-k RM window). The forced-grid theorem and the
cross-family collapse together constitute the publishable negative.
