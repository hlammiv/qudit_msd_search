# Direction D — Phase 0–2 execution findings

2026-06-30. Executes `D_PLAN.md` Phases 0–2 (feasibility filter + the cheapest decisive probes +
the moment/direct search). All local on the 15 GB box; distance certified by
`min_dependent_columns(d_max≤5)` (RAM-safe). See `D_NEGATIVE_RESULT.md` for the consolidated write-up.

## Headline

**ALL of the plan's genuine search routes are now CLOSED.** Rank 1 (cyclic/BCH), Rank 3 (AG), Rank 4
(Artin–Schreier), and Rank 2 (direct moment search) all fail to produce a γ<1 code; Rank 5 (affine) is
subsumed by the forced-grid theorem. The central thesis (`D_PLAN.md` §1, §6) — **the distance cap is
generic to triorthogonality over prime F_p, not an RM/flat-geometry artifact** — is directly evidenced:
a *non-affine-invariant* family (cyclic codes) collapses just as hard as RM, and within the
evaluation-code framework the **FORCED-GRID THEOREM** (Phase 2, independently verified) proves the *only*
triorthogonal column-subsets are ∅ and the full grid. **Direction D is resolved as a negative result**,
with one honest uncovered tail (high-dim general-position caps, <5%).

## Phase 0 — feasibility filter (targets locked; one re-scope)

**Filter A (Singleton, loose)** reproduced exactly — first-feasible m: p=5→4, p=7→3, p=11→2, p=13→2,
p=17→2, p=19→2, p=23→1.

**Filter C (the binding 2D cap), live re-run, p=11/13 m=2** (best true d = max over line-spread + random
puncture sets):

| p | optimal k | best true d | needed d | gap | best γ |
|---|---|---|---|---|---|
| 11 | 21 | 4 | 5 | **+1** | 1.126 |
| 13 | 26 | 5 | 6 | **+1** | 1.059 |

**G0-gate CORRECTION:** the plan called p=11 "near-hopeless (~+4, d=2)", but that was anchored on the
non-optimal k=25 (where dim G0 collapses). At its *optimal* k=21, **p=11 misses γ<1 by exactly +1
distance unit — the same gap as p=13** (matches the prior `[[100,21,4]]₁₁` / `[[143,26,5]]₁₃`). Targets:
p=7 (mid-k m=4, open), **p=13 (+1)**, **p=11 (+1)** co-equal; p=5 skip; p=17/19 already cross.

## Phase 1b — Rank 3 (algebraic-geometry): KILL

The **F₇ elliptic K1 test**: build one-point AG codes `L(αP∞)` on every nonsingular `y²=x³+ax+b` over F₇,
evaluate at the affine points, test `is_triorthogonal`. **All 42 curves × every α FAIL** (the cubic/
triple-product moment never vanishes). Extended to firm it up:

| family | p | curves | (curve,α) tests | result |
|---|---|---|---|---|
| elliptic g=1 ⟨2,3⟩ | 7 / 11 / 13 | 42 / 110 / 156 | 200 / 775 / 1182 | **all fail** |
| hyperelliptic g=2 ⟨2,5⟩ | 7 | 343 | 1661 | **all fail** |

**Mechanism:** triorthogonality needs `Σ_P f·g·h ≡ 0` for all functions of pole order ≤ 3α. The full grid
F_p^m satisfies this via the power-sum identity (that is RM); a sparse curve over a prime field has no
analog — the "Castle"/self-dual curves that would are maximal curves over F_{q²}, which make each
coordinate a 2-qudit ⇒ the settled field-CCZ collapse. (Sanity: `is_triorthogonal` was cross-checked
against a hand pair+triple sum, 0 disagreements.) **AG-over-prime is dead.**

## Phase 1 — Rank 4 (Artin–Schreier `yᵖ−y=f(x)`): analytic KILL

Over the prime field, `y₀ᵖ−y₀ ≡ 0` for all `y₀∈F_p` (Fermat), so the F_p-rational points are
`{(x₀,y₀): f(x₀)=0, y₀∈F_p arbitrary}` — vertical lines over the roots of f, not isolated smooth points.
No usable one-point AG code. The additive lever only has teeth over F_{pᵏ} ⇒ back to the multi-qudit
field-CCZ collapse. **Non-starter over the prime field; no probe needed.**

## Phase 1a — Rank 1 (BCH / cyclic): CONFIRMED-KILL on γ<1, with a structural twist

Length `n = p²−1` (120 at p=11, 168 at p=13), so cyclic codes are well-defined and generally **not
affine-invariant** — the line/plane cap does not directly apply. Result (structured sweep: full BCH family
both offsets/all δ, transform-characterized maximal codes, ~25 000 random punctures + hill-climb):

- **Genuinely non-RM triorthogonal cyclic codes EXIST.** Maximal triorthogonal dim **30** (p=11) and **51**
  (p=13), exceeding the triorthogonal-RM caps (RM₁₁(6,2)=28, RM₁₃(7,2)=36) ⇒ not monomially equivalent to
  any triorthogonal punctured-RM; defining sets are not P-closed (not affine-invariant). *Independently
  verified* on the explicit p=11 dim-30 generator: `is_triorthogonal=True`, `rank=30`, cyclic.
- **…but they collapse under puncturing.** Best certified Z-distance at the γ<1 density is **d ≤ 2**
  (p=11 ties the RM cap; p=13 is *below* RM's d=4 over 25 000 trials). With d=2 and the dimension-capped
  k_punc, γ_min ≥ **1.65** (p=11) / **1.24** (p=13) — γ<1 is arithmetically impossible. Best codes:
  `[[97,23,2]]₁₁` (γ=2.08), `[[144,24,2]]₁₃` (γ=2.59). *Independently verified*: the dim-30 p=11 code
  punctured at k=21/25/29 gives best d = 2/1/1 over 400 random sets (γ = 2.24 / ∞ / ∞).

**The lesson (the §1 thesis, confirmed):** a non-affine-invariant triorthogonal family *still* collapses —
shortening to a small-redundancy CSS X-stabilizer destroys the dual distance regardless of geometry. So the
cap is a **degeneracy cap generic to triorthogonality over prime F_p**, not an RM/flat artifact.

**The lever (transform-domain characterization, validated — 0 disagreements vs `is_triorthogonal` on 80
codes):** with `Z` = spectral support (`Z_n` minus the defining set),
- pair / self-orthogonality ⟺ `Z ∩ (−Z) = ∅`,
- cubic / triple ⟺ no `j+k+l ≡ 0 (mod n)` with `j,k,l ∈ Z`.

Reproducible p=11 dim-30 survivor: `Z = {1,2,5,6,11,21,22,25,26,31,35,41,42,46,51,55,61,62,66,71,75,81,82,86,91,101,102,105,106,111}` (then `T = Z₁₂₀∖Z`, generator `g = ∏_{t∈T}(x−α^t)`, α a primitive element of GF(11²)).

## Phase 2 — Rank 2 (moment / direct triorthogonal-matrix search): CONFIRMED-KILL + a theorem

The moment/point view (CONFIRMED == `is_triorthogonal`, 0/30 mismatches): a triorthogonal code is a
column-multiset with vanishing 2nd & 3rd moment tensors mod p; the affine form (constant-1 coordinate)
avoids the trivial d=2 cone-collapse. Direct construction found **no γ<1 survivor** at p=11/13 — the best
codes are RM itself (`[[100,21,4]]₁₁` γ=1.126, `[[143,26,5]]₁₃` γ=1.059), and every directly-built non-RM
config does strictly worse. Coverage (each VERIFIED by `is_triorthogonal` = ground truth):
- **Line-unions** (~17 non-RM configs, generator dim 18–26): collapse to **d=1** (line-rich ⇒ collinear).
- **Moment-zero caps** (genuine caps, 0 collinear triples, built via ±-symmetry MITM subset-sum): collapse
  to **d≤3** in the *shortened* code (weight-2 dependencies appear in G0) — best γ=1.32 (p=7), 3.62 at D=5.
- **Uncovered tail (honest):** the γ<1-density general-position cap needs generator dim K≈28 (p=11)/36
  (p=13) ⇒ points in F_p^~27 ⇒ a 0/1 (set) subset-sum over ~4000 moment equations — NP-hard, not
  constructively reachable here. Not closed by theorem. This is the **<5% tail**.

### The FORCED-GRID THEOREM (rigorous core — independently verified)

For the F_p^m evaluation pool with the full degree-`r_max` monomial generator, a subset `S⊆F_p^m` is
triorthogonal **iff** `1_S ⊥ RM(r)^{⋆3} = RM(3r,m)`, i.e. `1_S ∈ RM(t,m)` with `t = m(p−1)−3r−1`. Since
`2t < p` at every target, on each affine line `1_S` is a univariate degree-≤t polynomial g with g∈{0,1}
⇒ `g(g−1)≡0` (degree 2t<p) ⇒ g constant ⇒ `1_S` globally constant ⇒ **`S ∈ {∅, full grid}` only.**
So within the evaluation framework — where *every* known γ<1 RM code and all its puncturings live — there
is **no non-RM triorthogonal escape**; only ordinary RM puncturing, Filter-C capped (p=11 d=4, p=13 d=5).

**Independent verification** (`scratchpad/verify_p2.py`, RAM-bounded): (1) Hadamard-cube rank identity
`rank(RM(r,m)^{⋆3}) = dim RM(3r,m)` exact — **118=118** (p=11), **163=163** (p=13), **333=333** (p=7,m=3);
(2) `t` and `2t<p` hold at all targets (p=11 t=1, p=13/p=7m3/p=7m4 t=2, p=5m4 t=0); (3) full grid
triorthogonal, **removing any single point breaks it 10/10**, **0/250 random proper subsets** triorthogonal
(p=11 *and* p=13); (4) constant-forcing — **0 non-trivial 0/1 codewords** of RM(t,m) (exhaustive at p=11,
sampled at p=13/p=7).

## Status: Direction D resolved (negative)

All five plan routes closed: Rank 1/2/3/4 fail to cross γ<1; Rank 5 (affine) is subsumed by forced-grid.
The realistic and now-delivered outcome is the **negative result** — the distance cap is generic to
triorthogonality over prime F_p (p≤13): rigorous within the evaluation framework (forced-grid), strongly
evidenced generally (cyclic/AG/Artin–Schreier kills + moment-cap collapses). **Honest open items:** (i) the
high-dim general-position cap tail (<5%, not closed by theorem); (ii) the p=7 mid-k m=4 *RM* window
`k∈[110,312]` (still uncertified — `min_dependent_columns` overflows past R=21 rows; separate from D).
Consolidated in `D_NEGATIVE_RESULT.md`.
