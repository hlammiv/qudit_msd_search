# Direction D — Phase 0 & 1 execution findings

2026-06-30. Executes `D_PLAN.md` Phases 0–1 (feasibility filter + the two cheapest decisive probes).
All local on the 15 GB box; distance certified by `min_dependent_columns(d_max≤5)` (RAM-safe).

## Headline

**Three of the plan's escape routes are now CLOSED** — Rank 1 (cyclic/BCH) on the γ<1 objective,
Rank 3 (algebraic-geometry), and Rank 4 (Artin–Schreier). The central thesis (`D_PLAN.md` §1, §6)
— **the distance cap is generic to triorthogonality over prime F_p, not an RM/flat-geometry artifact**
— is now directly evidenced by a *non-affine-invariant* family (cyclic codes) that collapses just as
hard as RM. The **publishable-negative outcome is the front-runner**. One genuine search route remains
(Rank 2, moment-ILP) plus the affine falsifier (Rank 5) and the still-uncertified p=7 mid-k RM window.

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

## What remains

- **Rank 2 — moment-ILP direct search (Phase 2).** The last genuine search route. The triorthogonality
  moments are *linear* in column multiplicities ⇒ mod-p linear feasibility on a small pool + capset
  hill-climb. The §1 thesis predicts the same collapse, but this is the cleanest test of the *general*
  claim (and the <5% tail where a survivor could appear). **Recommended next phase.**
- **Rank 5 — PRM / cap-puncture.** Fast affine-route falsifier; expected to fail (inherits the flat cap).
- **p=7 mid-k m=4 RM window** `k∈[110,312]`. Still uncertified (`min_dependent_columns` overflows past
  R=21 rows). This is *RM*, not non-RM — closing it rigorously (overflow-safe distance) would tidy the RM
  no-go but is a separate task from Direction D.

**Assessment.** With Rank 1 (γ<1), Rank 3, and Rank 4 closed and the degeneracy-cap thesis directly
evidenced, the realistic deliverable is the **negative result**: the cap is generic to triorthogonality
over prime F_p (p≤13). Rank 2 is the one experiment that would make that airtight — or, against the odds,
produce the survivor. Per-p success odds revised down from the plan: p=7 mid-k (RM) is the only genuinely
open window; non-RM p=11/p=13 now look as caped as everything else.
