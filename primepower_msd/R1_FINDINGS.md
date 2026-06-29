# R1 findings — analytic ring-triorthogonality construction (and the consolidated d=4 obstruction)

*Status 2026-06-29. Code: `r1_search.py` (+ `ring_css.build_css_general`, `ring_transversal.stab1`).
Parallel + resource-capped (`parallel.py`).*

## What R1 does (construct, don't sample)

M2b left one honest gap: the genuinely-quadratic level-3 gates have trivial *global* translation
stabilizer, but the true per-coset transversality is weaker and could hold on a special code.  That
"special code" condition is **ring-triorthogonality**:

> pick X-stabilizers M_X and a logical-rep module L = M_Z^perp with M_X ⊆ L such that
> **(self)** Σᵢ φ(aᵢ) ≡ 0 (mod N) ∀ a ∈ M_X, and **(transversal)** Φ(c+a) ≡ Φ(c) ∀ a ∈ M_X, c ∈ L.

The transversality region **V = { c : Φ(c+a) ≡ Φ(c) ∀ a ∈ M_X }** is a *union of M_X-cosets*, so R1
**constructs** transversal codes instead of sampling: find a cyclic order-d coset generator g with
g, 2g, …, (d−1)g all in V, set L = M_X + ⟨g⟩, M_Z = L^perp.  This reaches the measure-zero solutions
the random searches (M2/M2b) cannot.

## R1 result (d=4, n ≤ 7)

- **19,691 self-orthogonal M_X** and **209,012 constructed per-coset transversal codes** in 567 s / 14 cores.
- **Quadratic gates CAN be per-coset transversal** — the construction builds them in bulk (this
  *resolves M2b's open question affirmatively*: the trivial-global-Stab gates are still transversal
  on special codes).
- **But every constructed code is distance 1.**  **No distance-≥2 ring-triorthogonal code for any of
  the 12 level-3 gates.**

## The consolidated obstruction (three independent methods agree)

| method | family | scale | distance-≥2 single-ququart MSD code? |
|---|---|---|---|
| **M2** | weakly-self-dual, random | 443 codes (n≤7) | **none** |
| **M2b** | general CSS, Stab(φ)-guided | 537,600 codes (n≤7) | **none** |
| **R1** | analytic ring-triorthogonal construction | 209,012 constructed transversal codes (n≤7) | **none** |

Across all three, **transversal level-3 single-ququart gates exist only on distance-1 (unencoded)
codes.**  The mechanism R1 exposes: the transversality region V is structurally *thin* — forcing the
logical operators (which must lie in L ⊆ V) to remain weight-1.  Transversal-level-3 and error
protection are in direct tension for **diagonal-transversal CSS codes over the cyclic ring Z_4**.

## Honest scope — this is strong evidence, not yet a theorem

Bounded to **d = 4, n ≤ 7, diagonal gates, CSS, k = 1**.  Two escape routes remain genuinely open:

1. **Larger n.** In M2, distance-≥2 codes *first appear* at n=7; transversal-and-protected codes may
   need n ≥ 9–11.  The brute d^n region/dual enumeration caps n (RAM), so this needs the **M3
   Howell-scaled ring kernel/distance** (replace the d^n grid with a scalable kernel).
2. **d = 8 (and 16, 32).** The substrate is already parametric; the d=8 Hoggar-lines structure
   (Zhu) behaves specially and is worth a separate R1 pass.

## Recommended next

- **M3**: build the Howell-scaled ring distance / kernel to push R1 to n ≥ 9–11 (the main caveat).
- **d = 8 R1**: re-run the construction at d=8 (parametric substrate already supports it).
- **No-go proof**: formalize *why* V is thin — prove that for a level-3 diagonal gate over cyclic Z_4,
  M_X ⊆ L ⊆ V with k=1 forces a weight-1 logical.  A clean theorem here is a publishable negative and
  would explain all three empirical nulls.
