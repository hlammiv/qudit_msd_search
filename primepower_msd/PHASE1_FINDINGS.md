# Phase 1 — the diagonal-carrier route hits a clean structural obstruction (with an M1 theorem)

*Status 2026-06-29.  Workflow (M1 derive ∥ M2 GR-construct → decisive experiment → verify → synthesize)
+ independent PI re-verification of every load-bearing claim the workflow flagged as not adversarially
re-derived.*

## Headline

**No γ<1 code, no distillable ring-magic code — but a sharp, mechanistically-explained OBSTRUCTION for
the DIAGONAL-carrier route, a solid M1 theorem, and two verified boundary witnesses that bracket it.**

The make-or-break (distance) is **not** the killer.  The killer is a squeeze between *hosting the carrier*
and *having distance*:

> **For a diagonal ω=i carrier over Z₄: hosting (M1) ⟺ pure-Z stabilizers (confined support) ⟺ d ≤ 2.
> Distance d ≥ 3 ⟺ no pure-Z stabilizers (full support) ⟺ the carrier is NOT hosted.  Mutually exclusive.**

## M1 — the ring-symplectic triorthogonality condition (DERIVED, adversarially confirmed)

For a diagonal level-3 carrier U=diag(ζ₁₆^{q(x)}) applied as U^{⊗n}, the codespace of a Z₄ stabilizer
code S is preserved **iff** for every a ∈ proj_X(S) the diagonal Clifford D_a=diag(ζ₁₆^{−Σᵢ Δ_{aᵢ}q(xᵢ)})
fixes the codespace, where Δ_a q(x)=q(x+a)−q(x).  Z-only stabilizers (a=0) give D₀=I.  The non-additivity
D_{a+a'}=D_a D_{a'} diag(ζ₁₆^{−Δ_aΔ_{a'}q}) (the *second* difference = CZ/symplectic coupling) is why the
condition is over the whole module proj_X(S), not just generators.  This is the genuine ring-level-up
analog of Bravyi–Haah CSS triorthogonality; the CSS case collapses to coset-constancy of the carrier
phase on M_X (the package's earlier-validated test).  **All bookkeeping is mod 8 (carrier mod 16)** — the
discriminating invariant Δ³q ∈ {4→ω=i, 8→ω=−1 soft, 12→ω=−i} is *invisible mod 4*.

**Verified:** foundation identity to 1.47e-15; reproduces the [[8,1,2]]₄ gate and the [[5,1,3]]₄ V7 null;
**adversarial re-derivation 0 disagreements / 312 trials** (0/48 null, 48/48 gate, 0/216 random), distinct
code paths.  **Adversarial verdict: PASS, strengthened.**  Scope: M1 assumes a *diagonal* U.

## M2 — GR(4,r) distance (the make-or-break): d≥3 at k≥2 IS achievable

A sharp dichotomy on the self-adjoint twist M:
- **Naive multiplicative/diagonal twist** (literal "ev Mf" RM): **d ≤ 2 for all k** — structural, not search-
  limited (the a-side code rank is capped by GR multiplicative structure, so its dual carries weight-≤2
  codewords ⇒ pure-Z logicals; larger r is monotonically worse). 0/4000 with d≥3.
- **Off-diagonal (coordinate-coupling) twist** breaks the d=2 cap: d≥3 at k=1..5.

**Best (PI-verified):** `[[12,2,3]]₄` non-CSS, isotropic, |S|=4¹⁰, k=2, distance 3 (brute symplectic
enumeration), **0 pure-Z stabilizers ⇒ full support**, γ = log(6)/log(3) = **1.63**.  Generators saved at
`M2_code_12_2_3.npy` (repo root).  So **the GR distance gate is NOT the binding constraint** — but γ<1
needs k≥3, N≤8 (e.g. [[8,3,3]]₄ γ=0.89) and random search found none (0/30000); existence open.

## The decisive experiment — OBSTRUCTION (both witnesses PI-verified)

The sharpened M1 "support criterion" (feasible at any n): U^{⊗n} preserves C ⇒ supp(C) ⊆ Z(a) =
{x : Σᵢ Δ_{aᵢ}q(xᵢ) ≡ 0 mod 16} for every a∈proj_X(S), where supp(C)=Wz^⊥ (pure-Z-stabilizer-orthogonal).
- ω=i carrier ⇒ θ(·,a) strictly non-constant (range hits 15) ⇒ Z(a) a proper subset ⇒ **hosting needs
  confined support ⇒ nontrivial pure-Z stabilizers**.
- d≥3 ⇒ maximally non-CSS ⇒ **no pure-Z stabilizers ⇒ full support (density 1.0)**.

Verified on both genuine d=3 codes (0/12 carriers hosted on [[5,1,3]]₄ and [[12,2,3]]₄), and a loophole-
closing hosting search at n=6..10.

**Two boundary witnesses (both independently PI-verified):**
- **`[[8,2,2]]₄` (closest to success):** non-CSS, k=2, **genuinely HOSTS the ω=i carrier** (M1 condition
  holds on all support points for c=1,2,3 — equivalently dense leakage 3.28e-15).  A genuine multi-ring-
  qudit code transversally carrying the genuine ring-CCZ carrier — **but d=2 (detection only, γ=2).**  It
  has 1 X-generator a₀=[1,1,1,1,1,1,3,3] and **6 pure-Z stabilizers** (the confined support that lets it
  host — and caps its distance).
- **`[[12,2,3]]₄`:** d=3, k=2, but **no pure-Z stabilizers ⇒ hosts NO carrier**, γ=1.63.

## Honest scope (workflow flags + PI checks)
- M1 is adversarially hardened.  The M2 codes, the two witnesses, and the obstruction were re-verified
  **by the PI this turn** (isotropy, k via dⁿ/|S|, distance by brute symplectic enumeration, hosting via
  the verified M1 condition) — they hold.  The obstruction is **computational (n≤10 + the two explicit
  codes + the mechanism), not a closed-form all-n proof.**
- γ<1 was never reached; a near-MDS [[≤8,≥3,3]]₄ code's existence is **open, not refuted.**

## The single best next step

**Abandon the diagonal carrier; extend M1 to a NON-DIAGONAL (Clifford-conjugated) level-3 carrier
C₁·D·C₂.**  The entire obstruction rests on one fact — a *diagonal* U gives D_a depending on the first
difference of q, whose non-constancy forces confined support incompatible with d≥3.  A non-diagonal
carrier breaks the "supp(C) ⊆ Z(a)" derivation at its root; M1 itself flags diagonal-only as its scope.
Concretely: re-derive the support/triorthogonality condition for C·diag·C† carriers and re-test against
the verified [[12,2,3]]₄ d=3 code that the diagonal carrier cannot host.  (Secondary: multivariate
GR(4,r)^m RM, the one surviving distance direction — but distance is not the binding constraint, so this
is subordinate to the carrier question.)


> **2026-10-19 update (UNIVERSALITY_CHECK.md):** the `[[8,2,2]]₄` carrier is hosted, but its INDUCED LOGICAL is a mod-4 phase (no A₂ escape, no odd third difference) — the logical gate is 2-torsion ("soft"), not universality-enabling. Transmission needs a logical-X axis with an odd number of odd entries (`rgt.axis_parity_report`).
