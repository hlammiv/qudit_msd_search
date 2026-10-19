# Phase 2 — γ<1 with a clean diagonal carrier: unreached; the "d≤2" obstruction is OPEN, not a theorem

*Status 2026-06-30.  Workflow (prove ∥ refute → critique → synthesize) + PI re-verification.  The user's
goal: a low-overhead (γ<1) genuine ring-qudit MSD code with a DIAGONAL ω=i carrier (clean injection).*

## Headline (three-valued, honest)

- **γ<1 with a clean diagonal carrier was NOT achieved.**  Best genuine clean-diagonal-injection code is
  still **[[8,2,2]]₄, k=2, d=2, γ=2** (PI-verified: hosts the diagonal ω=i carrier, induces a genuine
  non-Clifford level-3 logical on (Z₄)², anti-collapse passes).
- **"Diagonal-carrier hosting ⇒ d≤2" is NOT a closed-form theorem.**  The prove-track's distance lemma
  (hosting ⇒ rank(M_X)≤1, |M_X|≤4 ⇒ d_Z≤2) was only computational (n≤7) and **never searched the all-even
  M_X = 2U family** — exactly where the workflow's adversarial critique reported an all-even (two-Steane,
  n=14) code that hosts the carrier with |M_X|=64 and d=3.  *PI caveat:* I could **not** independently
  reproduce that hoster in a quick rebuild (the construction needs a specific Z-part Y with U⊥Y; my naive
  Y=U gave k=8 and hosts=False), so I treat the n=14 refuter as **workflow-internal, not PI-confirmed** —
  but the lemma being gapped (computational n≤7, all-even unsearched) stands independently.
- **For GENUINE (Z₄)^k level-3 codes, γ<1 is OPEN.**  Across thousands of genuine diagonal-carrier hosters
  (n≤13), **none reached d≥3** — every genuine hoster caps at d≤2.  The only known escapes to d≥3 are
  **all-even M_X = 2U codes, which are NON-genuine** (logical X has order 2 ⇒ (Z₂)^… content ⇒ the carrier
  induces a logical *Clifford*, not level-3 magic — and they're high-overhead anyway, e.g. n=14 γ=2.40).

**Net:** a genuine clean-injection ring MSD code exists (γ=2); γ<1 is **strongly obstructed but formally
open** — neither reached nor proven impossible with a diagonal carrier.

## What survived adversarial critique (the durable algebra)

- **Hosting criterion (exact):** with Δ_a q(x)=q(a)+B(x,a), B(u,v)=2(u mod 2)(v mod 2) mod 8 — the diagonal
  carrier hosts iff Q(a)≡0 mod 16 **and** the cocycle-orthogonality Σᵢ B(xᵢ,aᵢ)≡0 mod 16 on supp(C), for
  every X-stabilizer a.  Cross-validated both directions vs the dense projector.
- **Full support ⇒ no hosting** (the ω=i θ is non-constant); hosting ⇒ confined support ⇒ pure-Z stabilizers.
- **Exact support identity |supp(C)| = 4^k·|M_X|** (module duality) — held on every code incl. the refuter.
- **Shadow lemma L1:** Δ_a q mod 2 = parity(a) (x-independent) ⇒ hosting forces an even number of odd
  entries per X-stabilizer.
- **Genuine ⇒ odd-X**, but ONLY with the anti-collapse/odd-commutator check (the bare non-quadraticity
  helper gives even-X false positives — the genuineness oracle is incomplete alone).

## What broke / is gapped
- **REFUTED:** the rank cap |M_X|≤4 and the d_Z≤2 Singleton consequence (all-even M_X=2U escapes).
- **Gapped:** the one-line 2-torsion "weight-2 pure-Z logical" argument (needs equal parity across *all*
  a∈M_X, not one); the minwt(H^⊥\W_z)≤2 reduction is strong but computational and false without a
  genuineness restriction.

## Best codes per family (verified status noted)

| family | code | [[n,k,d]] | γ | genuine level-3? |
|---|---|---|---|---|
| non-CSS witness | [[8,2,2]]₄ | [[8,2,2]] | **2.0** | **YES (PI-verified)** |
| CSS over Z₄ | [[7,2,2]]₄ | [[7,2,2]] | 1.81 | hosts by M1; **genuineness NOT verified — do not cite as genuine** |
| all-even 2U (two Steane) | n=14 | [[14,1,3]] | 2.40 | **NO** (logical (Z₂)², carrier→Clifford); **PI could not reproduce hosting** |
| BH/RM triorthogonal lift | — | — | — | fails (qubit triorthogonality ↛ Z₄ Q≡0 mod16) |
| HGP / concat / ⊕ | — | — | — | HGP not isotropic over Z₄; concat → entangling logical; ⊕ keeps γ=2 |
| genuine sampler (n≤13) | thousands | all d≤2 | ≥1 | yes, but ZERO reached d≥3 |

## The single open gap + the recommended next step

**The live question:** does a GENUINE (odd-X) diagonal-carrier hoster with |M_X|>4 and d≥3 exist?  The only
route to d≥3 (all-even 2U) destroys genuineness; whether odd-X can coexist with a high-dual-distance X-code
at k≥1 is unresolved.  **Best next step:** template off the all-even refuter but force odd entries into M_X
(keeping Q≡0 mod 16 + B-orthogonality on support), and test whether genuineness (odd logical commutator)
survives |M_X|>4 at d≥3.  Either it yields a genuine d≥3 hoster (push to γ<1) or it provably cannot (odd-X
⇒ |M_X|≤4 or k→0), which finally closes the diagonal route as a theorem.  **Run the non-diagonal
(Clifford-conjugated) carrier in parallel as the hedge** — it is untouched by every Phase-2 obstruction, and
(being semi-Clifford) it keeps injection clean *iff* its induced logical is diagonal.

## Honest flags (per the synthesis + PI)
- The rate-route "Bravyi-Haah floor n/k≥4 ⇒ γ<1 closed at d=2" was **not adversarially vetted** — treat as
  un-confirmed (though qubit BH does floor at γ→log₂3≈1.58, so γ<1 at d=2 is implausible).
- [[7,2,2]]₄ γ=1.81 is **not** verified genuine — do not cite.
- The n=14 refuter is workflow-internal; **PI verified only the framework (M1, hosting criterion, support
  identity) and the genuine codes [[8,2,2]]₄ / [[5,1,3]]₄ / [[12,2,3]]₄**, not the n=14 hoster.
