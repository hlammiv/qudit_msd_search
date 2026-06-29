# M3 — THE RESULT: single-qudit prime-power (d=4) magic state distillation EXISTS

*Status 2026-06-29.  Code: `m3_search.py`, scalable kernel in `ringlinalg.right_kernel_howell`.
Verified witness frozen in `tests/test_m3_witness.py` (42 tests pass).*

## Headline

> **A single-ququart (d = 4) magic-state-distillation code exists.**
> `[[9, 1, ≥2]]₄` CSS code over the **cyclic clock ring Z₄**, with a **transversal genuinely-quadratic
> level-3 gate** (antidiff(S³Z¹)) that induces the **strict level-3 logical gate diag(1, −i, 1, i)**.
> Verified algebraically AND by dense codespace simulation in C^(4⁹) — **leakage = 1.3×10⁻¹⁵**.

This affirmatively answers the project's central question.  The earlier nulls (M2, M2b, R1 — all at
n ≤ 7) were a **finite-size artifact**: those searches enumerate the d^n grid, so they could not reach
the block lengths where witnesses live.  M3 removes the d^n ceiling and the obstruction dissolves.

## The verified witness (`[[9,1,≥2]]₄`)

- **Physical transversal gate:** antidiff(S³Z¹) = diag(ζ₁₆^[0,1,12,13]); **stab1 = {0}** ⇒ genuinely
  quadratic (the HARD case — not a qubit-like gate).
- **X-stabilizer:** M_X = ⟨(1,3,3,0,1,1,1,1,1)⟩ over Z₄ (with Σᵢ φ(aᵢ) ≡ 0 — the self condition).
- **Logical X generator:** g = (1,0,0,0,1,1,2,0,3), order 4 in M_Z^⊥/M_X ⇒ **cyclic Z₄ logical**.
- **k = 1**, |M_X|·|M_Z| = 4⁸ = d^{n−1}; X/Z commutation holds.
- **distance ≥ 2:** no weight-1 logical on either side.
- **Induced logical gate:** diag(1, −i, 1, i), oracle level = 3, not Clifford ⇒ **strict level-3**.
  (Independently: it is not in ⟨S,Z⟩ so not Clifford; UXU† = X·S²Z ∈ Clifford so it is C₃.)
- **Anti-collapse:** logical clock is cyclic Z₄ (shift order 4) ⇒ genuine single-ququart, not two qubits.
- **Dense simulation in C^(4⁹=262144):** D^{⊗9} preserves the codespace (leakage 1.3e-15), the induced
  logical matrix is exactly diagonal (off-diag norm 0) and re-certifies as strict level-3.

A second witness exists at n=13 (transversal antidiff(S¹Z³) ⇒ logical diag(1,i,1,−i)); witnesses are
**reproducible** (per-n counts {9, 13} across a sweep) — not a fluke.

## Why M3 reached it (the scalable substrate)

Nothing in M3 touches the d^n grid:

- **Scalable Howell right-kernel** `right_kernel_howell` (augmented-Howell `[Aᵀ | I]`), validated
  exactly against the brute kernel — gives M_Z = L^⊥ and M_X^⊥ as **generators**, no enumeration.
- **Transversality is checked only on the small module L** = M_X + ⟨g⟩ (|L| ≤ d·|M_X| ≤ 64): the
  CSS codewords are |M_X|-sparse, so "Φ constant on each M_X-coset of L" **is** the exact codespace
  condition (the dense C^(4⁹) sim confirms the equivalence at leakage 1e-15).
- **distance ≥ 2 == "no weight-1 logical"** — checked over the n·(d−1) weight-1 vectors via membership.

So n is limited only by O(n) bookkeeping; n = 9, 11, 13, 15 are all trivial (the full sweep runs in
~7 s on 14 cores).

## Honest scope — what is and isn't established

- **Established:** existence of a genuine single-ququart MSD code (transversal level-3 logical gate on
  a distance-≥2, k=1, cyclic-Z₄ CSS code), fully verified.  The hard genuinely-quadratic gate works.
- **NOT yet:** **low overhead.**  This minimal witness has distance ≈2 and γ = log(n/k)/log(d_code)
  ≈ 3 — high (the prime [[15,1,3]]₂ is 2.47).  Existence ≠ efficiency.
- Exact code distance and A_d are not yet computed (only distance ≥ 2 certified); d = 4 only so far.

## Next

1. **Optimize overhead** — now that existence is settled, search the (scalable) construction for higher
   distance / larger k / smaller n-per-distance to drive γ down (toward, and below, the prime baselines).
2. **Exact distance + A_d over the ring** (finish M3): a ring min-weight / MacWilliams engine for the
   real suppression δ_out ~ A_d δ_in^d and an honest γ.
3. **The distillation routine** — wire the code into `ring_distill` (yield/cost), confirm δ_out < δ_in.
4. **d = 8, 16, 32** — repeat with the parametric substrate; expect witnesses to appear at larger n.
5. **Characterize the witness family** — the self+transversal+level-3 conditions define a structured
   set; a closed-form (analytic) sub-family would replace search and likely give better parameters.
