# Phase 3 — the generalized-M1 (semi-Clifford) test CLOSES Track B on the best d=3 testbed

*Status 2026-06-30.  PI-built generalized-M1 normalizer test, cross-validated 192/192 against the
dense ground truth, then run exactly on [[12,2,3]]₄.*

## Headline

**`[[12,2,3]]₄` — the strongest verified d=3 ring-qudit code — hosts NO semi-Clifford ω=i carrier
`C₁ D C₂` in ANY of the 192 single-qudit Clifford frames (48 symplectic × 4 carrier shifts).
Exact, polynomial, no 4¹² object.**  This was the one decisive experiment left open at the end of
Phase 2 (the tandem A+B caveat: "[[12,2,3]] untested under the full carrier family, n=12
dense-infeasible").  It is now settled: **the diagonal/semi-Clifford route to γ<1 via [[12,2,3]] is
NULL, as a theorem on this code — not a computational gap.**

## The reduction that makes it polynomial (exact)

A semi-Clifford carrier is `W = C₁ D C₂` with `D` the diagonal ω=i carrier `diag(ζ₁₆^q)`,
`q=[0,1,4,13]`.  Key fact: **a diagonal level-3 `D` maps a stabilizer code to a *Pauli* stabilizer
code iff it M1-hosts it** (otherwise the conjugated stabilizers carry non-trivial diagonal-Clifford
corrections, so the image is not a Pauli code).  Therefore `D^⊗n (C₂·C) = C₁†·C` forces
`C₂·C = C₁†·C`, and the full family collapses to:

> `C` hosts some `C₁ D C₂`  ⟺  ∃ single-qudit Clifford `c` with `D` hosting `c·C`
> ⟺  the Clifford-conjugated carrier `W = c† D c` hosts `C`.

Sweeping `c` over the single-qudit Clifford group = **48 symplectic frames `SL(2,Z₄)`**
(the X↔Z-mixing part) **× 4 carrier shifts `q(x)→q(x+s)`** (the Pauli-translation part of `C₂`;
the Z-translation commutes with `D`).  For each `(c,s)`: transform every stabilizer generator's
symplectic vector `(aᵢ|bᵢ)` by `M_c`, then run the diagonal M1 hosting criterion with the s-shifted
carrier.  All operations are on the 10×24 generator matrix — polynomial in n.

## Why every frame fails (the sharp mechanism)

`0 ∈ supp(C)` always (`0` is in every right-kernel), so hosting **requires**
`θ(a,0) = Σᵢ Q[aᵢ] ≡ 0 mod 16` for *every* `a ∈ proj_X(S)`.  In all 192 combos `[[12,2,3]]₄` has at
least one X-stabilizer generator with `Σ Q[aᵢ] ≢ 0 mod 16` (e.g. identity frame Q-weights =
`[13,14,1,9,3,11,13,12,3,12]`).  **No Clifford rotation aligns all of [[12,2,3]]'s stabilizer
Q-weights to 0 mod 16.**  Reject-path tally: 192/192 fail at this x=0 condition — the full
support/cocycle machinery is never even needed.

## Validation (airtight)

- **Convention + full criterion** cross-validated on `[[5,1,3]]₄` (n=5, dense feasible): the
  polynomial pipeline `[transform by M_c → shifted M1 hosting]` vs the dense ground truth
  `transversal_test(c†D_s c)` — **192/192 agree, 0 disagreements.**  (Side-result: `[[5,1,3]]₄`
  hosts no semi-Clifford ω=i carrier either — the clean 0/192 closure, superseding the earlier
  0/9216 dense sweep.)
- The n=12 NOs come via the **sound** cheap filter (`θ(a,0)≠0 ⟹ not hosted`), which is a logical
  subset of the dense-validated full criterion.  Identity frame reproduces the known Phase-1
  result (`[[12,2,3]]` hosts no plain diagonal carrier).
- Code: `scratchpad/genM1.py` (group + crossval), `scratchpad/n12run.py` (sweep).
  Clifford group closes at 768 = 48 symplectic × 16 translations.

## What this does and does NOT prove

- **DOES:** the best d=3 ring-qudit testbed cannot carry a genuine ω=i magic gate under *any*
  single-qudit semi-Clifford dressing.  Both Phase-2 routes (diagonal + non-diagonal/semi-Clifford)
  are now closed *on this code*, exactly.
- **Does NOT:** prove γ<1 impossible in general.  It closes one code (the strongest we have).
  The x=0 condition `Σ Q[aᵢ] ≡ 0 mod 16 ∀ a∈proj_X(S)` is now a **cheap, polynomial, frame-by-frame
  necessary screen applicable to ANY candidate code** — the natural engine for a broad d≥3 sweep
  (the recommended next step) toward either a hosting d≥3 code (→ γ<1) or a general theorem.

## Population screen (the broad sweep)

The x=0 kill condition (`Σᵢ Q[aᵢ] ≡ 0 mod 16 ∀ a∈proj_X(S)`, checkable per-frame in ms) was pointed
at a **population of ~170 genuine d≥3 ring codes** (`scratchpad/pop_screen*.py`), RAM-capped at 4 GB:
- 26 random isotropic d≥3 codes (n=6,7,8; k=1,2);
- 145 codes spanning the **full local-Clifford orbits** of `[[12,2,3]]` (k=2) and random d≥3 seeds
  (independent SL(2,Z₄) per qudit — preserves d/k/genuineness but is inequivalent w.r.t. *transversal*
  hosting, so each orbit samples a whole local-Clifford class).

**Result: 0 host; not one code passes even the x=0 screen in ANY of 192 frames.**  Two independent
code families, entire orbits, k=1 and k=2 — uniformly killed by the simplest necessary condition.

### The emerging lemma (toward a theorem)
`Q=[0,1,4,13]` decomposes as: mod 2 → `Q[a]≡a&1`; mod 4 → `Q[a]≡a&1`; mod 8 → `Q[a]≡(a&1)+4(a>>1)`.
So the necessary condition `Σ Q[aᵢ]≡0 mod 16` implies, in the hosting frame, for every X-stabilizer:
(i) **even # of odd entries** (shadow L1, mod 2), and (ii) **#(odd entries) ≡ 0 mod 4** (mod 4).
No d≥3 code tested admits a uniform Clifford frame where the *whole* X-module simultaneously meets
this.  This is the candidate general no-host theorem.

### Honest scope of the screen
- **k≥3 d≥3 codes — the actual γ<1 regime (n<3k) — were NOT tested:** they are hard to construct
  (random/GR sweeps yield ~none at integer k≥3), so the population is k≤2.  The obstruction is
  established broadly for *buildable* d≥3 codes; the γ<1 codes themselves remain unconstructed, and
  the lemma is *conjectured* to extend to them, not verified.
- Result is empirical (≈170 codes), not yet a proof.  Recommended next step: **prove the no-host
  lemma** (even-#odd + #odd≡0 mod4 incompatible with d≥3 over all Clifford frames) rather than
  sample more codes.

## Champion unchanged

Best genuine clean-injection ring MSD code remains **[[8,2,2]]₄, k=2, d=2, γ=2** (PI-verified).
