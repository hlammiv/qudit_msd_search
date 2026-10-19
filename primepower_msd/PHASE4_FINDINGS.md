# Phase 4 — the multi-qudit pivot: γ<1 genuine ring codes EXIST (verified); gate-hosting is the crux

*Status 2026-06-30.  Ideation workflow (7 lenses → 21 ideas → adversarial vetting → synthesis) on the
new relaxation "the magic gate MAY be MULTI-qudit, as long as it is a genuine same-dimension (Z₄) ring
gate (anti-collapse), not a single-qudit carrier."  Plus PI re-verification of the load-bearing claim.*

## Headline (two-part, honest)

1. **VERIFIED MILESTONE — γ<1 genuine d=3 ring-qudit CODES exist.**  Independently certified by PI
   checkers (isotropy, k from two engines, d=3 exact = no weight≤2 logical + weight-3 present, genuine
   anti-collapse = odd logical commutator):
   - **`[[10,4,3]]₄`, γ = log(10/4)/log 3 = 0.834** — `primepower_msd/code_10_4_3_gamma0834.npy` (6 gens).
   - **`[[12,5,3]]₄`, γ = log(12/5)/log 3 = 0.797** — `primepower_msd/code_12_5_3_gamma0797.npy` (7 gens).
   These OVERTURN the earlier program belief "γ<1 needs k≥3, random found 0/30000 (open)".  That was an
   **under-powered-search artifact**: a *targeted GR(4,2) graph-code* search (`gr_gamma.py`, rows=(c|Mc),
   M symmetric over GR(4,2), specific (n,k) target, ~10²–10³ tries) finds them readily
   (`[[12,5,3]]` at trial 49, `[[10,4,3]]` at trial 2524).  `[[8,3,3]]`/`[[8,4,3]]` were NOT found in
   40k tries — both saturate the Z₄ Singleton bound k≤N−2(d−1) and are almost certainly nonexistent
   (the synthesis correctly flagged them mythical).

2. **THE CRUX IS NOW THE GATE, NOT THE CODE.**  A γ<1 *code* is necessary but NOT sufficient for MSD.
   The magic gate must HOST on it: a transversal (or k-local) genuine multi-qudit **ring-CCZ** that
   induces a non-Clifford LOGICAL with an **ODD cubic associator** `m = Σᵢ vᵢ³ mod 4`.  The single-qudit
   program died exactly at hosting (4160 random codes, 0 hits — "hosting must be designed, not found").
   Whether `[[10,4,3]]`/`[[12,5,3]]` host such a gate is **OPEN**.

## The workflow's unifying principle (the durable idea)

**Genuineness lives at the LOGICAL level, not the physical level.**  The multi-qudit ring-CCZ
`diag(i^{Tr(uvw)})` / `diag(i^{xyz})` is structurally genuine *physically* — its order-4 3-body term is
irreducible; no single-ququart Clifford (SL(2,Z₄), the 768-element group) can factor it into a GF(4)/
qubit gate (Clifford preserves clock order + acts linearly).  So the whole program reduces to ONE
tension: **can a code be self-orthogonal enough to HOST the gate, yet keep the induced logical cubic
form `Σ vᵢ³ mod 4` ODD (not pushed to even = a CZ/sign collapse)?**  Corollaries: (i) the substrate must
be the true ring with a NON-DEGENERATE trace/Gram form (GR(4,2): det = unit) — supplies real distance
AND order-4 phase, unlike GF(4) (Z₂ collapse) or full-ring Z₄^m (Gram ≡ 0, distance 1); (ii) γ
accounting is convention-load-bearing — per-encoded-qudit vs per-magic-unit can differ by ~1.

## Top-3 surviving ideas (21 generated, 3 survived adversarial vetting)

1. **GR(4,2) RM + transversal trace-CCZ `i^{Tr(uvw)}`** (top pick).  The only idea whose code half is
   already realized (the verified γ<1 codes).  Trace trilinear form in Z₄ (6/8 tensor entries odd ⇒
   order-4), decisively distinct from the GF(4) Z₂-trace collapse.  Remaining gap = the logical-hosting
   question above.
2. **Z₄ cubic-triorthogonal CSP** — directly SMT/Hensel-solve a native-Z₄ code (n<3k) + a 3-uniform
   hypergraph E so `∏_E CCZ_ring` closes the codespace with odd associator.  Strongest "dodges-closed-
   route" (genuinely 3-local, attacks the exact k≥3 regime).  Existence now partly relieved by the
   verified GR codes.
3. **Tri-block transversal physical ring-CCZ** — three Z₄ CSS blocks, ring-CCZ across corresponding
   triples.  Ranked last: honest γ accounting per magic-unit ≈ 1.63, and the base code is unverified.

## The decisive next experiment (shared crux of all three)

**Test whether the trace-CCZ induces a genuine LOGICAL ring-CCZ.**  On a small dense-feasible code
(n_phys = 6, 4⁶ = 4096) that hosts a ring-CCZ triple: apply `G = diag(i^{Tr(u₁u₂u₃)})`, measure leakage
‖(I−Π)GΠ‖; if 0, EXTRACT the induced logical and certify (i) `level_of = 3` and (ii) the logical
associator `Σ vᵢ³ mod 4` is an odd unit, using the hardened odd-commutator oracle ON THE INDUCED LOGICAL.
Physical-gate genuineness is established; **logical-gate genuineness has never been tested.**  This one
test either validates the multi-qudit ring-magic premise or kills it via logical-level collapse.

Then (the program's flagged make-or-break IOU): a **polynomial / Heisenberg multi-qudit-CCZ hosting
criterion** (the mod-8/16 ring-triorthogonality condition) so the actual γ<1 codes — `[[10,4,3]]`,
`[[12,5,3]]`, dense-infeasible at 4¹⁰/4¹² — can be tested for hosting without dense simulation.

## Open questions (from synthesis)
1. Does any γ<1 GR code actually HOST a transversal trace-CCZ (leakage ~0)?  (Single-qudit died here.)
2. Is there a multi-qudit analog of "hosting ⟺ confined support ⟹ d≤2" — does the distance cap reappear
   one level up?
3. Can the mod-8/16 hosting condition be made polynomial (symbolic) for 4¹²-scale codes?
4. Does self-orthogonality generically force the logical associator EVEN (collapse), or is there
   logical-coset slack that keeps it odd?
5. CSP route: minimal n at which a hypergraph CCZ closes at k≥3.

## Honest flags
- Code genuineness = the hardened anti-collapse cert (≥1 logical pair with odd commutator).  Whether ALL
  k logical qudits are genuine order-4 (vs a soft mix) is a refinement not yet checked — matters for the
  exact γ accounting if any logical is a qubit.
- The whole γ<1 *MSD* claim is contingent on the gate-hosting question, which is unsolved.  As of now:
  **γ<1 genuine ring CODES are verified; a γ<1 genuine ring MSD CODE (code+gate) is not.**

## The decisive experiment — RAN (scratchpad/dexp_gate.py, decisive.py, decisive_search.py, witness_ccz.py)

**STEP 0 (gate genuineness, verified):** the ring-CCZ `CCZ_ring=diag(i^{xyz})` on three Z₄ qudits is
**level-3, order-4** (multi-qudit hierarchy oracle).  Cleanly separated from the collapsed `(-1)^{xyz}`
(level-3 but **order-2** = forbidden qubit-content magic) and from Clifford CZ (level-2).  So the
unifying principle's assumption "physical ring-CCZ is structurally genuine" is now CONCRETE, not assumed.

**STEP 1 (the central worry — does encoding collapse the logical? — REFUTED):**  For a transversal ring-
CCZ across 3 copies of a Z₄ code (logical-X gen g), the induced logical associator is `m = Σᵢ gᵢ³ mod 4`
(odd ⇒ genuine logical ring-CCZ; even ⇒ collapse to CZ/sign).  Hosting conditions (derived + **cross-
validated vs dense, 35/35 then 26/26 with Z-stabs**): `Σgᵢ²sᵢ ≡ Σgᵢsᵢs'ᵢ ≡ Σsᵢs'ᵢs''ᵢ ≡ 0 mod 4` over
M_X generators.  Search found **4840 codes that host with ODD associator** (n0=2–12).  **⇒ Self-
orthogonality does NOT force the logical to collapse.  Genuine multi-qudit ring magic SURVIVES encoding.**
The workflow's #1 feared failure mode is dead.  The multi-qudit direction is ALIVE.

**STEP 2 (the next obstruction — distance — INCONCLUSIVE, leaning hard):**  A random proper-CSS witness
search (M_X + M_Z, k0=1, real cap'd distance; 240k trials, n0≤8) found **NO hosting+genuine code with
d0≥2** — every CSS tri-block hoster is distance-1.  CAVEAT: random CSS codes are rarely d≥2 anyway, so
this is **suggestive, not proven** — it may be the program-wide hosting-vs-distance tension reappearing
"one level up" (synthesis open-Q #2), or just a weak random search.  Do NOT cite as a closed obstruction.

### Net verdict of the decisive experiment
The genuineness/collapse question — the one that would have killed the whole multi-qudit premise — is
answered: **NO collapse.**  The binding constraint has moved (again) to **distance + hosting
coexistence**.  The CSS tri-block route shows the mechanism but stalls at distance; the live target is
whether a **NON-CSS high-rate code (the verified γ<1 graph codes `[[10,4,3]]`/`[[12,5,3]]`) can host a
multi-qudit ring-CCZ** — which needs a **non-CSS multi-qudit hosting criterion** (the CSS H1-H3 derived
here do not apply to graph codes).  That criterion is the make-or-break next step.
