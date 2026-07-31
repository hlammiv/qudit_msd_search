# Direction D — The non-RM triorthogonal frontier (γ<1 small code at p ≤ 13)

**Status of the surrounding program (established, do not re-litigate):**
The flat-cap no-go has CLOSED punctured Reed–Muller for p ≤ 13 (p=11,13 airtight; p=7 closed
rigorously only for high-k k≥313 — the **mid-k m=4 window k∈[110,312] is NOT certified** and is
one of D's prime targets, see §2.3;
p=5 only crosses at m=4 via the flagship `[[519,106,5]]`, γ=0.987, with no smaller code;
the first genuine γ<1 *small* codes appear only at p=17,19 RM). Direction D asks the one
remaining question: **is there a non-RM triorthogonal family that crosses γ<1 with a SMALL
code at some p ∈ {5,7,11,13}?** This is high-risk and, if it fires, its own paper.

This plan fuses three scout briefs (AG/non-RM evaluation codes; direct triorthogonal-matrix
search; obstruction theory + feasibility filter). It is **kill-criteria-first**: rule out the
hopeless p before searching, and structure every probe to *trigger* a kill fast.

> **EXECUTION STATUS (2026-06-30, see `D_PHASE01_FINDINGS.md` for full results).**
> Phases 0–1 run. **Three routes CLOSED:** Rank 1 (cyclic/BCH — non-RM triorthogonal codes exist but
> collapse to d≤2 under puncturing, γ≥1.65), Rank 3 (AG — not triorthogonal over prime F_p, genus 1&2),
> Rank 4 (Artin–Schreier — degenerate over the prime field). The §1/§6 thesis (**cap is generic to
> triorthogonality over prime F_p, not RM-special**) is now directly evidenced. **Re-scope:** Filter C
> shows p=11 misses γ<1 by only **+1** at its optimal k=21 (not the "+4 / near-hopeless" below), co-equal
> with p=13. **Next:** Rank 2 (moment-ILP, Phase 2) — the last genuine search route. The realistic
> deliverable is now the publishable negative.

---

## 1. Goal and success criterion

**Definitions** (repo conventions): `n = (#eval points) − k`, `k` = #logical qudits (#punctures),
`d` = certified minimum distance, `γ = ln(n/k)/ln(d)`. γ<1 ⟺ **d > n/k**.

**Primary success (a result):** a *certified* triorthogonal CSS code over prime F_p
(p ∈ {5,7,11,13}) that is **not monomially equivalent to a punctured Reed–Muller code**, with

- `is_triorthogonal(basis, p) == True` on its defining generator (the gate is logical), AND
- a non-Clifford logical action (at least one logical row has non-vanishing cubic self-moment
  `Σ_i ℓ_i^3 ≢ 0 mod p` — see §3), AND
- distance `d` certified by `min_dependent_columns(X_stab, p, d_max)` (exact, adversarially
  validated for d ≤ 6), AND
- `γ = ln(n/k)/ln(d) < 1`, with the code **smaller** (in n) than anything known at that p.

A *weaker but still publishable positive* is a non-RM triorthogonal family that ties or
slightly beats RM in A_d / n at fixed γ (use `structured_ad` / `weightcount` for A_d).

**Publishable negative result:** a clean structural no-go that the **flat/2D cap is generic to
triorthogonality over prime F_p**, not special to RM — i.e. any triorthogonal evaluation code
over F_p (p ≤ 13) at the γ<1 puncture density admits a low-dimensional flat/subvariety the
puncture set over-covers, forcing `d ≤ n/k`. This would close the entire evaluation-code
paradigm for p ≤ 13 and reduce the program to "RS at p ≥ 23 + the p=17/19 RM codes." It mirrors
the now-airtight RM flat-cap closure and is a real contribution.

---

## 2. The obstruction to beat, and the per-p feasibility filter (KILL FIRST)

### 2.1 The flat-cap, stated crisply
Min-weight codewords of a triorthogonal evaluation code's relevant dual concentrate on a
low-dimensional **affine flat** (line/plane/hyperplane). For RM over AG(m,p) the dual
min-weight words are products of affine forms, supported on flats, giving
`d = d_RM − max_flat |S ∩ flat|`. Because AG(m,p) is tiled by codim-1 flats, any size-`k`
puncture set over-punctures some flat by `~k/p^{m−1}`, collapsing d. The controlling tension
(AG brief): **triorthogonality is a character-sum / group-annihilation identity; the flat-cap
is a consequence of the same affine-group structure. Breaking affine-invariance to escape the
cap destroys the automatic character-sum vanishing that gives triorthogonality.** Demonstrated:
every n≡0-mod-13 elliptic AG code fails `is_triorthogonal` on the higher power-sums, not just
the constant row.

### 2.2 Three nested filters — only the third binds
- **Filter A (Singleton, necessary, LOOSE):** γ<1 impossible unless `D > 2√N − 2`
  (`N = p^m`, `D = dim RM(r_max,m)`, `r_max = ⌊(m(p−1)−1)/3⌋`). First-feasible m: p=5→m=4,
  p=7→m=3, p=11→m=2, p=13→m=2. Predicts "lots of room" at p=11,13 — **wrong obstruction.**
- **Filter B (line+arc, m=2, tighter):** d ≤ d_RM − t with the puncture set forced to be a
  t-arc (≤ (t−1)p+t points; no maximal arcs over prime F_p). Kills p=5 m=2 and p=7 m=2 but
  still (mis)predicts p=11/13 are easy. **Also not binding.**
- **Filter C (the ACTUAL constraint, 2D/flat cap, VERIFIED via `min_dependent_columns`):**

  | p | m | k (γ<1 density) | line bound d | **true d (2D)** | γ |
  |---|---|---|---|---|---|
  | 11 | 2 | 25 | 5 | **2** | 1.94 |
  | 13 | 2 | 29 | 6 | **4** | 1.14 |
  | 17 | 2 | 52 | 7 | 6 | **0.847 ✓** |

  The discriminator p≤13 vs p≥17 is whether line-supported words are the *global* minimum at
  the γ<1 density — first true at p=17 m=2. p=11 collapses catastrophically (d=2); p=13 misses
  by exactly **one distance unit** (true d=4, need 5).

### 2.3 Per-p go/no-go verdict from the feasibility filter
- **p=5 — LOW priority / SKIP.** Already crosses via RM (m=4). Non-RM only buys a smaller code;
  modest payoff, weak obstruction. Probe only if a generic engine exists at zero marginal cost.
- **p=7 — BEST TARGET.** No RM code crosses (m=3 plane-capped, m=4 hyperplane-capped — rigorous
  only for k ≥ 313). The **mid-k window k ∈ [110,312] at m=4 is NOT closed** (`min_dependent_columns`
  overflows there). Smallest field where a non-RM lever could open new ground; Hasse–Weil over F_7
  gives small, MITM-exact test codes (N ≤ 13). **GO.**
- **p=11 — near-HOPELESS.** 2D cap catastrophic (d=2 at γ<1 density); a non-RM family must beat
  RM by ~4 distance units. **NO-GO unless p=7/p=13 produce a transferable mechanism.**
- **p=13 — SECOND TARGET, borderline.** Misses by exactly +1 distance unit. The smallest gap to
  close; any non-RM family recovering +1 over the RM 2D cap crosses. **GO (after p=7 signal).**

---

## 3. The hard constraint: qudit triorthogonality / CSS-T condition mod p

Any candidate generator `B` (K×n over F_p), whose row-span we want triorthogonal so a transversal
diagonal third-level gate is logical (`qmsd/triorthogonal.py::is_triorthogonal`, NOTES Def 1/Thm 1),
must satisfy, **including repeated indices**:

- **Pair (self-orthogonality):** `Σ_i B[a,i]B[b,i] ≡ 0 (mod p)` for all a ≤ b  — `K(K+1)/2` eqns.
- **Triple (cubic):** `Σ_i B[a,i]B[b,i]B[c,i] ≡ 0 (mod p)` for all a ≤ b ≤ c  — `C(K+2,3)` eqns.

These are GL(K,p)×S_n-gauge-invariant as a *space* property (multilinearity), so the search gauge
group is GL(K,p) on rows × S_n on columns.

**Non-triviality (magic) caveat — the puncture path hides this:** pair=triple=0 makes the gate
logical but *Clifford-trivial* unless some logical row carries a non-vanishing cubic self-moment.
Split rows into stabilizer rows G0 (pair=triple=0 against everything) and k logical rows ℓ with
`Σ_i ℓ_i^3 ≢ 0 mod p` (the cubic self-product is the logical phase) while pair=0 and triple=0
against G0. The CSS code is CSS(G0→X, G'^⊥→Z). Use `structured_pe` logical-phase normalization
sanity (Campbell–Anwar–Browne convention) when asserting the gate is non-Clifford.

**The one genuinely new tractable handle (direct-search brief):** read columns of B as points
x_i ∈ F_p^K; the pair/triple conditions say the symmetric 2- and 3-tensor *moments* of the column
multiset vanish, and **the moments are LINEAR in the column-multiplicities n_x**:
`Σ_x n_x (x_a x_b) ≡ 0` and `Σ_x n_x (x_a x_b x_c) ≡ 0`. So "is there a triorthogonal code on
this column pool" is an **integer-linear feasibility problem mod p**, not a cubic-variety search.
RM is the full-grid solution (symmetry kills all degree-≤(p−2) moments). Direct search = find a
*smaller, non-grid* multiset that still annihilates the 2nd+3rd moment tensors yet spreads support
to keep distance high. This is the lever that makes Phase 1 cheap.

---

## 4. Candidate approaches, ranked by expected value

All approaches reuse the **family-agnostic** core: `is_triorthogonal` (constraint oracle),
`puncture_matrix` / `shorten_matrix` / `dual_matrix` (CSS assembly — note `build_triorthogonal_code`
accepts a prebuilt `G=`, but its `3r < m(p-1)` assert is RM-specific, so for non-RM generators call
the lower-level `shorten_matrix`/`dual_matrix` directly), and `min_dependent_columns(X_stab, p, d_max)`
(exact distance certifier, d ≤ 6; `_parallel` variant for n_jobs), plus `structured_pe.line_punctured_distance`
(cheap flat upper bound + witness — the cap diagnostic), `structured_ad.structured_ad` / `weightcount`
(A_d for survivors). Distance always capped at d_max ≤ 6 (MITM budget; do NOT run the OOM d=6 balanced count).

### Rank 1 — BCH / cyclic / duadic codes, n | p^s − 1  (HIGHEST EV, cheapest, 100% existing toolkit)
**Why:** retains a *multiplicative* group identity (Σ over n-th roots of unity = 0 when n | p^s−1),
so the zero-sum subspace is explicit (defining-set / cyclotomic-coset structure), yet for most
lengths the code is **not affine-invariant** ⇒ min-weight words are NOT flat-supported ⇒ the
`structured_pe` line/plane cap argument simply does not apply. Triorthogonality becomes a
**closure condition on the defining set under addition of exponents mod (p^s−1)** — finite,
checkable, and a much larger design space than RM at the same length (this is the CSS-T-from-cyclic
line of work).
**First experiment:** p=11, n = p²−1 = 120 (and p=13, n = 168). Enumerate BCH defining sets by
designed distance δ; for each build the generator, run `is_triorthogonal(B, p)`, and if it passes
certify d via `min_dependent_columns(B-derived X_stab, p, d_max≤6)`. Sweep k = n − dim, read off γ.
**Tooling:** `triorthogonal.is_triorthogonal`, `mindist.min_dependent_columns(_parallel)`,
`structured_ad` for A_d; generator construction is a small new helper (cyclotomic cosets over F_p)
— no curve theory.
**Obstacle:** the mod-p triple-divisibility tends to force the code into the dual-containing /
RM-like corner; genuinely non-RM cyclic survivors may have small d and land at γ ≥ 1.
**Kill:** if a full defining-set sweep at n ∈ {120,168} yields no triorthogonal code with γ<1, and
every triorthogonal one is affine-invariant (RM-equivalent — compare weight enumerators via
`weightdist`), the cyclic route collapses to RM ⇒ kill Rank 1.

### Rank 2 — Moment-ILP direct search on a structured column pool  (HIGH EV, the new lever)
**Why:** §3's moment-linearity turns the cubic variety into linear feasibility mod p over a chosen
*small* candidate pool. RM is one solution; seek a smaller non-grid multiset.
**First experiment:** p=11 (or 13), pick K = dim G0 = 5–6 (small — pool size p^K/(p−1) blows up past
K≈7). Candidate pool = degree-≤2 monomial evaluation vectors in F_p^K **minus the full grid** (a
sub-RM pool), or a 200–400-point cap in F_p^K (`sampling.random_cap`, `sampling.is_cap`). Solve the
linear moment system (pair+triple ≡ 0 mod p) for a 0/1 selection n_x (ILP feasibility — fast because
linear). For each feasible selection assemble B, assert `is_triorthogonal`, build CSS, certify d via
`min_dependent_columns(X_stab, p, d_max=4)`, compute A_d via `structured_ad`. Hill-climb with column
swaps that preserve the moment kernel, maximizing d — reuse `search.random_search(sampler="capset_climb",
n_jobs=...)` and `sampling.cap_extends`/`collinear`.
**Tooling:** `sampling.*`, `search.random_search`, `triorthogonal.is_triorthogonal`,
`mindist.min_dependent_columns`, `structured_ad.structured_ad`. ILP/feasibility solver is a small
new helper (mod-p linear system; do NOT reject-sample raw matrices — P ≈ p^−[K(K+1)/2+C(K+2,3)]).
**Obstacle:** the same flat/2D cap is **generic** — any moment-feasible multiset risks an affine flat
over-concentrating support; and at K ≳ 7 the pool forces you back onto structured (RM/AG) pools.
**Kill (the likely structural one):** if every Rank-2 survivor at the γ<1 density exhibits a surviving
weight-d 2-flat/3-flat codeword (detect with `min_dependent_columns` + `structured_pe` witness), the
2D-ceiling argument generalizes to arbitrary triorthogonal F_p codes ⇒ **publishable no-go** (§1) ⇒
kill Rank 2 (and the "stay-affine" sub-route). Budget kill: 10⁵–10⁶ feasible configs, best γ ≥ 1.

### Rank 3 — One-point AG / Goppa codes on low-genus curves over prime F_p  (MEDIUM EV, the real "escape")
**Why:** min-weight dual words are zero-divisors of functions in L(βG) on a genus-g curve — supported
on curve points, NOT affine flats ⇒ no "line absorbs many punctures" ⇒ the `d = d_RM − max_flat` cap
has no direct analog. Optimistic Goppa+Hasse–Weil gives γ<1 at tiny n for all p≤13 (e.g. p=7 g=1
`[[9,4,9]]` γ≈0.37) — **but this is the same dream Filter A sold and will likely be defeated by the
AG analog of the flat-cap** (rational points re-cluster on lines/conics over small prime fields).
**The load-bearing unknown:** generic AG over prime F_p FAILS triorthogonality on higher power-sums
(demonstrated for all F_13 elliptic curves). The fix requires a **constant-residue differential** /
self-orthogonal "Castle" structure (Hermitian/Suzuki/Ree-type) so all triple-product sums vanish at
once — but those maximal curves live over F_{q²}, where each coordinate is a 2-qudit ⇒ the settled
**field-CCZ collapse**. Over prime F_p such curves are scarce.
**First experiment (cheapest, MITM-exact, do FIRST as the decisive K1 test):** elliptic curve over F_7,
`y² = x³ + ax + b`, pick (a,b) maximizing N (≤ 13 by Hasse–Weil), one-point code L(aP_∞) with a=4
(cubic budget 3·4=12 < 13). Build the Riemann–Roch evaluation generator (new helper), then:
(1) `is_triorthogonal(B, 7)` — **if this fails, the whole AG-over-F_7 route is dead, not just this curve;**
(2) `min_dependent_columns(X_stab, 7, d_max)` at γ<1 density (k≈4–6) vs optimistic Goppa d = N−a;
(3) cap diagnosis — if d collapses, locate the low-weight support and check whether it lies on a
line/conic ∩ curve (the AG cap recurring).
**Pre-screen that avoids curve theory (AG brief, family A):** for a given point set compute the mod-p
**zero-sum monomial subspace** `Z = {monomials μ : Σ_P μ(P) ≡ 0 mod p}` from the eval matrix, and test
whether a low-degree L(αG) exists with `L(αG)^{·3} ⊆ Z`. Pure linear algebra over the eval matrix —
feed any survivor to `is_triorthogonal`, then `shorten_matrix`/`dual_matrix` → `min_dependent_columns`.
**Tooling:** new Riemann–Roch / point-enumeration helper; everything downstream is existing toolkit.
**Kill:** (K1) the F_7 elliptic code fails the triple-product test ⇒ AG route dead. (K2) the AG sub-curve
cap recurs (true d ≈ N − a − (curve ∩ low-degree intersection)) for genus 1–3 over F_7 and F_13 ⇒
non-RM gives no advantage ⇒ p ≤ 13 closed for evaluation codes. (K3) genus needed for γ<1 grows faster
than Hasse–Weil allows over F_p (genus penalty cancels distance gain) ⇒ asymptotically infeasible.

### Rank 4 — Artin–Schreier / additive-action curves `y^p − y = f(x)` over F_p  (MEDIUM novelty, RM-guard needed)
**Why:** the additive automorphism y→y+c and `y^p − y ≡ 0 ∀y∈F_p` give an affine point count
`p·#{x: f(x)=0} ≡ 0 mod p` automatically — restoring a *group* identity that elliptic/hyperelliptic
curves lack, while keeping genus ≥ 1. The (F_p,+)-action telescopes power-sums, enlarging the zero-sum
subspace Z (the thing missing in Rank 3).
**First experiment:** p=5,7, build `y^p − y = f(x)` for several f, evaluate L(αP_∞), run `is_triorthogonal`,
then **guard against RM-in-disguise**: compare weight enumerators (`weightdist`) to the matching RM code.
**Obstacle:** the fibration may make the code decompose as a tensor/sum of RM pieces (re-importing the
flat cap), or have too few F_p-rational points (n ≈ p).
**Kill:** if every triorthogonal Artin–Schreier code is RM-equivalent (same weight distribution) or has
n ≈ p, the additive lever buys nothing ⇒ kill.

### Rank 5 — Projective/generalized RM or cap-puncture (LOW EV, run only as a fast falsifier)
**Why/Obstacle:** stays affine-invariant-adjacent; PRM inherits flat-supported min-weight structure
(Sørensen/Lachaud) ⇒ same cap one dimension up. The repo's 2D ceiling (p=17,19 refutation) already
shows plane codewords cap d regardless of line-spread.
**First experiment:** p=11,13 m=3, replace flat-puncture with a cap/Hermitian-surface puncture
(`sampling.random_cap`), certify with `min_dependent_columns`.
**Kill:** exhibit a surviving weight-(d_RM − cap bound) 2-flat/3-flat codeword (as in p=17/19) ⇒ the
flat ceiling is dimension-intrinsic to affine-invariant codes ⇒ formally closes the stay-affine escape.

---

## 5. Phased milestone plan (compute budget + go/no-go gates)

All re-runs wrapped as `timeout 150 python3 -c '...'`; on timeout fall back to logic-audit + the relevant
pytest file. Never run anything marked >2 min / lenore / OOM (the ~5h A_6 stream, the 57-min balanced d=6
count). Distance certification capped at d_max ≤ 6.

- **Phase 0 — Feasibility filter (LOCAL, minutes).** Recompute Filters A/B/C per p (the §2.2 table) and
  fix the target (p, m, k) windows: confirm p=7 (mid-k m=4 window k∈[110,312]) and p=13 (m=2, +1 gap) as
  the live targets, p=11 hopeless, p=5 skip. **Gate G0:** if Filter C verification (rerun `min_dependent_columns`
  at the γ<1 density) contradicts the table, re-scope before any search. Deliverable: locked target list.

- **Phase 1 — Cheapest decisive probes (LOCAL, hours).** Run in parallel because they share zero code:
  (1a) **Rank 1 BCH sweep** at n=120 (p=11) and n=168 (p=13) — `is_triorthogonal` filter then
  `min_dependent_columns`. (1b) **Rank 3 K1 test** — the F_7 elliptic triple-product test (the single most
  decisive cheap move; if it fails, AG-over-F_7 dies immediately) plus the `L^{·3} ⊆ Z` linear screen for
  p=11,13 g≤2. **Gate G1:** any γ<1 certified non-RM survivor ⇒ jump to Phase 3 (validate + write up). If
  both 1a and 1b kill ⇒ AG-over-prime and cyclic are dead; proceed to Phase 2 only for the no-go paper.

- **Phase 2 — Direct moment-ILP + structural no-go (LOCAL hours; lenore only if a survivor needs n_jobs
  distance certification).** Run Rank 2 (moment-ILP + capset_climb) at p=11,13, K≤6. Primary purpose: either
  find a small non-grid survivor OR **trigger the §1 cap-kill** (surviving weight-d 2-flat witness via
  `structured_pe`). Run Rank 5 as the fast affine-route falsifier. **Gate G2:** survivor ⇒ Phase 3; uniform
  cap collapse across Rank 2/5 ⇒ promote the structural no-go to the publishable negative result.

- **Phase 3 — Validate + characterize survivor (LOCAL; lenore for A_d if n large).** For any γ<1 survivor:
  (a) re-certify d with `min_dependent_columns_parallel`, (b) confirm non-RM via weight-enumerator inequality
  vs RM (`weightdist`), (c) confirm non-Clifford logical action (§3 cubic self-moment), (d) compute A_d via
  `structured_ad` / `weightcount`, (e) regression against `oracle.py` + `tests/` to ensure the build still
  reproduces the paper codes. **Gate G3:** all checks pass ⇒ new-code result; else demote to "ties RM".

- **Phase 4 — Rank 4 Artin–Schreier (LOCAL, only if Phases 1–3 leave the additive lever unexplored and a
  small gap remains at p=7/13).** RM-guard via weight enumerators. **Gate G4:** RM-equivalent or n≈p ⇒ stop.

**Compute:** Phases 0–2 and most of 3 are local (system python3, numpy 2.2.6, qmsd importable; seconds–hours).
Lenore (32-core, ssh -p 60022) is reserved for: parallel d-certification of a large survivor
(`min_dependent_columns_parallel`, n_jobs) and any A_d count that is large but NOT in the OOM/57-min class.

---

## 6. Risks and honest per-p likelihood of success

- **Dominant risk — the cap is generic, not RM-special.** The square-code theorems
  (Couvreur / Márquez-Corbella–Pellikaan) say small-Hadamard-square codes are forced toward GRS/AG, i.e.
  exactly the families whose min-weight words concentrate on subvarieties. Triorthogonality (small *triple*
  Hadamard square) plausibly cannot be had without paying the algebraic structure that re-creates a cap.
  Most likely outcome across all ranks: cap recurs ⇒ the **negative result** (§1) is the realistic deliverable.
- **AG triorthogonality is unproven and demonstrably fragile** over prime F_p (generic curves fail higher
  power-sums; the fix forces F_{q²} → 2-qudit field-CCZ collapse). Rank 3/4 are gated by K1 before distance
  even matters.
- **Pool/dimension blowup** caps Rank 2 at K ≲ 6, pushing the search back toward structured (RM/AG) pools —
  the central tension of direct search.

**Per-p likelihood of a γ<1 small non-RM code:**
- **p=5:** N/A for a *new* crossing (RM already crosses); ~low value, mild chance of a smaller code.
- **p=7:** *best, still low* (~10–20%). Open mid-k m=4 window + small MITM-exact AG tests give the cleanest
  genuine opening; the K1 test is decisive and cheap.
- **p=11:** *near-zero* (<5%). 2D cap catastrophic (d=2); would need ~+4 distance over RM.
- **p=13:** *borderline, second-best* (~10%). One-unit gap; an AG/cyclic family recovering +1 over the RM
  2D cap crosses.

**Expected program outcome:** most probable is a rigorous **negative** (cap generic to triorthogonality over
prime F_p, p ≤ 13). The high-value tail is a p=7 or p=13 BCH/AG survivor. Either is publishable; the plan is
ordered so the cheapest probes (Rank 1 BCH, Rank 3 K1) decide the question first.
