# Direction (c): structured / analytic families for γ<1 at p=7 — scope

**Status:** free-form search is EXHAUSTED at p=7 m=4 (established this session). γ<1 is now a
**construction** problem, not a search problem. This doc scopes the construction attack.

## 1. Goal, and why search is dead

γ = log(n/k)/log(d) < 1  ⟺  **d > n/k**. At p=7 m=4 (N=2401, n=N−k):

| target | needs |
|---|---|
| γ<1 | d=6 @ k>343, or d=7 @ k>300, or d=8 @ k>267 … |
| beat our best (γ=1.19) | d=6 @ k>189, or d=7 @ k>? — modest |

**Free-form is proven exhausted:**
- d=5 ceiling PROVEN at k=307 (full parallel b=3 scan found a weight-5 witness — not an artifact).
- distance collapses at high k (uniform k=320 → d=2).
- d≥6 certification is ~2.5 hr for ONE set (memory-bandwidth-bound full scan) → a d=6 **hunt** is infeasible regardless.

So d≥6 must be **built in by structure**, and certified **geometrically** (not by MITM).

## 2. The reframe: a covering-design problem

Distance is `d(S) = min over dual codewords c ∈ RM₇(r̃=16,4) of |supp(c) \ S|`. Hence

    d(S) ≥ 6   ⟺   for EVERY dual codeword c:  occupancy(c) := |S ∩ supp(c)|  ≤  wt(c) − 6.

The min-weight dual codewords have **wt = d_RM = 21** (r̃=16 = 2·6+4 ⇒ a=2,b=4 ⇒ (7−4)·7¹ = 21),
and are **flat-supported** (Delsarte–Goethals–MacWilliams / Leducq). So the min-weight constraint is:
**occupancy ≤ 15 on every weight-21 support.**

**The encouraging observation.** For k punctures over 2401 points, the *average* occupancy of a
21-point support is `21·k/2401`: at k=343 that's **≈ 3**, at k=400 ≈ 3.5 — far below the budget of 15.
So the budget is generous on average; γ<1 (d=6 at k>343) is not obviously impossible. The binding
constraints are the two things a *structured* S must control:
1. the **MAX** occupancy over all weight-21 supports (a large-deviation / spread problem), and
2. the **higher-weight / full-span** dual codewords — the full-span crux (§5).

A d=5 code (what we keep finding) has ≥1 weight-21 codeword at occupancy exactly 16. The design goal
is an S that holds every weight-21 support to ≤15 (and higher weights proportionally).

## 3. Directions, ranked

1. **C-OCC (primary): occupancy-bounded structured puncture.** Build S that minimizes the max
   flat-occupancy on the weight-21 (and higher) binding flats; certify geometrically. Even short of
   γ<1 this can yield **d=6 at moderate k** — d=6 @ k=280 → γ=1.13, @ k=250 → γ=1.20 — a real
   improvement over 1.19 that free-form can't reach.
2. **C-GEO (enabler): a p=7-capable geometric certifier.** `structured_distance` computes distance
   from flat-occupancy with NO MITM — it sidesteps the 2.5 hr wall and is EXACT when the binding
   codeword is flat-supported. Currently p=3-tuned and its full-flat enumeration OOMs at p=7 m=4;
   needs a **point-restricted** version (enumerate flats spanned by puncture-point subsets, not all
   flats — the analogue of the qutrit fast path).
3. **C-MAN (secondary): Manhattan-family rate-boost.** The analytic Manhattan family gives d=6 only
   at k=126 (γ=1.615) — worse than free-form. Check whether a Manhattan d=6 code can be perturbed to
   higher rate while tracking distance analytically. Lower promise (the family sits at high γ).
4. **C-LIT: literature / NOTES scan.** Look in Saha–Prakash (arXiv:2510.10852) and the local NOTES
   for other closed-form-distance puncture families reaching high rate at d≥6 (subfield/subcode RM,
   product/coordinate-subspace structures). AG and cyclic are already closed (see
   D_NEGATIVE_RESULT / nonrm-triorthogonal-frontier).

## 4a. GATING RESULT (2026-08-02): NEGATIVE for j=2 — the crux is active at p=7 m=4

Built the point-restricted p=7 geometric certifier (`geometric_distance_dual` in
`structured_distance.py`): per 2-flat, min punctured weight = `min_dependent_columns(dual_matrix(
RM₇(4,2) restricted to F∖S))`, so it works where the codeword-enumeration OOMs. **Validated EXACT
on all four qutrit m=5 oracle codes** (reproduces d = 6,5,4,3).

On `[[2093,308,5]]` (true d=5): **`d_geo(j=2) = None`** — none of the 117 heavy 2-flats (|S∩F|≥15)
carries a codeword capturing the weight-5 binding. So **the binding codeword is NOT 2-flat-supported
at p=7 m=4** — the full-span / higher-flat crux is active, unlike qutrit m=5 where it was exact.

**Deeper reason (connects to the known p=7 hard core):** caps stall at **k≈75** in AG(4,7), so there
is no flat-binding structured substrate at high k — high-k codes are effectively random and
crux-bound, which is *why* geometric certification is loose here. This is the "2D-codeword crux"
(`nonrm-triorthogonal-frontier`, memory). C-OCC needs high k for rate but flat-structure fails past
k=75: **high-rate and geometric-structure are incompatible at p=7 m=4**, so C-OCC is largely blocked.

**Consequence for the plan:**
- C-OCC (2-flat occupancy) does not get traction at p=7 m=4. Optional deeper check: `j=3` (is the
  binding 3-flat-supported? ~1–3 hr, not run — the cap-stall argument closes C-OCC regardless).
- Pivot weight to **C-LIT** (other analytic families) and **higher m** (does the crux relax at
  p=7 m=5, where caps may reach higher absolute k?), or accept γ<1 at p=7 m=4 is crux-blocked.
- `geometric_distance_dual` itself is a keeper: a p-general point-restricted certifier that IS exact
  for flat-binding codes — the right tool if a flat-binding high-rate family is ever constructed.

## 4b. C-LIT RESULT (2026-08-02): all structured/analytic alternatives CLOSED — (c) converges on the crux

Scanned the paper (`literature/2510.10852/`, incl. NOTES) + the full internal theory corpus
(D_NEGATIVE_RESULT, D_CRUX_REDUCTION, D_PROOF_MAP, P7_MIDK_PROPOSALS, D_P7M4_WINDOW, ...).

**Paper-grounded facts:**
- γ<1 at p=7 IS solved analytically — but only by the Manhattan family at **m=13**:
  `[[96448935471, 440074936, 231]]₇` (~10¹¹ qudits, d=231, paper Table 2). γ₀(7)=0.508 asymptotic.
- The Manhattan γ decreases slowly with m: m=4→1.595, m=7→1.248, m=8→1.224, reaching our 1.19 only
  at **m≥9 (block ≥40M)**. So **`[[2093,308,5]]₇` (γ=1.19, block 2401) is the best known p=7 code for
  every block size up to ~40M** — it beats the analytic family across the whole practical range.
- Certified γ<1 exists at moderate block only for **p≥17** (`[[237,52,6]]₁₇` γ=0.847; RS for p≥23).

**Every structured/analytic alternative to Manhattan is CLOSED** (Forced-Grid Theorem + family-by-family,
D_NEGATIVE_RESULT / D_PHASE01_FINDINGS): RM-puncture (flat-capped), Reed–Solomon (inapplicable at p=7,
r_max=1), cyclic/BCH (exist but collapse to d≤2 punctured), AG/Goppa & Artin–Schreier (not triorthogonal
over prime F_p; extension-field versions → multi-qudit field-CCZ collapse), direct-moment/cap
(crux-blocked, caps stall k≈75), projective/GRM (subsumed by forced-grid). The nominally "untried"
C-LIT families reduce to closed cases too: **subfield-subcode RM / product / coordinate-subspace** are
evaluation codes on a (product) grid → forced-grid → only RM puncturing; **concatenation** gives a
γ that is the *mediant* of its components' γ's → can never beat the best component.

**Verdict for (c):** the only genuine open residuals are (a) the p=7 m=4 **mid-k window k∈[127,~250]**
where a d≥6/d≥7 survivor could hide but CANNOT be certified without cracking the crux, (b) the <5%
high-dimensional general-position moment tail, and (c) the crux itself. All three are gated by ONE
open theorem: the **full-span DGM second-weight-class weight-hierarchy** result (D_CRUX_REDUCTION §4–6).
This is a coding-theory NO-GO direction, not a construction — proving it CLOSES all of p≤13 (incl. p=7
m=4) as a corollary; the relevant math is Leducq's second-weight GRM classification (arXiv:1203.5244)
and Dang–Ghorpade min-weight enumeration (arXiv:2504.21816). No construction reachable from here gives
γ<1 at moderate-block p=7.

## 4. Gating experiment (original plan — now executed, see §4a)

**Does the geometric law `d = d_RM − max_flat_occupancy` hold EXACTLY at p=7 m=4?** i.e. is the
distance-binding codeword flat-supported there, as it was for qutrit m=5?

- **Test:** on the known d=5 codes `[[2093,308,5]]` / `[[2094,307,5]]`, does the geometric certifier
  return d_upper = 5 (matching the MITM), with a weight-21 flat at occupancy 16?
- **YES** → geometric cert works at p=7 → C-OCC is viable: build the occupancy-minimizing sampler,
  target max-occ ≤ 15, and certify candidates geometrically (no 2.5 hr MITM).
- **NO** (d_upper > 5; a full-span codeword binds below the flat bound) → the full-span crux dominates
  at p=7 → C-OCC is blocked; pivot to C-LIT / C-MAN or accept the d=5 ceiling.
- **Blocker to clear first:** `structured_distance` full-flat enumeration OOMs at p=7 m=4. Build the
  **point-restricted** flat-occupancy (the binding flats pass through puncture points) so the gating
  test is even computable. This is the first concrete code task.

## 5. The central risk (be honest)

- **The full-span crux** (`D_CRUX_REDUCTION.md`): higher-weight, full-support dual codewords can sit
  below the flat-occupancy bound, so the geometric certifier is only an *upper* bound in general and
  was exact only up to d≤6 for qutrit. If it's not exact at p=7, C-OCC can propose a candidate but
  can't cheaply certify it (falls back to the 2.5 hr MITM for a single confirm — acceptable once, not
  for a hunt).
- **p=7 may have no γ<1 code at all.** Prior search returned 0 γ<1, and the m=2 crux is a known
  obstruction. This program could conclude "γ<1 impossible at p=7 m=4" — which would itself be a
  clean result (a certified no-go via the occupancy argument), not a failure.
- **Realistic near-term payoff** is d=6 at k≈250–300 (γ≈1.13–1.20), a modest beat of 1.19, *if*
  occupancy-bounding works and the binding codeword is flat-supported. γ<1 is the longer shot.

## 6. Concrete first steps (all design / light compute — no 2.5 hr grinds)

1. **Point-restricted p=7 geometric occupancy** (code): enumerate the weight-21 binding flats through
   puncture-point subsets; compute max-occupancy. Validate against the MITM d on the d=5 codes.
2. **Gating check** (§4): is d = d_RM − max_flat_occ exact on `[[2093,308,5]]`?
3. **If viable:** occupancy-minimizing sampler (drive max weight-21 occ toward ≤15) at k=250→343;
   certify hits geometrically; the first d≥6 or clean no-go.
4. **In parallel, cheap:** C-LIT scan of Saha–Prakash / NOTES for higher-rate d≥6 analytic families.

## 7. Tooling inventory (what exists vs. what's needed)

- **Exists:** `structured_distance` / `max_flat_occupancy` (p=3-tuned, full-flat enum), the geometric
  no-go machinery, `distance_is_ge6_parallel` (single-candidate d≥6 confirm, ~2.5 hr — a backstop).
- **Needed:** point-restricted flat-occupancy for p=7 m=4 (weight-21 flats); occupancy-minimizing
  sampler; the weight-21 support characterization (which flat family carries the min weight at p=7).
