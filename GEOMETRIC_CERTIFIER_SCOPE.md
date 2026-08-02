# Geometric / flat-occupancy distance certifier — scope + Phases 1–2

> **Headline result.** The certifier proves **d(S) ≤ 6 for every qutrit m=5 code with k ≥ 3** —
> so `[[230,13,6]]₃` is *distance-optimal* and **d ≥ 7 is impossible** at m=5. The original "cover
> the d≥7 regime" goal is answered: there is no d≥7 code there to find.

**Goal.** Compute the distance of punctured-RM triorthogonal codes from their *geometry*, not
column enumeration — covering the **d≥7** regime where the meet-in-the-middle certifier (capped
at d≤6 by its C(n,⌈d/2⌉) weight tables) and Brouwer–Zimmermann both die. (BZ was benchmarked and
ruled out: it scales with code dimension K=n−R, huge for our small-redundancy codes; GUAVA's fast
`MinimumWeight` is GF(2)/GF(3)-only and stalls >8 min on `[[230,13,6]]` where our MITM does
seconds.)

## The identity

Distance = minimum *punctured* weight of a dual RM codeword:

    d(S) = min_{c ≠ 0 ∈ RM_p(r̃,m)}  |supp(c) \ S|  =  min_c ( wt(c) − |supp(c) ∩ S| ).

RM min-weight codewords are **geometric** (Delsarte–Goethals–MacWilliams / Leducq: products of
affine hyperplanes → indicators of affine flats), so the low-weight part of this min is a finite,
structured computation. For qutrit m=5 the capping code is **RM₃(6,5), d_RM=9**, and its
min-weight codewords sit on **2-flats** (r̃=6 ⇒ a=3, span j=m−a=2).

## Phase 1 — upper bound / fast screen (BUILT + VALIDATED)

`qmsd/structured_distance.py`: enumerate flat-supported codewords of affine span ≤ jmax, take the
min punctured weight → **d_upper ≤ d(S)** (a *certified* upper bound — it exhibits an explicit
codeword). **Exact** whenever the distance-binding codeword is flat-supported.

For the qutrit m=5 codes the j=2 term is exactly

    d = d_RM − max over 2-flats of |plane ∩ S|

— the flat-occupancy law we found empirically. **Validated exact on all four qutrit m=5 oracle
codes** (test_structured_distance.py), agreeing with the independent MITM:

| code | d_RM − max-2flat | d_upper | true d (MITM) |
|---|---|---|---|
| `[[230,13,6]]` | 9 − 3 | 6 | 6 |
| `[[215,28,5]]` | 9 − 4 | 5 | 5 |
| `[[206,37,4]]` | 9 − 5 | 4 | 4 |
| `[[200,43,3]]` | 9 − 6 | 3 | 3 |

**Uses today:** rigorous no-go bounds (d ≤ X ⇒ γ ≥ 1); a millisecond pre-screen in the search
(reject low-d codes before the MITM); and exact distance for the flat-binding arc codes, with no
d≤6 cap. Cost = (#flats) × (tiny per-flat), no C(n,d) scan.

**Correction it already forced.** Running Phase 1 refuted an earlier claim that "d≥7 codes are
abundant at small k." Those plane_spread/flat_spread codes have max-2flat = 3 ⇒ d_upper = 6, and
the MITM agrees (d=6). d≥7 requires **max-2flat ≤ 2**, a stronger 2-flat condition our samplers do
not target; whether such codes exist at m=5 is open (and the certifier is the tool to check them).

## Phase 2 — certified LOWER bound (BUILT) → and the no-go it proves

`weight_hierarchy` + `flat_lower_bound` in `qmsd/structured_distance.py`. The distance splits over
the RM weight classes: the min-weight (span-2) codewords give the **exhaustive** term
`d_RM − max_2flat_occupancy(S)`; every heavier codeword has RM-weight ≥ the **second weight** w₂,
so its punctured weight is ≥ w₂−k. Hence the certificate

    d(S) ≥ min( d_RM − max_2flat_occupancy(S),  w₂ − k ),   w₂ = 12 for RM₃(6,5).

**What Phase 2 actually settled — a rigorous NO-GO, not a d≥7 certificate.** The d≥7 target does
not exist at qutrit m=5: **any 3 puncture points lie on a common 2-flat** (they span ≤ 2 dims),
whose weight-9 indicator is a codeword, so `|supp(c)\S| ≤ 9 − 3 = 6`. Therefore

    d(S) ≤ 6  for EVERY puncture set with k ≥ 3.

So `[[230,13,6]]` is **distance-optimal** and **d≥7 is impossible** at qutrit m=5. This is the
Phase-1 upper bound applied *structurally* (`max_2flat_occupancy(S) ≥ 3` whenever k≥3), and it
definitively closes the "is there a d≥7 qutrit m=5 code?" question — NO. It also retires the old
"d≥7 abundant at small k" claim for good.

**Where the lower bound pins the exact distance (no MITM).** In the small-k regime the two bounds
meet: for `w₂ − k ≥ d_RM − max_2flat`, `flat_lower_bound = d_upper` and the distance is exact.
Cross-checked against the MITM on small sets (e.g. k=5, max-2flat=3 → lower = upper = MITM = 6).

**Weight hierarchy computed.** `weight_hierarchy(3,5,3) = {2: 9, 3: 12}` (span-2 = d_RM = 9;
span-3 = 12). The general **large-k** lower bound (proving w₂ is a true bound for *all* spans 4,5
— no weight-10/11 full-span codeword) still inherits the **full-span crux** (`D_CRUX_REDUCTION.md`),
but the no-go above makes it moot at m=5: there is nothing above d=6 to certify.

## Build status & next steps

- **Done (Phase 1):** upper bound / screen (`structured_distance`, `max_flat_occupancy`), exact for
  the qutrit m=5 arc codes. Generalizes over `structured_ad`'s flat enumeration.
- **Done (Phase 2):** lower bound (`flat_lower_bound`), weight hierarchy (`weight_hierarchy`), and
  the **d ≤ 6 no-go** for all k≥3 at qutrit m=5. 8 tests in `test_structured_distance.py`.
- **Open (other m/p):** the no-go argument is m-specific (it uses that 3 points span a 2-flat and
  the 2-flat indicator sits in RM_p(r̃,m)). For higher m the analogous ceiling is
  `d_RM − (points forced onto the lowest binding flat)`; recompute per (m,p) before claiming a
  distance ceiling elsewhere.
- **Cheap next:** per-flat min-weight for j≥3 via a small restricted-code MITM (tighter upper
  bound when the min-weight term isn't binding); wire `structured_distance` as a pre-screen into
  `qmsd.search`.
