# Geometric / flat-occupancy distance certifier — scope + Phase 1

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

## Phase 2 — certified d≥7 LOWER bound (research bet, NOT built)

To certify d ≥ w, only codewords of unpunctured weight **< w+k** can violate it. Classify + verify
exactly the window [d_RM, w+k). For w=7: k=6 → 4 weight-levels, k=8 → 6, k=13 → 11. If every
codeword in the window is flat-classified (DGM/Leducq), the lower bound is exact. **Ceiling = the
full-span crux** (`D_CRUX_REDUCTION.md`): a full-span codeword below w+k escapes the geometry.
Small-k / short-window is where d≥7 can actually be certified.

## Build status & next steps

- **Done:** Phase 1 upper bound / screen (`structured_distance`, `max_flat_occupancy`), exact for
  the qutrit m=5 arc codes; 5 tests. Generalizes over `structured_ad`'s flat enumeration.
- **Cheap next:** per-flat min-weight for j≥3 via a small restricted-code MITM (tighter upper
  bound when the min-weight term isn't binding); wire `structured_distance` as a pre-screen into
  `qmsd.search`.
- **Research:** Phase 2 lower bound — enumerate the [d_RM, w+k) weight classes; gated on the
  full-span crux.
