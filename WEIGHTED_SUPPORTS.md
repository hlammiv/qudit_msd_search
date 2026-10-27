# Weighted Reed–Muller supports at p = 7: the first γ < 1 search code

*2026-10-27. Data under `runs/weighted/`; search log in `NEXT_STEPS.md` §H; cost tables in
`PERFORMANCE.md`.*

## Result

| code | support | γ | certificate |
|---|---|---|---|
| **[[1756,260,7]]₇** | sphere charge, m = 4, r = 7 | **0.982** | d = 7 exact (three-part proof + witness, below) |
| [[1746,270,6]]₇ | same | 1.042 | d = 6 exact |
| [[1751,265,6]]₇ | same | 1.054 | d = 6 exact |
| [[263,24,5]]₇ | ellipsoid charge, m = 3 | 1.49 | d = 5 (MITM) |
| [[2093,308,5]]₇ (unit space, free-form, for comparison) | F₇⁴ | 1.19 | d = 5 (MITM) |

No paper we reference has a γ < 1 code at p = 7: Saha–Prakash (arXiv:2510.10852) tabulate
k = 1 codes (γ 2.6–3.3) and an analytic [[2331,70,9]]₇ (γ 1.595); the cyclic codes of
arXiv:2610.08012 have γ 1.77–2.1; the doubled [[20,1,4]]₇ of arXiv:2610.29457 has γ 2.16;
Krishna–Tillich reach γ < 1 only from p = 23. The asymptotic optimum γ₀(7) = 0.508 of
arXiv:2510.10852 is an n → ∞ limit, not a code.

On cost (leading-order surrogate against the p = 7 standard library of RS, RM, cyclic and
doubled codes): [[1746,270,6]]₇ saves 4.3× at 1e-2 → 1e-30 and 3.6× at 1e-3 → 1e-30,
[[2093,308,5]]₇ 4.1× / 2.9×, [[263,24,5]]₇ 2.6× / 2.1×. The cost of [[1756,260,7]]₇ is
pending its weight-7 census (todo in `NEXT_STEPS.md`).

## Construction

Charge picture of arXiv:2610.29457 §4: a polynomial g of degree ≤ m(p−1)−1−3r gives a
weighted triorthogonal space RM(r, m) on X = supp g with charges μ = g(x). At p = 7, m = 4,
r = 7 the budget is degree 2. We take the sphere g = x₁²+x₂²+x₃²+x₄²: |X| = 2016, rank 326.
Outputs must sit at points whose charge is in one cube class; we use the class g ∈ {1, 6}
(672 points) and k random points of it as the puncture set S. Inputs at non-cube charges
consume U^μ magic states, counted as one raw state each (the convention of
arXiv:2610.29457 for its own [[20,1,4]]₇). `qmsd/weighted.py` builds (G0, H1, μ); explicit
weighted triorthogonality of the m = 4 support was checked (all triple, pair and single sums
vanish).

## The class-block lemma (why these codes work and how they are certified)

Outputs confined to one charge class lie on the degree-4 surface (g−1)(g−6) = 0, so the 35
polynomials (g−1)(g−6)·RM(3,4) are stabilizers that vanish on every class point. Hence:

1. the rank of RM(7,4) on the class points is 326 − 35 = 291, the ceiling for k;
2. the 672 − k non-output class points ("the class block") see only 291 − k independent
   checks (21 at k = 270, 31 at k = 260), so the low-weight words of the code live there:
   k = 290 gives d = 1 (56 unchecked class columns), k = 280 gives d = 4, k = 270 gives d = 6, and one k = 260 draw gives
   d = 7 (another gives 4: d depends on the draw as much as on k);
3. a word that touches a non-class column has its non-class part forced, by the 35 kernel
   rows, into RM(3,4)^⊥ = RM(20,4), whose minimum distance is 5.

Point 2 makes certification cheap: a meet-in-the-middle on the ≤ 412 class columns replaces
scans on 1756 columns (the 3+3 weight-6 witness takes minutes; the 2+4 scan on the full code
would take weeks at the measured probe rate). Point 3 makes it rigorous.

### Certificate of d = 7 for [[1756,260,7]]₇

1. No weight ≤ 5 in the full code: serial numba MITM for weights 2–4; thread-parallel full
   2+3 scan for weight 5 (`mindist_nb.weight5_exists_parallel`, 12243 s on lenore).
2. No weight-6 word on the 412 class columns: `mindist_nb.weight6_witness_3x3`.
3. No weight-6 word touching a non-class column: by point 3 its non-class part u is a
   weight-5 or weight-6 word of RM(20,4) (weight ≤ 4 is excluded by d = 5) and its class part
   has at most one column. Every full-support weight-5/6 word of RM(20,4) is supported on a
   line: verified by enumerating point sets up to affine equivalence (a frame of the span plus
   the remaining points; 76,242 sets of span 2–4; none carries a full-support annihilator of
   the cubics). All 1,197,504 line words inside the non-class points were then tested and none
   extends to a word of ker G0 (`scratchpad/mixed_w6_check.py`, 428 s).
4. A weight-7 word exists: G0 columns [523, 1045, 1117, 19, 490, 1260, 1556] with
   coefficients (3, 3, 1, 1, 1, 3, 1), found by `mindist_nb.weight7_witness_3x4` (5.8 h) and
   verified locally (rank 6, every 6-subset of full rank).

The puncture set, the charge, both witnesses and the certificate text are in
`runs/weighted/p7_m4_sphere_1756_260_cols.json`.

## Engines added for this (all in `qmsd/mindist_nb.py`, numba, thread-parallel)

| function | split | use |
|---|---|---|
| `count_weight5_per_output` | 2+3 | per-output leading order (surrogate) at n ~ 2000 |
| `weight6_witness_parallel` | 2+4 | weight-6 witness on the full code (slow; superseded) |
| `weight6_witness_3x3`, `count_weight6_per_output` | 3+3 | weight-6 witness / census on a ≤ 450-column block |
| `weight7_witness_3x4` | 3+4 | weight-7 witness with early exit |

Each is validated against the serial kernel or brute force (`tests/test_mindist.py`,
`tests/test_objective.py`). The per-output surrogate (`objective.surrogate_from_matrices`)
matches the exact curves to within 4% where the codes are applied (n·ε ≤ 0.26).

## Geometry of the low-weight words

The witnesses of the certified codes are **planar**: their supports lie in an affine 2-flat of
F_7^4, and the coefficient is constant on each of the two spheres (g = 1 and g = 6) that make
up the cube class.

| word | points | affine span | collinear triples | spheres (g-values) | coefficients |
|---|---|---|---|---|---|
| [[1746,270,6]]_7 weight-6 witness | (3,1,4,6) (1,4,4,6) (2,0,6,6) (0,1,2,6) (5,1,3,6) (3,2,6,6) | 2 | 2 | 6 6 6 6 1 1 | 3 3 3 3 1 1 |
| [[1756,260,7]]_7 weight-7 witness | (3,4,0,2) (2,1,1,4) (4,0,3,4) (5,3,0,0) (3,4,6,1) (0,2,0,5) (0,2,1,6) | 2 | 5 | 1 1 6 6 6 1 6 | 3 3 1 1 1 3 1 |
| k = 260 sample 0 (d = 4), its unique weight-4 word | (6,3,4,1) (1,6,3,3) (1,2,2,5) (6,4,3,6) | 3 | 0 | 6 6 6 6 | 1 3 3 2 |

Why planes are special: the block checks are the 291 - k residual rows of RM(3,4) restricted
to the class points, and RM(3,4) restricted to a plane is RM(3,2) (dimension 10 of 49). A
plane carries on average 672/49 = 13.7 class points, of which (672 - k)/49 = 8.4 survive at
k = 260, so inside one plane the surviving columns are 8-21 vectors in a space of dimension
<= 10, and dependencies among them are far more likely than among a random set of columns in
F_7^31. The weight-7 witness has five collinear triples: three points on each of two lines of
its plane plus one more point. But planarity is not forced: the k = 260 sample-0 draw (d = 4)
has a single weight-4 word (up to scalar) whose four points span a 3-flat and all lie on the
g = 6 sphere.

**Plane-restricted distance.** d_plane(cols) = min over the 139 650 affine planes of the minimum
weight of a code word supported inside that plane (each plane is a <= 47 x (<= 21) sub-problem
with nullity <= 6, solved by enumerating its null space; ~2 min per draw on one core,
`scratchpad/flat_distance4.py`, `runs/weighted/p7_m4_sphere_plane_distance_check.json`). It is
a rigorous upper bound on d, because a word supported in a plane is a word of the code.
Against the nine draws with known distance:

| draw | d (certified / block) | d_plane |
|---|---|---|
| k = 290 | 1 | 1 |
| k = 280 | 4 | 4 |
| k = 275 sample 0 | 4 (block) | 4 |
| k = 275 sample 1 | 5 (block) | 5 |
| k = 270, [[1746,270,6]]_7 | 6 | 6 |
| k = 265 sample 0 | 3 (block) | 3 |
| k = 265 sample 1 | 6 | 6 |
| k = 260 sample 0 | 4 (block) | **6** |
| k = 260 sample 1, [[1756,260,7]]_7 | 7 | 7 |

Eight of nine exact, one over-estimate (the non-planar weight-4 word above). So d_plane is a
cheap, rigorous *rejection* filter (d_plane < 7 rejects a draw in two minutes, against ~0.5 h
on 30 cores for the full 3+4 scan at n = 1756) and a good but not exact predictor; a draw that
passes it still needs the full certificate. Note the generic MITM kernel `weight_d_exists_fast`
only certifies d <= 6 (3+3 halves); it now raises for d >= 7 instead of silently scanning
weight 6, which is how the first version of this check mis-reported the k = 260 draws.

## Beyond the sphere: degree-5 charges at r = 6

The sphere is the only degree-2 charge at m = 4; the next budget, deg g <= 5, allows r = 6
(deg g <= m(p-1) - 1 - 3r = 23 - 18 = 5). The census of degree-5 charges (products of five
affine-linear forms, products of the sphere with linear forms, diagonal quintics and a few
mixed monomials; `scratchpad/census_r6.py`) gives:

| charge g | \|X\| | rank | cube-class points | rank on class | stabilisers vanishing on class |
|---|---|---|---|---|---|
| sphere (deg 2, r = 7) | 2016 | 210 | 672 | 195 | 15 |
| five general hyperplanes | 1110 | 205 | 370 | 205 | 0 |
| x1 x2 x3 x4 (x1+x2+x3+x4-1) | 1111 | 205 | 340 | 205 | 0 |
| sphere * x1 x2 x3 | 1272 | 207 | 424 | 205 | 2 |
| sphere (sphere+1) x1 | 1442 | 209 | 504 | 193 | 16 |
| sum x_i^5 + 1 | 2058 | 210 | 686 | 210 | 0 |
| x1(x1-1)...(x1-4) | 686 | 140 | 686 | 140 | 0 |

At r = 6 no stabiliser vanishes on the class for the product-of-hyperplane charges, so there is
no class block: the k outputs see all 205 - k residual checks. The attraction is the small
support (n ~ 1110 - k against 2016 - k on the sphere), so n/k ~ 4.5 at k ~ 200 and gamma < 1
would need only d >= 5 at k >= 190.

First probes on the five-hyperplane charge (random cube-class output sets, full weight <= 5
certificates, `runs/weighted/p7_m4_r6_hyperplanes_probe.json`):

| k | n | residual checks | d | gamma |
|---|---|---|---|---|
| 180 | 930 | 25 | 5 | 1.020 |
| 185 | 925 | 20 | 5 | 1.000 |
| 190 | 920 | 15 | 4 | 1.138 |
| 200 | 910 | 5 | 2 | 2.19 |

So the r = 6 family reaches gamma = 1.02 at k = 180 and exactly 1.000 at k = 185 (n/k = 5 = d) in single draws but has no block lemma to lean
on: with only 205 - k checks over ~930 columns, the distance collapses quickly above k = 180
(d = 4 at 15 checks, d = 2 at 5 checks). It is a near miss, not a route below the sphere's
0.982, unless a charge with a vanishing-stabiliser block at r = 6 is found (the sphere * x1 x2 x3
charge has 2 such stabilisers; sphere (sphere+1) x1 has 16 but a larger support).

## Caveats

- Cost figures are leading-order surrogates on the standard-library baseline; the exact
  per-output curves (7^{dim G0} words) are out of reach, and the Lemma-17 subcode bounds are
  too loose at this size.
- The weighted codes consume three species of magic state (U^μ for μ in the three cube
  classes). Under a stricter single-species accounting the comparison would need a conversion
  cost.
- d of a draw varies strongly at fixed k; the certified codes are single draws, not optima.
