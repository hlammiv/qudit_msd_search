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
| [[1746,270,>=7]]_7 (sweep k = 270 seed 407, d_plane = 8) | F_7^4 | <= 0.959 (0.898 if d = 8) | d >= 7: block weight-6 scan + full weight-5 scan + mixed-word exclusion (2026-10-29, `CERTIFY_P7.md`); no planar or 3-flat weight-7 word found; d = 8 open |
| [[1746,270,7]]_7 (sweep k = 270 seed 19) | F_7^4 | **0.959** | d = 7 exact (2026-10-30): block weight-6 scan, full weight-5 scan, mixed-word exclusion, planar weight-7 witness (`CERTIFY_P7.md`) |
| [[1751,265,7]]_7 (sweep k = 265 seed 8) | F_7^4 | 0.970 | d = 7 exact (2026-10-30): block weight-6 scan + mixed-word exclusion (the full weight-5 scan is implied, see `CERTIFY_P7.md`), planar weight-7 witness |
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

## Multi-seed sweep with the cheap filters (2026-10-28)

The two cheap filters (plane distance, then the class-block weight <= 5 check for draws with
d_plane >= 7) were run on random cube-class output sets of the sphere charge on a 4-core box
(`scratchpad/plane_sweep.py`; results `runs/weighted/plane_sweep/`, summary
`summary.json`). A *survivor* has d_plane >= 7 and no block word of weight <= 5; it is a
candidate for the full certificate (weight-6 block scan, full weight <= 5 scan, weight-7
witness), not a certified code.

| k | n | draws | d_plane 3 / 4 / 5 / 6 / 7 / 8 | plane pass, block fail | survivors | rate | CPU h |
|---|---|---|---|---|---|---|---|
| 250 | 1766 | 15 | 0 / 1 / 1 / 9 / 4 / 0 | 0 | 4 | 27% | 1.3 |
| 255 | 1761 | 22 | 0 / 1 / 3 / 14 / 4 / 0 | 0 | 4 | 18% | 1.3 |
| 260 | 1756 | 23 | 0 / 1 / 5 / 13 / 4 / 0 | 0 | 4 | 17% | 1.2 |
| 265 | 1751 | 600 | 3 / 14 / 91 / 410 / 82 / 0 | 3 | 79 | 13% | 21.3 |
| 270 | 1746 | 600 | 6 / 17 / 102 / 424 / 50 / 1 | 0 | 51 | 8% | 15.3 |
| 275 | 1741 | 63 | 1 / 0 / 9 / 53 / 0 / 0 | 0 | 0 | 0% | 0.5 |

Findings:

- d >= 7 candidates are not rare: one draw in 8 at k = 265 (gamma 0.960 if d = 7) and one in 12
  at k = 270 (gamma 0.959). The [[1756,260,7]]_7 draw was a 1-in-6 event, not a fluke.
- The plane filter does almost all the work: only 3 of 1323 draws passed it and then failed the
  block weight <= 5 check (all at k = 265).
- One draw, k = 270 seed 407, has d_plane = 8 (no planar word of weight <= 7); if its full
  distance were 8 it would be [[1746,270,8]]_7 at gamma = 0.898. The random-count model below
  expects ~27 non-planar weight-7 words at k = 270, so d = 7 is the likely truth.
- **Two populations of low-weight words.** Planar words (geometric, caught by the plane filter)
  grow slowly with k: the survival rate 27% -> 13% -> 8% from k = 250 to 270 is a Poisson mean of
  1.3 -> 2.0 -> 2.5 planar words per draw. *Accidental* non-planar words among the 672 - k class
  columns in F_7^(291 - k) grow like 7^5 per five outputs: the random-count estimate
  C(672 - k, w) 6^(w - 1) / 7^(291 - k) of block words of weight w is

  | k | checks | weight 5 | weight 6 | weight 7 |
  |---|---|---|---|---|
  | 265 | 26 | 1e-8 | 5e-6 | 2e-3 |
  | 270 | 21 | 2e-4 | 0.08 | 27 |
  | 275 | 16 | 3 | 1200 | 4e5 |

  so the window closes between k = 270 and 275. Checked directly (`scratchpad/block_test.py`,
  `runs/weighted/plane_sweep/block_test_k275.txt`): all six k = 275 draws have a block word of
  weight <= 5, four of them non-planar weight-5 words (affine span 3 or 4), two planar weight-4
  words. At k = 275 the plane filter therefore cannot see the words that kill the draw, and the
  k = 275 sweep was stopped.
- **Witnesses.** Every survivor with d_plane = 7 contains an explicit planar weight-7 word; these
  were extracted and verified against the full G0 (`scratchpad/flat_witness.py`, about 2 min
  per draw; `runs/weighted/plane_sweep/flat_witness_k*.jsonl`: G0 columns, coefficients,
  points). So d <= 7 is proven for 129 of the 130 survivors at no lenore cost; the remaining
  one is k = 270 seed 407 (d_plane = 8), whose 3-flat witness scan is still running (2026-10-28; a hit would settle d = 7, a miss leaves d = 8 open).
- Consequence for certification: at k = 265 a survivor can go straight to the full weight <= 5
  scan; at k = 270 the ~8% chance of an accidental weight-6 block word means the 3+3 weight-6
  block scan (6 GB table, ~10 min on lenore) should come first.

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

| k | n | residual checks | draws | d | gamma |
|---|---|---|---|---|---|
| 180 | 930 | 25 | 1 | 5 | 1.020 |
| 185 | 925 | 20 | 3 | 5, 5, 5 | 1.000 |
| 190 | 920 | 15 | 2 | 4, 3 | 1.138, 1.436 |
| 200 | 910 | 5 | 1 | 2 | 2.19 |

So the r = 6 family reaches gamma = 1.02 at k = 180 and exactly 1.000 at k = 185 (n/k = 5 = d, three of three draws) but has no block lemma to lean
on: with only 205 - k checks over ~930 columns, the distance collapses quickly above k = 180
(d = 4 at 15 checks, d = 2 at 5 checks). It is a near miss, not a route below the sphere's
0.982, unless a charge with a vanishing-stabiliser block at r = 6 is found (the sphere * x1 x2 x3
charge has 2 such stabilisers; sphere (sphere+1) x1 has 16 but a larger support).

## Do weighted supports reopen p = 3, 11, 13? (2026-10-28, exploratory)

The no-go of `D_NEGATIVE_RESULT.md` is for unweighted evaluation codes; the p = 7 code is weighted.
A first pass at the other small primes (`scratchpad/weighted_census_general.py`,
`scratchpad/flat_cap_check.py`; logs in `runs/weighted/census_p_le_13/`):

1. **Census + random-count model.** For every (p, m, r) with a non-negative charge budget and
   every diagonal quadric / linear / diagonal-cubic charge and charge class: |X|, class rank
   rho, and the largest k at which accidental block words of weight <= w are absent, with the
   gamma that d = w + 1 would give. The model is necessary, not sufficient: it ignores
   geometric (flat-supported) words, which are the binding constraint at p <= 13.
2. **Exact flat check** on random draws at the model's k: minimum weight of a word supported on
   an affine line / plane (a rigorous upper bound on d).

| p, m | r | best model window | needed d | exact flat words at that k | verdict |
|---|---|---|---|---|---|
| 3, 6 | 3 (budget 2) | none (rho = 78, |X| >= 432: n/k >= 4.5 at k = rho) | - | - | no window |
| 5, 4 | 4 (budget 3; r = 5 has budget 0) | none (rho = 70, |X| >= 475) | - | - | no window; the unweighted r = 5 [[519,106,5]]_5 stands |
| 13, 2 | 7 | parabola x1^2 + x2, k = 19, gamma 0.95 | 8 | lines: weight 4-5 (5/5 draws) | dead |
| 11, 3 | 8 | cubic sum x_i^3 + c, k = 138, gamma 0.93 | 9 | lines: weight 5 (2/2 draws) | dead |
| 13, 3 | 11 | quadric x1^2+x2^2+x3^2+1, k = 221, gamma 0.96 | 9 | lines: weight 5 (1 draw, 5267 dependent lines) | dead |
| 3, 7 | 4 (budget 1, linear charges) | x7 + c (X = two parallel hyperplanes, 1458 pts, rho = 246), k = 196, gamma 0.957 at d = 7 | 7 | full-code MITM: d = 5, 4 (2/2 draws; `scratchpad/p3m7_probe.py`) | dead |

**Why lines kill p = 11 and 13.** RM(r, m) restricted to a line is RM(r, 1), the full function
space on the p points only when r >= p - 1. For r < p - 1 every line carries dual words (weight
p - 1 on a full line, less after removing the outputs and the charge zeros on it), and a
weighted support leaves lines essentially full. At p = 7, m = 4 we had r = 7 >= 6, so no line
words; at p = 11 and 13 the condition r >= p - 1 fails at m = 3 (r_max = 9, 11) and holds only at
m = 4, where n > 10^4 and no exact certificate is in reach. So weighted supports do not reopen
p = 11, 13 at any size we can certify. At p = 3 and 5 the class ranks are too small for a window
at m <= 6 / m <= 4 (every unit is a cube at p = 3, 5, 11, so there is a single class and no
class-block structure either).

## Caveats

- Cost figures are leading-order surrogates on the standard-library baseline; the exact
  per-output curves (7^{dim G0} words) are out of reach, and the Lemma-17 subcode bounds are
  too loose at this size.
- The weighted codes consume three species of magic state (U^μ for μ in the three cube
  classes). Under a stricter single-species accounting the comparison would need a conversion
  cost.
- d of a draw varies strongly at fixed k; the certified codes are single draws, not optima.
