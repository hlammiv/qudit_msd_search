# Weighted Reed–Muller supports at p = 7: the first γ < 1 search code

*2026-10-27. Data under `runs/weighted/`; search log in `NEXT_STEPS.md` §H; cost tables in
`PERFORMANCE.md`.*

## Result

| code | support | γ | certificate |
|---|---|---|---|
| **[[1746,270,7]]₇** (sweep seeds 19 and 407) | sphere charge, m = 4, r = 7 | **0.959** | d = 7 exact, both draws (block6 + mixed6 + weight-7 witness; `CERTIFY_P7.md`) |
| [[1751,265,7]]₇ (seed 8) | same | 0.970 | d = 7 exact |
| [[1756,260,7]]₇ (the first) | same | 0.982 | d = 7 exact (three-part proof + witness, below); weight-7 census complete |
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
[[2093,308,5]]₇ 4.1× / 2.9×, [[263,24,5]]₇ 2.6× / 2.1×, and [[1756,260,7]]₇ 4.1× / 3.4× (weight-7
census 2026-10-31: A_7 = 24, L_7 = 1.29, c_ε = 4.6e-6; `PERFORMANCE.md`). The γ < 1 code does not
beat the γ = 1.04 code on cost at these targets: both act as the last round at ε ≈ 1e-6, where
d = 6 already overshoots, so the 270-output rate wins. The k = 270 d = 7 draws ([[1746,270,7]]₇,
seeds 19 and 407) have the same n, k as [[1746,270,6]]₇ and one more unit of distance, so at
leading order they can only be cheaper (unless L_7 exceeded ~1e5, against a census-scale
expectation of ~1); their own weight-7 censuses (~33 h each) would make that exact.

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
   extends to a word of ker G0 (`scratchpad/mixed_w6_check.py`, 428 s). The same holds at
   weight 7 (section "Weight-7 words of RM(20,4)" below), which is what a d = 8 claim needs.
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
| [[1746,270,7]]_7 (sweep k = 270 seed 407, d_plane = 8) | F_7^4 | 0.959 | d = 7 exact: block weight-6 scan + full weight-5 scan + mixed-word exclusion at weights 6 and 7, then a weight-7 block word from the 3+4 scan (14 min with the indexed kernels, 2026-10-30); plane distance 8 but a non-planar block word of weight 7 exists |
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

- d >= 7 candidates are not rare: one draw in 8 at k = 265 (gamma 0.970 if d = 7) and one in 12
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

## Weight-7 words of RM(20,4) (2026-10-30)

*Scripts `scratchpad/rm20_w7_classify.py`, `scratchpad/rm20_w7_sample_check.py`,
`scratchpad/mixed_w7_check.py`; library `qmsd/pointsets.py`; data
`runs/weighted/rm20_w7_classification.json`, `runs/weighted/rm20_w7_sample_check.json`,
`runs/certify/p7_k270_seed407_mixed7.json`; tests `tests/test_pointsets.py`.*

The mixed stage of the certificate (point 3 of the class-block lemma) rests on knowing the
supports of the low-weight words of RM_7(20,4) = RM(3,4)^⊥ over the full grid F_7^4: a word of
ker G0 that touches a non-class column has its non-class part, reweighted by (g−1)(g−6), in that
code. For weight ≤ 6 the supports are lines (`mixed_w6_check.py`). The weight-7 case is now
settled the same way, and the answer is again "lines only".

**Lemma (words of weight ≤ 7 of RM_7(20,4)).** Let u ≠ 0 satisfy Σ_x u_x f(x) = 0 for every
polynomial f of degree ≤ 3 on F_7^4, and let |supp u| ≤ 7. Then supp u lies on an affine line
ℓ = {a + t v}, and u is the restriction to ℓ of a polynomial q(t) of degree ≤ 2 (RM(3,1)^⊥ =
RM(2,1) on the 7 points of ℓ): weight 5 when q has two roots, weight 6 when it has one (q linear
or a square), weight 7 when it has none (q a non-zero constant or an irreducible quadratic). Per
line there are 21 / 14 / 22 projective words of weight 5 / 6 / 7; with 137 200 lines the weight
enumerators are A_5 = 17 287 200, A_6 = 11 524 800, A_7 = 18 110 400.

Equivalently (a finite-field Cayley–Bacharach statement for d = 3, n ≤ 2d + 1): *a set of ≤ 7
points of F_7^4 fails to impose independent conditions on the cubics iff it contains 5 collinear
points, and it carries an annihilator of full support iff it is 5, 6 or 7 points of one line.*
In particular no planar configuration (two lines 4 + 3 or 3 + 3 + 1, a parabola, 7 of the 8
points of an ellipse) and no 3- or 4-dimensional 7-point set carries a weight-7 word. The lemma
stops at weight 7: the 8 points of an ellipse x² + y² = 1 in a plane impose only 7 conditions on
the cubics and carry a full-support weight-8 word, so the weight-8 mixed stage would have to
include conics.

**Proof (enumeration, `rm20_w7_classify.py`, 7 min on one core, < 1 GB).** A 7-point set of affine span
j contains a frame 0, e₁..e_j and 7 − (j + 1) further points of its j-flat, so up to AGL(4,7)
every set of span 2, 3, 4 is one of C(46,4) + C(339,3) + C(2396,2) = 163 185 + 6 435 689 +
2 869 210 frame sets (span 1 is the full line). For each set the 35 × 7 evaluation matrix of
RM(3,4) is tested for rank deficiency by a batched F_7 elimination (`pointsets.batched_rank_modp`,
50 000 matrices per 2.3 s); the dependent sets get their null space, a full-support test (no
coordinate is zero on the whole null space; over F_7 this is equivalent to the existence of a
single all-non-zero null vector) and an affine-invariant classification (collinearity profile,
`pointsets.classify_pointset`). An independent second method row-reduces the frame columns to
unit vectors and decides dependence from the residuals of the extra points (span 4: equal
residual directions; span 3: c′ ∈ ⟨a′, b′⟩, six directions per pair; span 2: batched rank of the
residual blocks); the two methods agree on every span.

| span j | frame sets | dependent | types (nullity) | full-support |
|---|---|---|---|---|
| 2 | 163 185 | 1 515 | 1 500 × five collinear + 2 (nullity 1), 15 × six collinear + 1 (nullity 2) | 0 |
| 3 | 6 435 689 | 60 | 60 × five collinear + 2 (two frame points and three extras on one of the C(4,2) frame lines) | 0 |
| 4 | 2 869 210 | 0 | – | 0 |

The counts are exactly the sets with ≥ 5 collinear points (span 2: 3 · C(5,3) · 41 + 3 · 6 · C(6,4)
= 1 500 and 3 · C(5,4) = 15; span 3: 6 · C(5,3) = 60; span 4: a line meets the frame in ≤ 2 points,
so 5 collinear is impossible), which is the "iff" of the lemma.

**Cross-check (`rm20_w7_sample_check.py`, 6 × 10⁵ random 7-sets, brute-force null spaces).**
Strata: uniform on the grid, inside a random plane, inside a random 3-flat, line-heavy sets
(4–7 points of a random line plus random points off it), 7-subsets of the union of two random
lines, 7-subsets of random plane conics (random quadrics, degenerate ones included). In every
stratum every dependent set contains ≥ 5 collinear points and every full-support word is a full
line; 0 violations in 600 000 sets (16 min, dominated by the classification of the dependent sets):

| stratum (10⁵ sets each) | dependent | five collinear + 2 | six collinear + 1 | full line | full-support words |
|---|---|---|---|---|---|
| uniform on F_7^4 | 0 | – | – | – | 0 |
| inside a random plane | 1 172 | 1 148 | 24 | 0 | 0 |
| inside a random 3-flat | 3 | 3 | 0 | 0 | 0 |
| 4–7 points of a line + random points | 74 894 | 24 972 | 25 029 | 24 893 | 24 893 (all full lines) |
| 7 of the union of two lines | 28 878 | 25 948 | 2 872 | 58 | 58 (all full lines) |
| 7 of a plane conic | 7 086 | 5 262 | 693 | 1 131 | 1 131 (all full lines) |

The nullities are 1, 2, 3 for the three configurations, as the lemma predicts (two, one, zero
missing points of RM(2,1) on the line).

### The mixed stage at weight 7 (`scratchpad/mixed_w7_check.py`)

For a weight-7 word e of ker G0 with non-class part e_N ≠ 0, u = (g−1)(g−6)·e_N is a word of
RM(20,4) of weight 5, 6 or 7 inside the non-class points N, hence a line word, and the class part
has 2, 1 or 0 columns. Writing s for the syndrome G0 e_N in the 21 block coordinates (the
coordinates in which the 402 class columns live after the 35 kernel rows are split off by a row
reduction of the class block; the kernel-row part of s vanishes automatically and is asserted),
the three cases are

| non-class part | class columns | condition on s | candidates (k = 270 seed 407) |
|---|---|---|---|
| weight-7 line word, all 7 points in N | 0 | s = 0 | 1 536 lines × 22 = 33 792 |
| weight-6 line word, its 6 points in N | 1 | s ∝ c_i | 1 536 × 14 + 14 112 × 2 = 49 728 |
| weight-5 line word, its 5 points in N | 2 | s ∈ ⟨c_i, c_j⟩ | 1 536 × 21 + 14 112 × 6 + 32 928 × 1 = 149 856 |

(lines with 7 / 6 / 5 points in N: 1 536 / 14 112 / 32 928; the seventh point of a line with 6 in N
may be a class column, an output or a g = 0 point, it does not enter). All three conditions are one
hash lookup: the directions of the zero vector, of the 402 single class columns and of the
6 · C(402,2) = 483 606 combinations c_i + β c_j are tabulated once as base-7 integers
(`pointsets.pair_direction_table`, 484 009 keys, no duplicates — a duplicate would be a weight-3
block word), and the normalised s of every candidate is binary-searched
(`pointsets.lookup_directions`). A hit is verified explicitly against G0 and reported with its
weight, so lower-weight mixed words (which the mixed6 stage excludes) would be caught as well.
A positive control plants the syndromes of real candidates as fake single and pair class columns
and checks that the lookup finds them. Cost: table 0.5 s, 233 376 candidates in 3 s (the weight-5/6
part reproduces the 1 197 504 = 6 × (149 856 + 49 728) non-projective candidates of the mixed6
stage); the stage is now cheaper than the plane filter. If non-line configurations had appeared in
the lemma, each type would have needed its own enumeration inside N (planar types: the ~10⁵
images per plane under AGL(2,7) over 139 650 planes, pruned by the N-mask line by line); none is
needed.

**Result for k = 270 seed 407 ([[1746,270,≥7]]₇, d_plane = 8).** 233 376 candidates, 0 extensions
(`runs/certify/p7_k270_seed407_mixed7.json`). Together with the earlier stages (no block word of
weight ≤ 6, no word of weight ≤ 5 anywhere, no mixed weight-6 word) this proves: *every weight-7
word of this code, if any, lies entirely in the class block.* d = 8 (γ = 0.898) is therefore
decided by the 3+4 weight-7 scan of the 402 class columns alone (the deferred lenore job); no
further open-ended scan is needed. The same two lookups (mixed6 + mixed7, ~10 s) apply to every
sphere-charge draw, so for any survivor the certificate of d = 8 is "class block has no word of
weight ≤ 7, full code has no word of weight ≤ 5". (The candidate set is the same for every draw,
because N = {g ∈ {2,3,4,5}} does not depend on the outputs; only the class columns, hence the
direction table, change.) Run on the certified [[1756,260,7]]₇ draw as well
(`runs/weighted/p7_m4_sphere_1756_260_mixed_w7_check.json`, block rank 31, 484 009 keys): 0 of
233 376 candidates extend, so every weight-7 word of that code is a class-block word and the
queued weight-7 census of its class block counts all of them.

## Weight-8 words of RM(20,4) and a proved weight-8 bound (2026-11-11)

*Scripts `scratchpad/rm20_w8_count.py`, `scratchpad/rm20_w8_brute.py` (cross-check), `scratchpad/w8_proved.py`;
data `runs/weighted/rm20_w8_count.json`, `runs/weighted/w8_proved.json`.*

**Why.** The rigorous curve (`qmsd/rigorous_curve.py`) uses the exact weight-7 census and the Singleton bound
N_w ≤ C(n,w) p^(w−d+1) for every heavier weight. At the input errors the protocols use (10⁻⁶ to 10⁻⁵) the
weight-8 term of that bound, about 10²³ words, is what keeps the proved savings below the estimate. Replacing
it by any bound below about 1.5 × 10²¹, with Singleton kept from weight 9 on, already gives the optimistic
saving. An exact weight-8 count is not needed, only an upper bound.

**Bound.** Let e be a weight-8 word of ker G0 for [[1746,270,7]]₇ (class block of b = 402 columns), and e_N its
part on the non-class points N. Then u = (g−1)(g−6)·e_N is a word of RM_7(20,4) with the support of e_N, so
|e_N| ∈ {0, 5, 6, 7, 8}.

- |e_N| = 0: a block word. Its 8-set carries a subcode of length 8 and distance ≥ 7, so of dimension ≤ 2. A
  [8,2,7] code is MDS and has no word of weight 8, so each 8-set carries at most p − 1 = 6 such words:
  at most 6·C(402,8) = 9.46 × 10¹⁶.
- |e_N| = 5, 6, 7: u is a line word (weight-7 lemma), and the block part has weight 3, 2, 1. Two completions
  of the same e_N differ by a block word of weight ≤ 6, which does not exist (d = 7). So there are at most
  A₅ + A₆ + A₇ (RM(20,4)) = 4.69 × 10⁷ such words.
- |e_N| = 8: u is a weight-8 word of RM(20,4). There are exactly A₈ = 436 272 744 600 of these (below).

So N₈ ≤ 9.46 × 10¹⁶ for [[1746,270,7]]₇ and 1.15 × 10¹⁷ for [[1756,260,7]]₇ (b = 412), and L₈ ≤ N₈.

**A₈ of RM_7(20,4), by frame enumeration (`rm20_w8_count.py`, 3.5 min, one core).** This is the method of the
weight-7 lemma with 8 points: frame 0, e₁..e_j plus 8 − (j+1) extras of the coordinate j-flat, with
dependence decided from the extras' residuals after the frame columns are row-reduced.

| span j | frame sets | method | dependent | full-support sets | weight-8 words |
|---|---|---|---|---|---|
| 2 | C(46,5) = 1 370 754 | batched rank | 38 382 | 1 521 | 436 272 744 600 |
| 3 | C(339,4) = 5.4 × 10⁸ | 60 dependent triples + pair-direction collisions | 23 310 | 0 | 0 |
| 4 | C(2396,3) = 2.3 × 10⁹ | pair directions against single directions | 100 | 0 | 0 |

Span 4 has no proportional residual pair (weight 7), so a set there is dependent iff c′ ∝ a′ + βb′. The 100
dependent sets are exactly the "2 frame points + 3 extras on their line" sets (10 · C(5,3)). Span 3 adds
the 4-sets containing one of the 60 dependent triples and the collisions a′ + βb′ ∝ c′ + δd′. None of
these carries a full-support annihilator. Every weight-8 word is planar, and the annihilator of its
support is unique up to scale (nullity 1):

| support | frame sets | ordered frames per set | count |
|---|---|---|---|
| 8 points of an ellipse | 21 | 336 | 139 650 planes × 6 174 ellipses × 6 = 5 173 194 600 |
| two lines, 4 + 4 points (intersection point not used) | 1 500 | 288 | 139 650 × 514 500 × 6 = 431 099 550 000 |

Orbit counting: #(8-sets of span j) = Σ over frame sets R of N_j / F(R), where N_j is the number of
ordered affinely independent (j+1)-tuples of F_7⁴ and F(R) is the number of ordered frames inside R.
Summed over all C(46,5) span-2 frame sets this gives exactly 139 650 · C(49,8) (asserted). Both rows
also agree with independent closed forms:
- **Ellipses:** a fixed line is external to q³(q−1)²/2 = 6 174 conics of PG(2,7).
- **Line pairs:** 1 372 intersecting pairs × C(6,4)² + 168 parallel pairs × C(7,4)² = 514 500 per plane.

This is the finite-field Cayley–Bacharach picture at 2d + 2 = 8 points: the 8-point sets that fail to
impose independent conditions on cubics and carry a full-support relation are conics.

**Cross-check (`rm20_w8_brute.py`).** A brute-force rank test of all 5.4 × 10⁸ span-3 frame sets is running
on loranne (4 workers).

**Effect on the proved numbers (`w8_proved.py`).** These are the curves with the census at weight 7, the bound
above at weight 8, and Singleton from weight 9:

| | [[1746,270,7]]₇ | [[1756,260,7]]₇ |
|---|---|---|
| f_hi(10⁻⁶), before → now (surrogate) | 6.2e-32 → 1.4e-35 (6.9e-46) | 6.4e-32 → 1.5e-35 (4.6e-48) |
| rigorous threshold | 0.340 % → 0.355 % | 0.338 % → 0.353 % |
| proved saving 1e-3 → 1e-30, all species | 3.09× → **3.60×** (= optimistic) | 2.96× → **3.44×** (= optimistic) |
| proved saving 1e-2 → 1e-30, raw U only | 1.00× → **2.00×** (= estimate) | 1.00× → **1.92×** (= estimate) |
| proved saving 1e-3 → 1e-30, raw U only | 1.78× → **1.90×** (estimate 2.05×) | 1.71× → **1.83×** (estimate 1.97×) |

On the default grid the proved all-species savings now equal the optimistic bound for both codes. The
proved 1e-3 protocol becomes [[6,1,?]]₇ RS → code, at 39.2 states and 171 gates per output, against
[[7,1,3]]₇ → code at 45.6 and 209. The remaining raw-U-only gap at 1e-3 sits at ε = 7.5 × 10⁻⁶, where the
weight-9 Singleton term binds; the weight-9 bound below closes most of it. The argument needs d = 7 and the class-block structure. It does not apply to
[[1746,270,6]]₇ or to the unit-charge [[2093,308,5]]₇.

### Weight 9 (2026-11-11, `scratchpad/rm20_wW_count.py 9`, `runs/weighted/rm20_w9_count_v2.json`)

The weight-9 bound does not need unique completions everywhere. With crude counts where completions are not unique,
it gives

N₉ ≤ 126·C(b,9) + A₅·6⁴·C(b,4) + A₆ + A₇ + A₈ + A₉(RM₇(20,4)).

- **Block words.** The subcode on a 9-set has dimension ≤ 3 and distance ≥ 7. The [9,3,7] MDS code has
  342 − 9·48 + 36·6 = 126 words of weight 9, and smaller dimensions have fewer.
- **Weight-5 line part plus 4 class columns.** Each completion is one of the 6⁴·C(b,4) weight-4 block vectors, so no
  uniqueness is needed.
- **6, 7 or 8 non-class points.** The completion has weight 3, 2 or 1 and is unique.
- **9 non-class points.** These are the weight-9 words of RM₇(20,4).

The bound gives N₉ ≤ 1.11 × 10²⁰ for [[1746,270,7]]₇ and 1.35 × 10²⁰ for [[1756,260,7]]₇, against the Singleton
1.4 × 10²⁶.

**A₉ of RM₇(20,4) = 9 202 538 394 000, all planar.** The counter generalises the weight-8 enumeration, and
rerun at W = 8 it reproduces A₈ exactly. In span 3 and span 4 it takes every candidate set, including those
containing a dependent smaller subset whose last point lies in its span, plus the extras that form a circuit
(collisions between pair directions and pair or triple directions). Neither span yields a full-support set
(71 397 and 29 750 candidates). Span 2 (C(46,6) sets) has 42 881 full-support frame sets in 8 types: 6 types
of nullity 1, plus 2 100 sets of nullity 2 with 12 full-support vectors and 1 200 sets of nullity 2 with 30.
Orbit check: 139 650 · C(49,9).

| with the weight-8 and weight-9 bounds | [[1746,270,7]]₇ | [[1756,260,7]]₇ |
|---|---|---|
| f_hi(10⁻⁶) | 5.9e-38 | 7.2e-38 |
| rigorous threshold | **0.370 %** | **0.368 %** |
| raw U only, 1e-3 → 1e-30 | 1.90× → **2.00×** (estimate 2.05×) | 1.83× → **1.92×** (estimate 1.97×) |
| 30 % species spread (r = 1.3), 1e-3 → 1e-30 | 3.18× → **3.58×** (= estimate) | – |

All other grid entries are unchanged, because weight 8 had already made them optimal. With up to 30 % species
spread, proved now equals the estimate at every entry of the sensitivity table. At r = 2–3 the 1e-2 entries
(1.31×, 1.08× proved against 4.04×, 3.68× estimated) would need weights ≥ 10. Data:
`runs/weighted/w89_proved.json`, `species_cost_w89.json`, `species_sensitivity_w89.json`,
`runs/coverage/total_cost_w89.json` (unchanged from weight 8).

## Orbit-union output sets (2026-10-31, in progress)

Output sets that are unions of orbits of a group Γ ≤ O(4, F₇) (so g(Mx) = g(x) and the cube
classes are preserved) are a different population from the random cube-class draws, and a
Γ-symmetric code's low-weight words come in Γ-orbits. `scratchpad/orbit_rank_probe.py` first
asks which groups allow a full-rank set of size k at all: the output columns must be independent
in RM(7,4)|_X, and a Γ-invariant set respects the isotypic decomposition of that space, so a
component that is too small caps k. Any group containing −I is capped at 2·dim(even part) = 260
(rank deficiency 17 at k = 272 for shift+neg, rot12, …); S₄-type and rotation groups are also
deficient. Full-rank unions past k = 270 exist only for small groups without −I: the coordinate
shift x₁→x₂→x₃→x₄ (order 4, full rank to k ≈ 285), the 3-cycle (order 3) and neg₁₂ (order 2).

`scratchpad/orbit_sweep.py` (same cheap filters as `plane_sweep.py`: plane distance, then the
class-block weight ≤ 5 check) on the shift group, 60 draws per k:

| k | n | population | d_plane 3 / 4 / 5 / 6 / 7 / 8 | survivors | rate |
|---|---|---|---|---|---|
| 270 | 1746 | random (600 draws) | 6 / 17 / 102 / 424 / 50 / 1 | 51 | 8% |
| **272** | 1744 | **shift orbits** | 0 / 5 / 5 / 23 / 19 / **8** | **18** | **30%** |
| 275 | 1741 | random (63 draws) | 1 / 0 / 9 / 53 / 0 / 0 | 0 | 0% |
| **275** | 1741 | **shift orbits** | 0 / 5 / 14 / 29 / 11 / 0 | 2 | 3% |
| 275 | 1741 | 3-cycle orbits | (3 draws with d_plane ≤ 2) 5 / 9 / 6 / 29 / 2 / 0 | 1 | 2% |
| 280 | 1736 | shift orbits | 0 / 4 / 13 / 42 / 0 / 0 (one d_plane 1) | 0 | 0% |
| 280 | 1736 | 3-cycle orbits | 0 / 2 / 11 / 39 / 0 / 0 (8 with d_plane ≤ 2) | 0 | 0% |
| 285 | 1731 | shift orbits | 1 / 18 / 29 / 12 / 0 / 0 | 0 | 0% |

The shift-symmetric sets are clearly the better population at k = 272 (30% survivors against 8%
for random sets at k = 270, and 8 of 60 with plane distance 8 against 1 of 600) and they keep a
few survivors at k = 275 where random draws have none; by k = 280 the window is closed for them
too. The 3-cycle is worse than random (many sets with plane distance ≤ 2: a 3-cycle-invariant set
contains fixed points of the cycle, whose orbits are single points, and the small orbits cluster).

**Certificates: all d = 6.** Every one of the 18 shift-symmetric k = 272 survivors (including the
8 with plane distance 8) and both k = 275 survivors has a weight-6 word inside the class block
(`certify_candidate.py --file`, block stage, 6 min each on lenore). The killer words are non-planar
(their six points span a 3-flat or all of F₇⁴, so the plane filter cannot see them) and mostly not
symmetric (four of them meet six different shift-orbits). The 15 survivors of the order-2 group
neg₁₂ at k = 272 are in the same chain.

**Why, and where the window really ends.** The class-block lemma gives the block 291 − k
independent checks on its 402 columns (the 402 is fixed by the charge: 672 class points minus
~270 outputs). The accidental-word count E[projective weight-w words] ≈ C(402, w)(p−1)^{w−1}/p^{291−k}
then reads

| k | checks | E[weight-5] | E[weight-6] | P(no weight-6 block word) | E[weight-7] |
|---|---|---|---|---|---|
| 269 | 22 | 3e-5 | 0.011 | 0.99 | 3.8 |
| 270 | 21 | 2e-4 | 0.079 | 0.92 | 27 |
| 271 | 20 | 1.4e-3 | 0.55 | 0.58 | 187 |
| 272 | 19 | 1e-2 | 3.9 | 0.02 | 1300 |
| 273 | 18 | 0.07 | 27 | 2e-12 | 9200 |
| 275 | 16 | 3.3 | 1300 | 0 | 4.5e5 |

This is the window edge: at k ≤ 270 a draw that passes the plane filter almost always has
d = 7 (both k = 270 survivors we certified did; the seed-19 census is finding ~30 weight-7 classes,
matching the 27 predicted), at k = 271 about half do, at k = 272 one in fifty, and from k = 273 none.
The 0 of 18 (shift) at k = 272 is what the table says for existence. **Correction (2026-11-01):** a
sampled census of one shift-symmetric k = 272 block (`scratchpad/routes/isd/`, Stern 2+2, 1552
iterations) finds about 115 distinct weight-6 words against the law's 3.7, concentrated on four
columns that very likely form one shift orbit. So symmetry does change the accidental count, and
upward: orbit-union sets are worse than random draws at the block stage, not merely no better.
The law matches random blocks (sampled N_7 = 32 ± 6 on the k = 270 seed-19 block, law 26.7, exhaustive
census about 30) and exact counts at p = 11 (below).

**k = 271 measured (2026-11-01).** 160 random draws give 13 plane-filter survivors (8%); the block
stages on lenore certified all thirteen: **9 are d = 7 exact** and 4 have a weight-6 block word:
[[1745,271,7]]₇, γ = **0.957**, the new p = 7 record (catalog `imported-p7-n1745-k271-*`,
`runs/certify/p7_k271_bundle.json`). Measured pass rate 9/13 = 69%, against the law's 58%.

**The law at p = 11 (exact).** On restrictions of the p = 11, m = 4 candidate to 13 or 14 checks and
100 or 120 columns, where an exhaustive scan fits, the mean number of weight-6 words is 3.67 (law
3.45, 12 subsets) and 1.20 (law 0.96, 10 subsets); the fraction with none is 0.30 against the law's
0.38. Data: `scratchpad/routes/isd/runs/law_restricted_*.json`.

What symmetric sets *would* buy is a cheaper certificate (scans over orbit representatives only),
which does not matter once the window is this narrow. Data: `runs/weighted/orbit_sweep/`,
`runs/certify/p7_k272_*`, `runs/certify/p7_k275_*` on lenore.

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
  classes). Their costs hold only if raw states of all three species are available at the same
  ε₀; see "Species accounting" below.
- d of a draw varies strongly at fixed k; the certified codes are single draws, not optima.

## Species accounting (2026-11-04)

At p = 7, U^a and U^b are Clifford-equivalent only if a/b = ±1, so there are three species: U (μ = 1, 6),
U² (2, 5) and U³ (3, 4). A sphere-charge code with k outputs in the U class takes 672 − k U inputs, 672 U²
inputs and 672 U³ inputs. The unit-charge [[2093,308,5]]₇ takes only U. The library [[20,1,4]]₇ also
mixes species (10 U, 2 U², 8 U³ relative to its output class). All other library codes are unit-charge.

**Rescaling.** Multiplying all charges by c leaves G0 and H1 unchanged and multiplies the logical exponent
by c, so the same code with the same curve f(ε) makes U^c from rescaled inputs. Checked on the
[[20,1,4]]₇ triple sums at c = 1, 2, 3. Each species therefore has an identical distillation tower.

**Conversion.** Teleporting U^b onto |U^a⟩ gives |U^(a+b)⟩: map |s⟩|x⟩ → |s⟩|s−x⟩, measure y = s−x, and
correct with the Clifford phase ω^(3b s²y − 3b s y²). Z errors add, so errors compose as
ε₁ + ε₂ − ε₁ε₂ p/(p−1), and the result is again uniform on nonzero Z powers. The state-vector check over
all a, b, y and error powers has worst infidelity 7e-16 (`tests/test_species.py`). From raw U alone,
|U²⟩ costs 2 U at about 2ε and |U³⟩ costs 3 U at about 3ε. One sphere-charge round then costs
(672 − k) + 2·672 + 3·672 U states: 3762 at k = 270, against n = 1746.

**Two scenarios**, from `scratchpad/species_cost.py` (`runs/weighted/species_cost.json`). It is a forward
Pareto search over (species, ε, cost) that allows conversion at any tier, and every round takes all inputs
at the worst input error. Saving is cost(library) / cost(library + code) within the same scenario.

| Code | ε₀ → target | All three species raw (proved / estimate) | Raw U only (proved / estimate) |
|---|---|---|---|
| [[1746,270,7]]₇ | 1e-2 → 1e-30 | 4.32× / 4.32× | 1.00× → **2.00×** / 2.00× |
| [[1746,270,7]]₇ | 1e-3 → 1e-30 | 3.09× → **3.60×** / 3.60× | 1.78× → 1.90× → **2.00×** / 2.05× |
| [[1746,270,7]]₇ | 1e-2 → 1e-15 | 1.08× / 1.08× | 1.00× / 1.00× |
| [[1756,260,7]]₇ | 1e-2 → 1e-30 | 4.14× / 4.14× | 1.00× → **1.92×** / 1.92× |
| [[1756,260,7]]₇ | 1e-3 → 1e-30 | 2.96× → **3.44×** / 3.44× | 1.71× → 1.83× → **1.92×** / 1.97× |
| [[2093,308,5]]₇ | 1e-2 → 1e-30 | 1.27× / 4.11× | 1.00× / 4.11× |
| [[2093,308,5]]₇ | 1e-3 → 1e-30 | 2.94× / 2.94× | 3.65× / 3.65× |
| [[2093,308,5]]₇ | 1e-2 → 1e-15 | 1.03× / 1.03× | 1.03× / 1.03× |

No code saves anything at 1e-3 → 1e-15 in either scenario. "Proved" uses f_hi, pacc_lo and "estimate"
uses the leading-order surrogate. The library baselines are 52.6, 210.3 and 28.8 in both scenarios, and at
1e-3 → 1e-30 they are 141.0 (all species) and 175.1 (raw U only), because [[20,1,4]]₇ stops being worthwhile.

- **All species raw:** the search reproduces `best_protocol` and every earlier saving exactly. The old
  "one raw state each" count is correct at every tier by rescaling, so the proved 4.32× and 3.09× stand. The
  assumption is physical: U² and U³ can be prepared at the same ε₀ and cost as U. Charging w raw U-equivalents for
  each U² or U³ raw state, the [[1746,270,7]]₇ saving reaches 1× at w ≈ 5.3 (1e-2 → 1e-30) and w ≈ 4.7
  (1e-3 → 1e-30), ignoring any change in ε.
- **Raw U only:** the sphere-charge codes keep about 2×, and γ loses its meaning as a raw-state exponent.
  Counted in U states, [[1746,270,7]]₇ has log(3762/270)/log 7 = 1.35, so γ < 1 at p = 7 depends on the
  all-species convention. The unit-charge [[2093,308,5]]₇ (γ 1.19 either way) becomes the best p = 7 code,
  with a proved 3.65× at 1e-3 → 1e-30.
- **The proved 1.00× at 1e-2 → 1e-30** was a limit of the bound. At the ε = 3.0e-6 this route needs, f_hi = 4e-28
  (the Singleton tail) missed the 1e-30 target, while the leading term is 2e-42. The weight-8 bound of 2026-11-11
  ("Weight-8 words of RM(20,4)" above) proves the 2.00× (arrows in the table; `runs/weighted/species_cost_w8.json`).

**Not done.** Per-column input errors are not modelled, though the low-weight words all lie in the U class
block, so U² and U³ columns could tolerate more noise. Conversions by a distillation code with
charge-changing outputs were not looked at. Clifford cost under conversion is not recomputed: the
2016 teleports per k = 270 round add about 7.5 per output, against 98 encoder gates per output.

### Testing the all-species assumption (2026-11-04)

Scripts `scratchpad/species_verify.py`, `species_snap.py`, `species_sensitivity.py`; data
`runs/weighted/species_verify.json`, `species_snap.json`, `species_sensitivity.json`.

| Test | Result | Status |
|---|---|---|
| Exact state-vector round of a mixed-species code (g = x on F₇, charges 2-6, output in the U class, n = 5), rescaled by 1, 2, 3, ε = 1e-3 to 5e-2 | P_acc and P_fail match the G0, H1 prediction to about 1e-12 | confirmed on this code; the rescaling itself is an identity |
| Pauli noise after the state (depolarizing, X, XZ², Z), and n-dephasing or coherent diagonal errors | same ε for all a, spread ≤ 1e-12 | confirmed |
| Loss before an ideal diagonal gate | ε = 6.053e-4 for all a | confirmed; any noise before an ideal diagonal gate is species-blind |
| Loss after the gate, coherent X rotation after the gate | ε spreads 30 % across species (loss g = 1e-3: U 2.76e-3, U² 2.86e-3, U³ 2.01e-3) | species-dependent |
| Pulse-level SNAP on Fock 0-6, unoptimized Gaussian pulses, phases calibrated per target, χ/2π = 1 MHz, cavity T1 1 ms, transmon T1 = Tφ = 100 µs | coherent ε differs by 2.3× (3.9e-3 to 9.0e-3 at σχ = 16); with decoherence ε is 0.20-0.21 at σχ = 16 (spread 4 %), 0.10-0.14 at σχ = 6 (spread 25 %) | species-dependent, small once decoherence dominates |
| Mana | U and U² 0.815, U³ 0.896 | conversion lower bound 1 U and 1.10 U, so the teleport's 2 and 3 are not shown optimal |

The SNAP ε values are far above any ε₀ we use, because the pulses are unoptimized and the noise values are
illustrative. Only the species dependence is the result, not the absolute ε. Landgraf et al. (arXiv:2310.10498) also report
target-dependent SNAP coherent errors. Wigner positivity is reached at ε = 0.49, 0.70, 0.54 for U, U², U³, so the
species do differ at large noise, far above any code threshold used here.

**Hardware.** On a transmon qudit, diagonal gates can be software phase updates: Champion et al. (arXiv:2405.15857, d = 3, 5, 8)
implement SNAP layers "by updating subsequent drive phases, producing virtual rather than physical phase rotations".
That paper does not measure a SNAP fidelity ("close to unity" is a statement, not a benchmark) and does not address
dependence on the phase values or d = 7. Its native displacement fidelity at d = 8 is 0.9895. With a virtual gate the
state is prepared by physical pulses followed by a frame update, which is the species-blind ordering in the table.
Nguyen et al. (arXiv:2312.17741, d ≤ 4) also use virtual Z gates. On a cavity, SNAP is a physical pulse whose duration is
set by χ rather than the phases, and the phases enter only as tone phases, so the cost is the same but the error is not.

**How much it matters** ([[1746,270,7]]₇, raw U at ε₀ cost 1, raw U² and U³ at r·ε₀ cost w; proved / estimate):

| | 1e-2 → 1e-30 | 1e-3 → 1e-30 |
|---|---|---|
| r = 1, w = 1 | 4.32 / 4.32 | 3.60 / 3.60 |
| r = 1.3, w = 1 | 4.24 / 4.24 | 3.18 → **3.58** / 3.58 |
| r = 3, w = 1 | 1.08 / 3.68 | 3.17 / 3.46 |
| r = 1, w = 2 | 2.44 / 2.44 | 2.52 / 2.52 |
| r = 1, w = 3 | 1.70 / 1.70 | 1.76 / 1.76 |
| r = 1, w = 5 | 1.06 / 1.06 | 1.09 / 1.09 |

(Proved column with the weight-8 bound of 2026-11-11, `runs/weighted/species_sensitivity_w8.json`; the arrow is the weight-9 bound,
`species_sensitivity_w89.json`, which changes no other row. Before it the proved
entries were 4.32, 1.49, 1.08, 2.44, 1.70, 1.06 at 1e-2 and 3.09, 3.09, 3.08, 2.17, 1.51, 1.00 at 1e-3.)

So error differences up to 3× barely move the estimates. A 30 % difference leaves the proved 1e-2 saving at 4.24×
(1.49× before the weight-8 bound) and lowers the proved 1e-3 saving to 3.18×; at 3× the proved 1e-2 figure falls to 1.08×,
because the proof tail is still steep from weight 9 on. Raw cost is what matters. Comparing like with like at w = 1, [[1746,270,7]]₇ leads the
unit-charge [[2093,308,5]]₇ by 18 % proved at 1e-3 → 1e-30 (39.2 vs 47.9; 5 % before the weight-8 bound) and by 5 % and 18 % estimated (48.7 vs 51.1 at
1e-2, 39.2 vs 47.9 at 1e-3), and falls behind it at w ≈ 1.1-1.3. Only the proved comparison at 1e-2 is large (48.7 vs 166.2),
because the d = 5 proof tail is weak. The species assumption buys proved γ < 1 and that proved figure, not much estimated cost.

**Still unverified.** (1) Raw states from anything other than a direct diagonal gate: no qudit cultivation is known to us,
and its checks could depend on the species. (2) The measured fidelity of a virtual diagonal gate at d = 7. (3) Uniform Z
powers: non-Pauli noise splits the orbits {1,2,4} and {3,5,6} by up to about 30 % (loss g = 1e-3 at a = 3: 0.83e-3 vs 1.18e-3), and the Clifford group
cannot twirl across them. Every curve here assumes uniformity; the cost effect is not computed and affects all species.
