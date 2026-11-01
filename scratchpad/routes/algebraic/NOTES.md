# Algebraic certificates for punctured-RM magic-state codes: p = 11, m = 4, r = 13

*2026-10-31. Scripts and logs in this directory (`scratchpad/routes/algebraic/`). All
computations ran locally in minutes under `ulimit -v 4000000`, `OMP_NUM_THREADS=1`.*

## 0. Summary

1. **Random output sets (the k = 2290 target).** The only dual words that the occupancy bounds
   exclude rigorously are the *planar* ones: every word of RM_11(26,4) supported in an affine
   plane has weight >= 55 and meets O in <= 43 points, so it leaves a kernel word of weight
   >= 12. The literature classifies the words of weight <= W_3 (the third weight, 61 <= W_3 <= 63)
   and all of them are planar, so **W\* = W_3 >= 61** (W\* = 63 if W_3 = 63, which is the
   hyperplane-arrangement value; the exact W_3 at q = 11 is open in Leducq's paper). Nothing
   heavier can be excluded from occupancy alone: a hyperplane carries up to 262 outputs, so
   3-flat words of weight 64..268 and all 4-span words are untouched. Those are the
   "accidental" population and must go to brute force (infeasible at weight 6, n = 12351).
   Structurally, *no* algebraic argument can help at k = 2290: the check space has dimension
   30 < 35 = dim RM(3,4), so O lies on no hypersurface of degree <= 10.
2. **Structured O (the research question) — solved by a one-line theorem.** Take O to be the
   *full* level set {g = c} of a polynomial g of degree delta. Then g·RM(13 − delta, 4) lies in
   the check space, so for every kernel word e the vector (g(x) − c)·e(x) is a word of
   RM(13 − delta,4)^perp = RM(26 + delta, 4), hence **d >= d(RM(26 + delta, 4))** with no
   brute force at all. Validated exactly by MITM at p = 7, m = 3 and p = 5, m = 4
   (`validate_small.log`: d equals the bound in every full-level-set trial, and
   d = min(bound, d_block) in every partial-level-set trial).
3. **Construction (d).** delta = 3: d >= d(RM(29,4)) = 22, k <= dim RM(13,4) − dim RM(10,4) =
   2320 − 1001 = 1319, and by Ax's theorem cubic level sets have size ≡ 0 (mod 11), so the
   largest admissible size is 1309. About one level in eleven of a random cubic has exactly 1309
   points (12 of 132 levels sampled), and all 12 such sets imposed independent conditions on
   RM(13,4). Result: **[[13332, 1309, >= 22]]_11, gamma <= ln(13332/1309)/ln 22 = 0.7508**,
   fully proven (theorem + a rank computation done two independent ways), against gamma = 0.9405
   for the unproven d = 6 of the k = 2290 draw and 0.959 for the best certified p = 7 code.
   Explicit g and O: `p11_cubic_best.json`; verification log: `verify_cubic_code.log`.

## 1. Setting and the exact duality

G = RM_11(13,4), dim 2320, on the 14641 points of F_11^4 (r = 13 = r_max since 3·13 < 40).
O = output set, |O| = k, required to have rank k in G. G0 = the words of RM(13,4) vanishing on
O, restricted to N = F^4 \ O; dim G0 = 2320 − k. The distance is d = min weight of ker G0.

**Check space.** Let I(O)_{<=13} = {h : deg h <= 13, h|_O = 0}. Then

    e ∈ ker G0  <=>  Σ_{x∈N} e_x h(x) = 0 for all h ∈ I(O)_{<=13}.                     (1)

Equivalently (the form in the task statement) e is the restriction to N of a word u of
RM(13,4)^perp = RM(26,4) with supp u ⊆ O ∪ supp e; a kernel word of weight w is a dual word u of
weight W = w + |supp u ∩ O|. Equivalently again, in interpolation language: d >= w + 1 iff for
every T ⊆ N with |T| <= w the set O ∪ T imposes independent conditions on polynomials of degree
<= 13 (a finite-field Cayley–Bacharach condition); d <= w iff some O ∪ T fails.

Two populations of low-weight kernel words: (i) geometric — u is a low-weight dual word (a
hyperplane-arrangement word) whose support is almost covered by O; (ii) accidental — u is a
heavy dual word with almost all its support in O; equivalently, a random-looking dependency
among the 30 = 2320 − 2290 check columns.

## 2. (a) Weight hierarchy of RM_11(26,4) and of its restrictions

Write r~ = 26 = t(q − 1) + s with q = 11: t = 2, s = 6. Two different restrictions matter.
The *restriction* of RM(26,4) to a 3-flat is RM_11(26,3) (degree <= 26 in 3 variables; this has
codimension dim RM(3,3) = 20 in the space of all functions on the flat, i.e. u|_H is orthogonal
to the cubics on H — the only local constraint a hyperplane imposes), and its restriction to a
plane or a line is unconstrained (every function on F_11^2 has degree <= 20 <= 26). The *words
supported inside* a flat are more constrained: u = 1_H · v with deg 1_H = 10 forces
v ∈ RM_11(16,3) (t = 1, s = 6) for a 3-flat, v ∈ RM_11(6,2) (t = 0, s = 6) for a plane, and
v ∈ RM_11(−4, 1) = 0 for a line: **no dual word lives on a line** (because r = 13 >= p − 1 = 10;
compare "Why lines kill p = 11 and 13" in WEIGHTED_SUPPORTS.md, which is the m = 3 failure of
this condition).

| code | t, s | W_1 (min) | W_2 | W_3 | first non-planar arrangement weight |
|---|---|---|---|---|---|
| RM_11(26,4) (words of F^4) | 2, 6 | 55 | 60 | c_6 ∈ [61, 63] | 88 (span 3), 132 (span 4) |
| RM_11(16,3) (words inside a 3-flat) | 1, 6 | 55 | 60 | c_6 | 88 |
| RM_11(6,2) (words inside a plane) | 0, 6 | 55 | 60 | c_6 | — |

* **W_1 = 55 = (q − s) q^{m−t−1}** (Kasami–Lin–Peterson 1968 / Delsarte–Goethals–MacWilliams
  1970, Theorem 1.1 in [L13]): the minimum-weight words are, up to scalar, the indicator of the
  union of q − s = 5 parallel affine subspaces of codimension t + 1 = 3 (lines) inside an affine
  subspace of codimension t = 2 (a plane): **5 parallel lines in a plane**.
* **W_2 = 60 = (q − s + 1)(q − 1) q^{m−t−2}** ([L13] Theorem 1.2 case 3, valid for m >= 3,
  q >= 3, 1 <= t <= m − 2, 2 <= s <= q − 2; proved by Rolland 2010 for s ≠ 1, Geil 2008 for
  t = 0, Cherdieu–Rolland 1996 / Sboui 2007 for the codeword structure when t = 0). The
  second-weight codewords ([L13] Theorem 2.3, q >= 4, 0 <= t <= m − 2, 2 <= s <= q − 2) are
  supported in a codimension-t subspace (a plane) and are either 6 parallel lines minus their
  intersection with a transversal line (6·11 − 6 = 60) or 6 concurrent lines minus the common
  point (60). **Planar.**
* **W_3 = c_6 · q^{m−t−2} = c_6**, where c_b is the third weight of RM_q(b, 2) ([L15]
  Theorem 4.10: m >= 2, q >= 5, 0 <= a ≤ m − 2, 2 <= b <= q − 2, b ≠ 3, and c_b < (q − b + 1) q;
  here a = t = 2, b = s = 6). Theorem 3.1 of [L15] gives c_b <= (q − 2)(q − b + 2) = 63 for
  q >= 5, 4 <= b <= floor(q/2 + 2) = 7, attained by 4 parallel lines plus 2 further parallel lines
  in another direction (or 2 lines through a point of one of the 4), so 61 <= W_3 <= 63; the
  equality c_b = (q − b + 2)(q − 2) is proved in [L15] Theorem 5.1 only for q >= 16, so the exact
  value at q = 11 is open there. The third-weight codewords ([L15] Theorem 6.6, q >= 7,
  4 <= b <= q − 2, c_b < (q − b + 1) q − 1 = 65, satisfied) are (1 − x_1^{10})(1 − x_2^{10})
  g(x_3, x_4) with |g| = c_6: **planar.**
* Beyond W_3 nothing is classified. Among hyperplane-arrangement words (k blocks of d_i
  parallel hyperplanes in independent directions, weight q^{m−k} Π (q − d_i), Σ d_i <= 26,
  [L15] Appendix Theorem A.1) the lightest word not contained in a plane has 2 blocks of 9 and 7
  (weight 11·2·4 = 88, support = 8 parallel lines spanning a 3-flat) and the lightest word with
  4-dimensional span has blocks 9, 9, 8 (weight 2·2·3·11 = 132). I do not know a theorem that
  these are the minima over *all* words of the given span; treat 88 and 132 as conjectural.

References. [DGM70] Delsarte, Goethals, MacWilliams, "On generalized Reed–Muller codes and
their relatives", Inform. Control 16 (1970). [KLP68] Kasami, Lin, Peterson, IEEE IT 14 (1968).
[CR96] Cherdieu, Rolland, "On the number of points of some hypersurfaces in F_q^n", FFA 2
(1996). [S07] Sboui, "Second highest number of points of hypersurfaces in F_q^n", FFA 13 (2007).
[G08] Geil, "On the second weight of generalized Reed–Muller codes", DCC 48 (2008). [R10]
Rolland, "The second weight of generalized Reed–Muller codes in most cases", Cryptogr. Commun. 2
(2010), arXiv:0902.0058. [L13] Leducq, "Second weight codewords of generalized Reed–Muller
codes", Cryptogr. Commun. 5 (2013) 241–276, arXiv:1203.5244. [L15] Leducq, "On the third weight
of generalized Reed–Muller codes", Discrete Math. 338 (2015) 1515–1535, arXiv:1401.5234. [Ax64]
Ax, "Zeroes of polynomials over finite fields", Amer. J. Math. 86 (1964) (the number of zeros of a
degree-d polynomial in m variables over F_q is divisible by q^{ceil(m/d) − 1}).

## 3. (b) What the occupancy bounds exclude for the k = 2290 draw

Measured on `runs/rm_m4/p11_r13_k2290_seed1.json` (`occupancy_p11.py`, 30 s):

| flat | number | mean outputs | max outputs | tail |
|---|---|---|---|---|
| lines | 1 948 760 | 1.72 | **9** | 9: 4 lines, 8: 72, 7: 747, 6: 5526 |
| planes | 1 963 346 | 18.9 | **43** | 43: 1, 41: 1, 40: 4, 39: 1, 38: 16 |
| hyperplanes | 16 104 | 208.2 | **262** | 262, 258, 257, 256, 254, 253, 252, 251 (one each) |

**Exact condition.** A dual word u of weight W yields a kernel word of weight
W − |supp u ∩ O|, so it is dangerous for d >= 7 iff |supp u ∩ O| >= W − 6. For a word supported
in a flat F this needs |O ∩ F| >= W − 6.

**Rigorous exclusions.**

1. Every dual word supported in a plane: W >= 55 and |supp ∩ O| <= 43, kernel weight >= 12.
   This covers *all* planar words of every weight, not just the first classes.
2. By §2, every dual word of weight <= W_3 (>= 61) is planar. Hence **W\* = W_3 >= 61**, i.e.
   all dual words of weight <= 61 are excluded, with margin 12 at weight 55 and >= 18 at weight 61;
   if c_6 = 63 then W\* = 63.
3. Words supported in a 3-flat with W >= 269: kernel weight >= W − 262 >= 7. Words in a 3-flat
   with 64 <= W <= 268 (e.g. the 88-weight arrangement of 8 parallel lines) are **not** covered
   by the hyperplane number alone; for the specific arrangement families one can do better (8
   parallel lines carry <= 8·9 = 72 outputs, kernel weight >= 16), but there is no classification
   to make this exhaustive.
4. Words of 4-dimensional span: only |supp ∩ O| <= k = 2290 is available, so W >= 2297 would be
   needed. Nothing is excluded here.

So the occupancy route certifies exactly the geometric population (planar words, provably all
words of weight <= 61) and leaves everything else — in practice everything of weight above ~63 —
to brute force. The remaining population is what the random-code model counts:
E[#projective kernel words of weight w] = C(12351, w)·10^{w−1}/11^{30}:

| w | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| E[words] | 1.8e-18 | 5.6e-14 | 1.4e-9 | 2.8e-5 | 0.50 | 7.7e3 | 1.1e8 |

so for this draw d = 7 or 8 is expected (d >= 7 with probability ~0.6 if weight-7 words are
Poisson with the mean above), but a certificate of d >= 6 or 7 needs the MITM scans that are
infeasible at this length. **Can a Nullstellensatz / footprint / Schwartz–Zippel argument
replace them for this O?** No: the check space I(O)_{<=13} has dimension 30 < 35 = dim RM(3,4)
(`random_model_p11.log` confirms rank 2290 on O), so no polynomial of degree <= 10 vanishes on O,
i.e. O is not contained in any algebraic hypersurface of degree <= 10 and the 30 checks have no
factorisable structure. Any structured argument must lower k.

## 4. (c) The algebraic route: output sets on a hypersurface

### 4.1 Theorem (full level set)

Let g be a polynomial on F_p^m of degree delta <= r and c ∈ F_p, and let O = {x : g(x) = c}.
Assume RM(r,m) has rank |O| on O. Then every nonzero e ∈ ker G0 satisfies

    |supp e| >= d( RM_p( r~ + delta, m ) ),     r~ = m(p − 1) − r − 1.

*Proof.* For every q ∈ RM(r − delta, m) the function h = (g − c)·q is a polynomial of degree
<= r (the reduced form of a product has degree at most the sum of the degrees) vanishing on O,
so h ∈ I(O)_{<=r} and by (1) Σ_{x∈N} e_x (g(x) − c) q(x) = 0. Hence the vector
e' = ((g(x) − c) e_x)_{x∈N}, extended by 0 on O, is orthogonal to RM(r − delta, m), i.e.
e' ∈ RM(r − delta, m)^perp = RM(m(p−1) − (r − delta) − 1, m) = RM(r~ + delta, m). Since g − c
does not vanish on N, supp e' = supp e, and e ≠ 0 gives |supp e| >= d(RM(r~ + delta, m)). ∎

The bound is the minimum distance of a *larger* RM code than the dual RM(r~, m), which is why it
beats the generic d(RM(r~,m)) − occupancy arithmetic: the output set is removed "for free"
because the checks factor through g.

### 4.2 Corollary (partial level set, class block)

If O ⊊ Z := {g = c} has rank |O| and B = Z \ O, then every kernel word meeting N' = F^m \ Z has
weight >= d(RM(r~ + delta, m)) (same proof, restricted to the N' part), and the remaining kernel
words are supported in B, where they are the dependencies of the residual check matrix (the
h_Z(r) − |O| checks of I(O)/I(Z) restricted to B, h_Z(r) = rank of RM(r,m) on Z). So
d = min( d(RM(r~ + delta, m)), d_block ). This *is* the class-block lemma of
WEIGHTED_SUPPORTS.md: at p = 7, m = 4, r = 7 the sphere's cube class is Z of the quartic
(g − 1)(g − 6), h_Z(7) = 326 − 35 = 291 (the cap below), and the mixed-word bound is
d(RM(16 + 4, 4)) = d(RM(20,4)) = 5 — exactly "point 3" of that lemma. (That code is weighted,
but the output-side statement is the unweighted one and matches.)

### 4.3 The cap on k and the design table at p = 11, m = 4, r = 13

Because (g − c)·RM(r − delta, m) ⊆ I(Z)_{<=r} and multiplication by g − c is injective on
RM(r − delta, m) unless a degree-(r − delta) polynomial vanishes on all of N (impossible when
|N| > (r − delta) p^{m−1}), the rank of RM(r,m) on Z is at most
dim RM(r,m) − dim RM(r − delta, m). The full-level-set construction therefore needs
|Z| <= cap(delta) := dim RM(13,4) − dim RM(13 − delta, 4), and then the question is only whether
the evaluation of the 2320 − dim RM(13 − delta, 4) dimensional quotient on Z has full rank |Z|
(generically yes: it is a random-looking matrix with at least as many rows as columns).

| deg g = delta | cap on k | d >= d(RM(26 + delta, 4)) | gamma at k = cap | typical \|Z\| and status |
|---|---|---|---|---|
| 1 | 520 (hyperplane: 530 = dim RM(13,3)) | 44 | — | \|Z\| = 1331 ≫ cap: block of >= 800 columns, dead |
| 2 | 959 | 33 | 0.760 | nondegenerate quadrics have 1210–1441 points ≫ cap: block, dead |
| **3** | **1319** | **22** | **0.748** | \|Z\| ≈ 1331 ± a few hundred, ≡ 0 mod 11 (Ax): **1309 realised, independent** |
| 4 | 1605 | 11 | 0.874 | \|Z\| ≈ 1331 < cap: generic quartic level sets independent, gamma ≈ 0.96 at k ≈ 1331 |
| 5 | 1825 | 10 | 0.846 | k ≈ 1331 → gamma 1.00 |
| 6 | 1990 | 9 | 0.842 | k ≈ 1331 → gamma 1.05 |
| 7 | 2110 | 8 | 0.857 | k ≈ 1331 → gamma 1.11 |
| 8 | 2194 | 7 | 0.892 | needs a degree-8 hypersurface with ~2194 independent points |
| 9 | 2250 | 6 | 0.952 | — |
| (random O, k = 2290) | — | — | 0.9405 if d = 6 | not provable |

(`random_model_p11.log` prints the table; d(RM(26 + delta, 4)) = 15 − delta for 4 <= delta <= 9
and 22, 33, 44, 55 for delta = 3, 2, 1, 0.) The cubic is the sweet spot: delta = 3 is the last
degree at which 26 + delta < 30 = 3(p − 1), so the bound jumps from the "p − b" regime (<= 11)
to the "2·11" regime. For gamma, k = 1309 with d = 22 beats every other row even at its cap:
delta = 4 would need k > 2080 > 1605 to reach gamma 0.75.

Why the usual tools do not apply directly to a random O, and what structure makes them work:
a Nullstellensatz/footprint argument needs a generating set of I(O) with controlled degrees;
for O = {g = c} the degree-<= 13 part of I(O) contains the whole module (g − c)·RM(10,4) and the
theorem is exactly "footprint of a principal ideal + duality". Schwartz–Zippel enters through
d(RM(29,4)) = 22 (a degree-<= 10 polynomial cannot vanish on all but <= 21 points of F_11^4 …
in its dual form: <= 21 points impose independent conditions on degree-<= 10 polynomials).
More general structures: O = Z(g_1, …, g_s) gives checks Σ (g_i − c_i) RM(13 − delta_i, 4) and a
kernel word must satisfy *all* (g_i − c_i)e ∈ RM(26 + delta_i, 4); but the cap drops by the
dimension of that sum (two quadrics: cap < 0), so complete intersections of codimension >= 2
are useless at this size. Unions Z(g_1) ∪ Z(g_2) = Z(g_1 g_2) are just the degree delta_1 + delta_2
row of the table. Within this family the cubic level set is optimal at p = 11, m = 4.

### 4.4 Ax's theorem fixes the attainable k

The number of zeros of a degree-3 polynomial in 4 variables over F_11 is divisible by
11^{ceil(4/3) − 1} = 11 [Ax64], so cubic level sets have sizes …, 1298, 1309, 1320, 1331, … and
the largest one below the cap 1319 is **1309**. Observed level-set sizes over the cubics sampled
(`p11_cubic_hunt.log`, 12 random cubics = 132 level sets): sizes from 1232 to 1408, mean
exactly 1331, all ≡ 0 mod 11; 34% of the level sets are <= 1319 and 9.1% (12 of 132, i.e. about
one level per cubic) are exactly 1309. **All 12 level sets of size 1309 were independent** (rank
1309 on RM(13,4), 20–60 s each with galois); the first random cubic of
`p11_hypersurface_scan.py` (|Z| = 1298) was independent as well, and the sets with
|Z| = 1331, 1353 > cap had rank exactly 1319 = cap, i.e. the quotient evaluation has full rank in
every case seen. Quartics behave the same way (`p11_quartic_scan.log`): three random quartic
zero sets of sizes 1373, 1359, 1327 (all below the cap 1605) were all independent, giving
provable [[13268, 1373, >= 11]]_11 (gamma <= 0.946) etc., inferior to the cubic as the table says.

## 5. (d) The construction: [[13332, 1309, >= 22]]_11, gamma <= 0.7508

**Code.** p = 11, m = 4, r = 13 (standard unweighted punctured RM_11(13,4), one species of
magic state). g = the cubic with coefficients `coeffs` (in `monomial_exponents(4, 3, 11)`
order) and O = {g = level} from `p11_cubic_best.json` (1-based columns in
`reedmuller.points` order), |O| = 1309, n = 14641 − 1309 = 13332.

**Proof of d >= 22.** By §4.1 with delta = 3: every kernel word e has (g − c)e ∈ RM(29,4), and
d(RM_11(29,4)) = (11 − 9)·11^{4−2−1} = 22 (29 = 2·10 + 9) [DGM70]. The hypothesis "rank 1309 on
O" was verified twice: galois `matrix_rank` in the hunt and an independent numpy Gaussian
elimination in `verify_cubic_code.py`. Nothing else is needed; in particular no MITM.

**gamma** = ln(13332/1309)/ln 22 = 2.3209/3.0910 = **0.7508** (upper bound: the true d may exceed
22, since only 1001 of the 1011 checks were used in the proof; a weight-22 kernel word would have
to be two full parallel lines avoiding Z(g) on which the remaining checks also vanish).

**Verification log** (`verify_cubic_code.log`, 7 min): (1) the stored O equals {g = 3} for the
stored cubic, whose degree-3 part is nonzero; (2) independent numpy elimination mod 11: rank of
RM(13,4) on O = 1309 = |O|; (3) G0 built as the left null space of G[:, O] times G, restricted to
N: shape 1011 × 13332, and G0 vanishes on O; smoke tests of the theorem on the actual G0: 0 of 200
random 21-column subsets dependent, 0 of 100 sets of 21 plane points dependent, and all 20
tested pairs of parallel lines lying entirely in N (22 columns) have rank 22 (so no weight-22
word of the two-parallel-lines type among those); (4) no kernel word of weight <= 3 (MITM on an
18-row random projection of G0, hits to be verified against G0 — none occurred); (5) gamma
0.7508 and d(RM(29,4)) = 22 recomputed.

**Comparison.** Target draw [[12351, 2290, 6?]]_11: gamma 0.9405 if d = 6, 0.866 if d = 7,
neither provable today. Best certified p = 7 code: [[1746, 270, 7]]_7, gamma 0.959. The m = 2
ladder: [[296,65,6]]_19 0.846, [[424,105,6]]_23 0.779, [[666,175,6]]_29 0.746; Reed–Solomon
gamma 0.924 at p = 23. So the cubic level-set code is the lowest certified gamma at any p <= 23
in the repository, with a proof that costs one rank computation. Its drawbacks are the size
(n = 13332) and the output rate k/n = 0.098.

**Other primes (same theorem, not run).** p = 7, m = 4, r = 7: cubic cap 326 − 70 = 256, bound
d(RM(19,4)) = 6, gamma >= 1.19 (the sphere/class-block route is better there, which is why p = 7
needed weights). p = 13, m = 4, r = 15: cubic cap 3816 − 1820 = 1996, Ax gives k <= 1989, bound
d(RM(35,4)) = 26, gamma 0.796 if a level set of <= 1989 points (mean 2197) exists and is
independent. p = 17, m = 4, r = 21: cap 12370 − 7295 = 5075 > 17^3 = 4913, bound d(RM(45,4)) = 68,
gamma = ln 16 / ln 68 = 0.657 for a generic cubic level set, n = 78608.

## 6. Case 1 sanity check (p = 7, m = 4, r = 7, weighted sphere support)

The theory above is for the unweighted grid; for the weighted sphere codes only the output-side
statements transfer, and they do: the cube class is the quartic level set structure Z =
{(g−1)(g−6) = 0}, the cap 326 − dim RM(3,4) = 291 is the class rank of WEIGHTED_SUPPORTS.md,
and the bound d(RM(16 + 4, 4)) = 5 is the mixed-word bound of the class-block lemma. The
certified codes there have d = 6, 7 > 5 only because the class block was scanned by MITM; the
theorem alone gives d >= 5 for any full-rank output subset of the class, and could not have given
7 (the block of 402 columns is where the accidental words live). On the unweighted grid at p = 7
the hypersurface route gives at best gamma 1.19 (cubic) / 1.22 (quadric: cap 204, d >= 7), so it
does not reopen p = 7; the jump to gamma 0.75 at p = 11 comes from d(RM(29,4)) = 22 versus
d(RM(19,4)) = 6.

## 7. Files

* `occupancy_p11.py` — line/plane/hyperplane occupancy of the k = 2290 draw (numbers in §3).
* `random_model_p11.py` / `.log` — rank of RM(13,4) on the k = 2290 set, random-code expectations, design table.
* `hyp_common.py` — helpers (random polynomial evaluation, rank, bound).
* `validate_small.py` / `.log` — brute-force validation of §4.1–4.2 at p = 7, m = 3 and p = 5, m = 4.
* `p11_hypersurface_scan.py` — random degree-delta polynomials at p = 11: |Z|, rank, independence.
* `p11_cubic_hunt.py` / `.log`, `p11_cubic_tested.jsonl`, `p11_cubic_best.json` — the cubic level-set hunt and the chosen code (coefficients, level, 1-based output columns).
* `verify_cubic_code.py` / `.log` — independent verification of the chosen code.
