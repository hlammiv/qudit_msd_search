# Level-set punctured RM codes: a proven γ < 1 code at p = 11

*2026-10-31. Found by the "algebraic route" agent (`scratchpad/routes/algebraic/NOTES.md`) while
looking for a certification method for n ~ 10⁴ codes; independently re-verified here
(`runs/rm_m4/p11_cubic_1309_bundle.json`, catalog `imported-p11-n13332-k1309-c7aae3253e15`).*

## The theorem

Let G = RM_p(r, m) (reduced degree ≤ r) with 3r < m(p−1), and let the output set O be the **full
zero set of a polynomial g of degree δ**, with the columns of G at O independent (|O| = k = rank of
G on O). Then the magic code [[pᵐ − k, k, d]]_p obtained by puncturing RM(r, m) at O satisfies

    d  ≥  d( RM_p( m(p−1) − r + δ − 1, m ) ),      and necessarily   k ≤ dim RM(r, m) − dim RM(r − δ, m).

*Proof.* For every h ∈ RM(r − δ, m) the function g·h has degree ≤ r and vanishes on O, so it is a
row of the shortened generator G0. A kernel word e (on the non-output points N) therefore satisfies
Σ_{x∈N} e_x g(x) h(x) = 0 for all h ∈ RM(r − δ, m); the vector e′ = g·e (zero on O) lies in
RM(r − δ, m)^⊥ = RM(m(p−1) − r + δ − 1, m). Since g ≠ 0 on N, supp e′ = supp e, so wt(e) ≥ the
minimum distance of that code, which is (p − b)·p^{m−a−1} for m(p−1) − r + δ − 1 = a(p−1) + b.
The cap on k is the dimension count: g·RM(r − δ, m) is a subspace of G0 of dimension
dim RM(r − δ, m), and dim G0 = dim G − k. ∎

This is the p = 7 class-block lemma with the block removed: there the outputs were a *subset* of a
quartic class set and the leftover class points formed a weakly checked block that had to be
scanned; taking O to be the *whole* zero set leaves no block, and the whole distance follows from
the dual RM distance with no scan. (Any O ⊆ {g = 0} gives the same bound for the part of a word
off {g = 0}; the leftover points of {g = 0} must then be certified separately, exactly as at p = 7.)

## The code: [[13332, 1309, ≥ 22]]₁₁, γ ∈ [0.648, 0.751]

p = 11, m = 4, r = 13 (= r_max), δ = 3. Dual code RM_11(29, 4): 29 = 2·10 + 9, distance
(11 − 9)·11 = **22**. Cap: dim RM(13,4) − dim RM(10,4) = 2320 − 1001 = 1319. Cubic zero sets in
F₁₁⁴ have size ≡ 0 (mod 11) (Ax); a random cubic's level sets have sizes 1232–1408 and about one in
eleven has exactly 1309 points; all twelve such sets found had independent columns. The chosen
one: `scratchpad/routes/algebraic/p11_cubic_best.json` (cubic coefficients, level, the 1309
output columns).

Verification done here, independently of the agent's code (`runs/rm_m4/p11_cubic_1309_bundle.json`
records it): the space of cubics vanishing on O is one-dimensional and its zero set is exactly O
(1309 points); the rank of RM(13,4) on O is 1309; g·h lies in RM(13,4) and vanishes on O for
random h ∈ RM(10,4); G0 has 1011 rows. Upper bound: a minimum-weight planar dual word (5 parallel
lines, 55 points) with 19 points in O gives an explicit kernel word of weight 36 (verified
G0 e = 0; the best of 4000 random planes, not exhaustive). So

| quantity | value |
|---|---|
| n, k | 13332, 1309 (rate 0.098) |
| d | 22 ≤ d ≤ 36 (proved / explicit word) |
| γ = ln(n/k)/ln d | **≤ 0.751** (0.648 if d = 36) |
| previous best at p = 11 | [[100,21,4]]₁₁ γ 1.126; the paper's analytic family 0.964 at n = 1.9·10⁷ |
| checks (dim A) | 1011: exact performance enumeration impossible, acceptance (1−ε)ⁿ |

This is the first γ < 1 code at p = 11 with a complete proof, at n = 1.3·10⁴. The output set is
not a random draw, so none of the sweep machinery is needed; the only computation is one rank.

## Context: why p = 11, m = 4 and not p = 7

At p = 7, m = 4 the budget allows δ ≤ 4 only through the class-block mechanism (δ = 4,
RM(3,4)^⊥ = RM(20,4), distance 5) and every full level set of degree ≤ 4 is larger than its cap, so
the theorem gives nothing beyond the class-block lemma. At p = 11, m = 4 the dual degree
26 + δ = 29 has residue 9 mod 10, which makes d(RM(29,4)) = 22 large, and a cubic zero set fits
under the cap. The random-draw window at p = 11 (`scratchpad/window_law.py`: k ≈ 2293, γ 0.865 if
d = 7) is better in γ but unprovable; the level-set code trades rate for a proof and a much larger
distance.

## The family (predicted, not built)

Same construction at r = r_max, δ = 3, O a cubic zero set of size ≈ p³ (always ≡ 0 mod p):

| p | r | dual code | d ≥ | cap on k | n ≈ | γ ≤ |
|---|---|---|---|---|---|---|
| 11 | 13 | RM(29,4) | 22 | 1319 | 13332 | 0.751 (built, proved) |
| 13 | 15 | RM(35,4) | 26 | 1996 | 26572 | 0.80 (needs a zero set ≤ 1996 points with full rank) |
| 17 | 21 | RM(45,4) | 68 | 5211 | 78608 | 0.657 |
| 19 | 23 | RM(51,4) | 76 | — | 123462 | 0.667 |
| 23 | 29 | RM(61,4) | 138 | — | 267674 | 0.627 |

Asymptotically n/k ≈ p and d ≈ p²/3, so γ → 1/2, the same limit as the paper's optimum, at block
size p⁴ (the paper's analytic family needs 10⁷–10¹¹). The p = 17 entry would beat our
[[237,52,6]]₁₇ (0.847) and the paper's [[77540,5981,15]]₁₇ (0.946) at the same block size. m = 5
and other δ are open.

## What the other routes said

- Integer programming / CP-SAT / SAT for "no word of weight ≤ w" (`scratchpad/routes/ilp/NOTES.md`):
  not competitive by 3–10 orders of magnitude at every size tried (cannot close weight 3 on the
  402-column p = 7 block in 300 s where MITM takes 0.45 s); the modular constraints give solvers no
  propagation. Closed.
- Information-set decoding (`scratchpad/routes/isd/NOTES.md`; the agent stopped before finishing its
  write-up, the data are complete): a finder and a sampled counter, never a certificate. Useful at
  n ≈ 400 (a weight-6 word in under a second, unbiased Monte-Carlo word counts); blind at n ≈ 12 000
  (success probability 4·10⁻⁹ per word per iteration). It found the 30x excess of weight-6 words in
  symmetric p = 7 blocks and confirmed the law on exact p = 11 restrictions.

## Next

1. Tighten the upper bound (exhaustive planar search; 3-flat words) and look for a better lower
   bound than 22 (the true distance is probably well above: the dual words of weight < 2·55 are all
   planar and meet O in ≤ 19 points).
2. Build the p = 13 and p = 17 members (rank computations at n = 2.9·10⁴ and 8.4·10⁴; the latter
   on lenore).
3. Cost: with d ≥ 22 and (1−ε)ⁿ acceptance, the rigorous curve of `qmsd/rigorous_curve.py` applies
   with the weight-22 term bounded by Singleton as well (no census needed) — compute the proved
   protocol saving at p = 11.
