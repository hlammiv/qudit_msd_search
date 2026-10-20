# Hail-mary hunt for power-of-2 ring-qudit MSD codes (n < 10⁴, γ < 2, universality-enabling)

*2026-10-19/20.  Workflow `wf_c8b0c8b7-d1c`: 12 explorers + 6 follow-ups completed; the verification and
synthesis stages failed (out of credits), so this document is the hand-written synthesis.  Per-agent
artifacts: `runs/diverse/hailmary/<tag>/`, scripts: `scratchpad/hailmary/<tag>/`.*

## 1. Verdict

**No family produced a universality-enabling protocol with γ < 2.**  The decisive criterion, derived during
the hunt (accounting-audit follow-up) and implemented as `rgt.escape_rank`:

> A hosted transversal layer with k logical outputs produces only **r** universality-enabling (T-like) states
> per round, r = rank_F₂ of the escape form B[i,l] = [Δ₂ₑᵢΔ₂ₑₗF ≠ 0] (for 0/1 lifts with a uniform carrier,
> B = O·diag(c)·Oᵀ, O = logical reps mod 2).  The overhead exponent is **γ_T = log(n/r)/log d**, not log(n/k)/log d.
> By the finite-cheap-group argument (UNIVERSALITY_CHECK §7) the other k − r outputs cannot substitute.

Every sub-2 candidate has r ≪ k:

| candidate | claimed γ (log n/k) | r | γ_T | status |
|---|---|---|---|---|
| simplex5 [[31,6,3]]₄ (non-uniform layer) | 1.49 | 1 | 3.13 | rank-1 (checked directly) |
| simplex-monomial [[127,17,3]] … [[1023,182,3]]₄ | 1.83 … 1.57 | 1 | ≥ 4.4 | rank-1 (only all-ones logical escapes) |
| flat family [[31+16F, 1+5F, 2]]₄ | → 1.70 | 1 | ≥ 5.5 | rank-1 |
| HGP(simplex+A, quadruply-even columns) [[8322,1309,3]]₄ etc. | 1.68–1.99 | 51–113 | 3.3–4.0 | B factorizes, first factor rank 1 |
| BCH-dual subcodes [[511,19,7]]₄, [[8191,287,15]]₄ | 1.69, 1.24 | 1 | 3.2–3.9 | magic-rank theorem |
| Gold [[511,28,5]]₄, [[511,10,5]]₄ | 1.80, 2.44 | 1 | 3.87 | rank-1 |
| **pRM(2,9) [[492,20,4]]₄** (pre-hunt) | 2.31 | **20** | **2.31** | full rank — still the best |
| [[15,1,3]]₄ non-uniform / [[31,1,3]]₄ | 2.46 / 3.13 | 1 | 2.46 / 3.13 | full rank (k=1) |

## 2. Theorems and mechanisms established by the hunt

1. **Escape form** (accounting audit, theory-bound): B = Δ₂Δ₂F/8 is the complete A₂-escape invariant; logical
   Cliffords act by B → AᵀBA; minimal number of single-qudit escapers = symmetric rank of B.
2. **Escape set** (theory-bound, Thm 1): the escape set E of any hosted layer satisfies E ∈ Σ^⊥, Σ = R_X∗R_supp∗R_supp
   (binary Schur cube of the residue codes); tight in every tested case.  Physical escape is impossible when Σ = F₂ⁿ
   (explains random CSS, Hamming, RM(2,6) punctured, GR evaluation codes).
3. **Magic-rank theorem** (rm-geometry): if all logical reps lie in a quadruply-even code punctured at P, then
   r ≤ |P|.  Enlarging L inside the RM space adds non-escaping logicals only.
4. **Puncture geometry**: ring-punctured RM(r,m) (4r < m) has d_Z = 2^{r+1} − (max (r+1)-flat occupancy of P) and
   r = |P| (full).  For r = 2: d ≤ 4 and |P| ≤ Sidon size ⇒ γ_T ≥ m/4 ≥ 2.25 — **RM(2,m) can never beat 2.**
5. **Non-CSS invariance** (noncss-solver): the hosted module and induced logical depend only on the CSS shadow
   (proj_X S, Z(S)); non-CSS twists can only raise X-type distance.
6. **Z₈/Z₁₆ reduction** (z8-z16): all level-≤3 single-qudit tables are q = a·y + b·C(y,2) + c·C(y,3) with the escape
   bit at the mod-8 level of c; free Z_{2d} lifts reduce to Z_d ⇒ **larger 2-power rings do not help** for lifted codes.
7. **Z₄-native codes** (Kerdock/Preparata/Goethals/Type-II QR/Klemm, ~450 CSS codes): hosted modules are
   Clifford-only or non-escaping whenever k ≥ 2; [[16,6,4]]₄, [[32,20,4]]₄ Kerdock CSS codes host nothing.
8. **Permutation / level-4 / semi-Clifford gates** (weird-gates): non-affine transversal permutations are hosted by no
   k ≥ 1 CSS code; level-4 two-step teleportation and semi-Clifford frames do not move γ.
9. **Across-block 2-body** layers on 0/1 lifts reduce to Q2/Q3 without Q1 plus a quadruple condition; high-distance
   Clifford-only codes stay Clifford-only.
10. **Concatenation** cannot beat the inner code (mediant theorem); HGP preserves hosting multiplicatively but B factorizes.

## 3. The one surviving route: ring-punctured RM(3,13) with Sidon puncture sets

From (4): RM(3,13) (n = 8192, quadruply even since 4·3 < 13, dual distance 16), P ⊂ F₂¹³, k = r = |P|,
d_Z = 16 − t, t = max 4-flat occupancy.  If P is a **Sidon set** (no a+b = c+d), its affine-dependency code is even
with distance ≥ 6, so s points in a 4-flat need an [s, s−5, ≥6] binary code — impossible for s ≥ 7 (Griesmer).
Hence t ≤ 6 and **d_Z ≥ 10 from flat-type words**, with γ_T = log((8192−k)/k)/log 10 < 2 ⇔ k ≥ 82.
The best known Sidon set in F₂¹³ has size **82** (Czerwinski–Pott, arXiv:2212.05887; known range 82–125):
⇒ [[8110, 82, 10]]₄, γ_T ≈ 1.995 (1.92 if t ≤ 5 gives d = 11).  Relaxed sets with t ≤ 7 (d=9) need k ≥ 100, t ≤ 8
(d=8) need k ≥ 127; k ≈ 250 at t = 8 would give γ_T ≈ 1.66.  In code language: t ≤ T ⇔ the dependency code of
(1,p)_{p∈P} has generalized Hamming weight d_{T−4} ≥ T+2.

Open certification issue: the flat argument bounds only words of weight < 32 of RM(9,13) (Kasami–Tokura forms);
a full d_Z certificate at n ≈ 8100 needs either a structural argument for heavier words or ISD evidence.

### 3a. Realized (2026-10-20): two constructions just below γ_T = 2

| construction | P ⊂ F₂¹³ | code | t (max 4-flat occ.) | γ_T |
|---|---|---|---|---|
| **Sidon-82** | parity-check columns of the even [82,68,6] code (codetables.de, BKLC n=82 k=68) — verified Sidon (3321 distinct pair sums) | [[8110, 82, 10]]₄ | 6 (witness: 6 points summing to 0; ≤6 by Griesmer) | **1.995** |
| **split-128** | S₁×{0} ∪ (A·S₁+b)×{1}, S₁ = {(x,x³) : x ∈ F₆₄} ⊂ F₂¹² (Sidon), A random invertible | [[8064, 128, 8]]₄ | ≤ 8 (a 4-flat meets each half in a 3-flat; Sidon ⇒ ≤ 4 per 3-flat) | **1.992** |

Verified for both: rank RM(3,13)|_P = |P| (so k = |P|, free lift, logical reps with l_p|_P = e_p);
escape form B = I ⇒ r = k (every logical output is universality-enabling); hosting of the uniform carrier is
automatic (RM(3,13) is quadruply even, 4·3 < 13; Q1–Q3 via |s∧u| = (|s|+|u|−|s+u|)/2 etc.; s vanishes on P).
Data: `runs/diverse/sidon82_F2_13.json`, `runs/diverse/split128_F2_13.json`.  With s(12) ≥ 66 the split gives
[[8060,132,8]]₄, γ_T = 1.977.

**Distance status.** d_Z ≤ 10 (resp. ≤ 8) by flat witnesses; the matching lower bound is proved for all
RM(9,13) words that are 4-flats (via t) — words of weight 24/28/30 (Kasami–Tokura forms, supported in ≤ 6-flats,
where a Sidon subset has ≤ s(6) = 9 points) also satisfy |c \ P| ≥ 15.  **Not yet certified:** heavier words
(32 ≤ |c| ≤ |P| + 9).  Equivalent statement: P ∪ R imposes independent conditions on cubics for every R of size
≤ d−1 off P.  If d dropped by one, γ_T would rise to ≈ 2.09 / 2.13.

**Why this paradigm saturates at 2.** For ring-punctured RM(r,m) one needs 4r < m (n = 2^m), d ≤ 2^{r+1} and
|P| ≲ 2^{m/2+1} (Sidon-type limits), so γ_T ≳ (m/2)/(r+1) ≥ 2r/(r+1) + O(1/r): at most slightly below 2 for every r,
and r = 3 (n = 8192) is the only option with n < 10⁴.  γ_T ≪ 2 needs a mechanism other than puncturing a
quadruply-even code.

## 4. Next steps (ranked)
1. Obtain/construct an explicit Sidon set |P| ≥ 82 in F₂¹³; compute t exactly (C(82,5) ≈ 2.6·10⁷ 5-subsets);
   verify rank(RM(3,13)|_P) = |P|, escape rank r = |P|, hosting (binary Q1–Q3, automatic), and run an ISD search for
   Z-logicals of weight < 16 − t.
2. Search for larger P with t ≤ 7 or 8 (GHW-constrained sets) — the only lever toward γ_T ≈ 1.7.
3. If 1–2 succeed: prove the d_Z lower bound for words of weight ≥ 32 (moment/occupancy argument).
