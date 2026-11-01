# Route: information-set decoding (ISD) over F_p as a low-weight-word finder

Code: `scratchpad/routes/isd/isd_fp.py` (`find_low_weight(G0, p, w, iters, seed, ...)`, plus
`mitm_exhaustive`, `stern_iteration`, `subset_mitm_iteration`, the success-probability formulas
and `cost_table`).  Runs and result JSONs: `scratchpad/routes/isd/runs/`.  All compute here was
single-thread under `ulimit -v 4000000` on the local box (2026-10-31).

## 1. The algorithm, and why this regime is degenerate

The code whose low-weight words we want is C = ker(G0): G0 (r x n, r = rank = 15-60,
n = 400-12 000) *is* the parity-check matrix H, and C has redundancy r and rate > 0.99.

Classical ISD picks a random information set of C, i.e. a random set B of r independent columns
of H, eliminates H -> [Q | I_r] (free columns A first; Q = B^-1 H_A is r x (n-r)), and uses the
fact that for any choice of e_A the redundancy part e_B = -Q e_A is free: a word is caught iff
only a small weight t of it lies on A and the remaining w - t positions lie on B.
- Prange: t = 0 (word entirely inside B).  Lee-Brickell: enumerate weight-t patterns on A.
- Stern/Dumer: split A into Z1, Z2, enumerate weight-t1 / t2 patterns on each, collide on ell
  coordinates of Q c (the word must vanish on those ell B-coordinates), complete on the rest.

ISD beats enumeration because a random r-set catches a constant fraction of the weight **when
r is a constant fraction of n**.  Here r/n = 0.002-0.05, so a weight-w word has on average
w r/n = 0.01-0.3 of its support in B, and with probability (1 - r/n)^w ~ 0.75-0.99 the
information set absorbs nothing at all.  The per-word, per-iteration success probability of
p-ary Stern with windows z1, z2 and ell collision rows for a word of weight w' = t1 + t2 + j is

    P = [w'! / (t1! t2! j!)] * z1_(t1) z2_(t2) r_(j) / n_(w')  *  C(r - ell, j) / C(r, j)

(x_(a) = falling factorial; multivariate hypergeometric for the split of the support into
(Z1, Z2, B, rest), times the chance that the j B-positions avoid the ell collision rows).
Lee-Brickell(t): P = sum_{j >= w'-t} C(w', j) C(n-w', r-j) / C(n, r).  Random-subset MITM on s
columns: P = C(s, w') / C(n, w').  Every factor r/n or z/n that enters is the price of not
looking at all columns at once, and since the list size grows like z^{t} while the coverage
grows like z^{2t}, **the exhaustive MITM on all columns (z = n/2, j = 0) is always the most
cost-efficient variant and is the only one that proves anything**; the sub-exhaustive variants
only make sense when the exhaustive list does not fit in memory, and then their success
probability per word is a product of tiny factors (table in section 3).

What was implemented and verified (fuzzed against `qmsd.mindist_nb.weight_d_exists_fast` on
random small matrices for p = 3, 5, 7, 11; every reported word is re-verified exactly):
- `mitm_exhaustive`: projective-key MITM (both halves with (p-1)^{t-1} normalised patterns,
  keyed by the ell-part of the syndrome up to scaling; collisions completed exactly).
- `stern_iteration`: p-ary Stern/Dumer as above (t1, t2 in {1,2,3}; ell chosen so that the number
  of projective keys is ~4x the list size; windows z1 = z2 = z chosen from a list budget).
- `subset_mitm_iteration`: random s-subset + the repo's exhaustive kernels on it.
- Lee-Brickell t = 1, 2 (dominated by Stern; kept for the table).
"Not found after N iterations" is reported as a per-word miss probability (1-P)^N and, with the
law's E_w, as P(all words missed) ~ exp(-E_w (1 - (1-P)^N)).  **This is a model statement about
the randomness of the information sets, not a proof of absence.**

Empirical check of the formula: planted weight-6 word in a random 19 x 400 matrix over F_7,
Stern 2+2 (z = 190, ell = 8): hit rate 3/531 iterations = 0.0056 vs P = 0.0033 (consistent at
3 hits); Lee-Brickell t = 2: P = 5.2e-5 (not found in 200 iterations, as expected); subset-MITM
s = 200: P = 0.015 but 60 s per iteration with the repo kernel single-thread.

## 2. Validation on the known instances (Part 1)

Stern 2+2 throughout, lists ~1e5 (p = 7, z = 190) or 2.8e6 (p = 11, z = 745), one thread.

| instance | n x checks | law E_5 / E_6 / E_7 | target | result | P/iter/word | time |
|---|---|---|---|---|---|---|
| p=7 class block k=272 seed 3 (shift orbits, d_plane 8) | 400 x 19 | 0.009 / 3.7 / 1262 | weight 6 | found in 2, 4, 2, 5, 1 iterations (seeds 0-4) | 3.3e-3 (w6) | 0.1-0.6 s |
| same, census (240 s, 1552 iterations) | | | all weight-6 words | **110 distinct weight-6 words**, 586 hits (5.3 hits/word, census ~99% complete): N_6 ~ 115, the known word [16,67,122,200,210,328] among them | | 240 s |
| p=7 class block k=270 seed 19 (random, certified d = 7) | 402 x 21 | 2e-4 / 0.079 / 26.7 | weight <= 6 | **not found** in 3000 iterations; per-word miss probability 1.5e-6 (w6), 1e-66 (w5) | 4.5e-3 (w6) | 412 s |
| same | | | weight 7 | found in 62 iterations: [17,38,40,160,174,189,254] | 2.9e-4 (w7) | 6.5 s |
| same, census (3000 iterations) | | | all weight-7 words | 28 hits, 21 distinct; Monte-Carlo estimate N_7 = hits/(iters P) = **32 +- 6**, law 26.7, repo census ~30 classes; no word of weight <= 6 seen | | 319 s |
| p=11 G0 k=2290 seed 1 | 12351 x 30 | 1e-9 / 2.8e-5 / 0.50 | weight <= 6 | not found, 40 iterations | **3.6e-9** (w6) | 206 s (5 s/iter) |
| same | | | weight <= 7 | not found, 40 iterations | **1.4e-11** (w7) | 196 s |

Notes.
- The brief's "E_6 ~ 0.03-0.4, E_7 ~ 5-500" for the p = 11 candidate correspond to 26-27 checks;
  this candidate has rank(G0) = 30 = dim A (verified), so the law predicts E_6 = 2.8e-5 and
  E_7 = 0.5 for it, i.e. probably nothing of weight <= 7 to find at all.  Either way the
  finder is blind here: with P = 3.6e-9 per word and 5 s per iteration, one weight-6 word
  would need ~2.8e8 iterations (~45 years single-thread); the exhaustive 3+3 list is 3.1e13
  entries.  P(all missed) after 40 iterations is 1 - 4e-12 under the law.
- The weight-7 census of the k = 270 block is a genuine quantitative test of the law by sampling:
  the estimator hits/(iterations x P) is unbiased for the number of words with no enumeration;
  it agrees with the law (32 +- 6 vs 26.7) and with the exhaustive census on lenore (~30).
- **The k = 272 shift-orbit block has ~115 weight-6 words, 30x the law's 3.7.**  The excess is
  concentrated: columns 133, 16, 210, 254 each lie in ~20 of the 110 words (four hot columns,
  very likely one shift-orbit), the other 299 columns in <= 7 each.  So for this symmetric set
  the accidental-word count is NOT the random-code count; the sentence in WEIGHTED_SUPPORTS.md
  "symmetry ... cannot change the accidental count" is contradicted at the level of counts
  (existence, 0/18 survivors, is of course consistent with both).  See section 4 for the
  discriminating check (other shift blocks vs fresh random k = 272 sets).

## 3. The law at p = 11, m = 4 (Part 2)

### 3a. What the sampler can and cannot see at n = 12 351

Cost table for the p = 11 code (n = 12351, r = 30; list budget 3e6 entries, ~5 s/iteration):

| variant | list size | P/iter/word, w = 6 | P/iter/word, w = 7 |
|---|---|---|---|
| exhaustive MITM (proof) | 3.1e13 (w6) / 9.7e17 (w7) | 1 | 1 |
| Lee-Brickell t = 1 | 1.2e4 | 3.6e-13 | 8e-16 |
| Lee-Brickell t = 2 | 7.6e8 | 4.2e-10 | 1.2e-12 |
| Stern 2+2, z = 745, ell = 8 | 2.8e6 | **3.6e-9** | **1.4e-11** |
| Stern 3+3, z = 56 | 2.8e6 | 1.6e-13 | 2e-15 |
| subset-MITM, s = 53 | 2.3e6 | 4.7e-15 | 1e-21 |

With 24-26 checks (k = 2294-2296) the law gives E_6 = 0.4-50 words, but even 50 words x 3.6e-9
x 10 iterations = 2e-6 expected hits per candidate.  The specified protocol (ISD with a fixed
budget on 8 + 8 candidates) therefore cannot distinguish P(no weight-6 word) = 0.66 from 0.0:
its outcome is 0/8 and 0/8 whatever the truth.  Run anyway, for the record:

PART2_TABLE

The exact weight-4 scan `weight_d_exists_fast(G0, 11, 4)` needs the a = 2 left table of
C(12351, 2)(p-1) = 7.6e8 entries (~15 GB): EXACT4_RESULT

### 3b. An exact law test at p = 11, m = 4 that fits in memory

The law is a statement about the kernel-polynomial evaluation matrix; it can be tested exactly
where the 3+3 scan fits: shorten the k = 2290 candidate's G0 further (extra random outputs)
down to 13 or 14 checks and restrict to random column subsets of s = 100 or 120 points, where
E_6 is O(1) and the exhaustive projective MITM (list 1.6e7 / 2.8e7) runs in ~1-2 min per
subset under 4 GB.  Structured (plane) words of RM(13,4)^perp have weight >= 55 - occupancy
>= 12 (RM(6,2) over F_11 has minimum weight 55), so weight 6 is purely accidental here, as in
the m = 4 candidates themselves.  Counts per subset vs Poisson(E_6):

RESTRICTED_TABLE

## 4. Symmetric vs random k = 272 class blocks (weight-6 census by sampling, 150 s each)

SYMMETRY_TABLE

## 5. Verdict: what ISD can and cannot do for certification

- **Certification: nothing.**  In this regime (r << n) every ISD variant is either the exhaustive
  MITM (then it is the repo's existing certificate, at the same cost) or a sampler whose
  per-word success probability carries factors (z/n)^{t} (r/n)^{j}; "not found" converts into a
  model-based miss probability only, never into d > w.  The one thing the model gives is a
  *quantified* non-finding: for the p = 7 blocks (n ~ 400) 3000 Stern iterations (7 min, one
  thread) leave a per-word miss probability of 1e-6 at weight 6, which is a useful cheap
  pre-filter before the 6-min/32-thread exhaustive scan; for n = 12 000 at p = 11 the same
  budget leaves a miss probability of 1 - 1e-7, i.e. nothing.
- **Finding words / counting them: useful at n ~ 400.**  Stern 2+2 finds a weight-6 word of the
  k = 272 block in < 1 s (vs the 8 GB, 6-min/32-thread exhaustive 3+3 that does not fit under
  4 GB here at all), a weight-7 word of the k = 270 block in 6.5 s, and the Monte-Carlo
  estimator hits/(iterations x P) measures the number of weight-w words without enumeration
  (N_7 = 32 +- 6 vs law 26.7 and exhaustive ~30).  This is the practical use: a cheap census of
  the accidental-word population of a candidate, which is exactly how the x30 excess of the
  shift-symmetric block was seen.
- **At n = 12 351, p = 11: useless for weight >= 6**, and so is every other sub-exhaustive
  method I know of (Wagner's k-tree needs lists of p^{r/3} ~ 2e8 pairs and only finds a p^{-ell}
  fraction of the solutions): finding one weight-6 word among ~1e4 candidates' worth of
  columns needs the full 3e13-entry 3+3 scan.  Law validation at p = 11 therefore has to be done
  on instances where the exhaustive scan fits (section 3b), or with the weight-4/5 scans at a k
  where E_4, E_5 = O(1) (k >= 2303; the a = 2 table then needs ~15 GB, i.e. lenore, not here).
