# Faster Minimum-Distance Algorithms for Punctured Reed-Muller / Triorthogonal Codes

**Target problem.** For the qudit magic-state-distillation code search (arXiv:2510.10852) over
F_p (focus p=3), the single blocking computation is the **exact minimum distance** — and ideally
A_d, the number of minimum-weight codewords — of a **punctured Reed-Muller dual code**.

Concrete blocked case (p=3, m=7): the triorthogonal space is RM_3(4,7), length N0 = 3^7 = 2187,
dimension 274. Puncturing/shortening at k columns gives a small generator
G0 = `[2187-k, 274-k]` whose **dual** G0^perp = PRM_3(9,7;S) = `[2187-k, 1913]` carries the quantum
distance we need. Unpunctured dual distance d_RM(9,7) = 18; distances of interest are d ~ 6..18, and
gamma<1 needs d > n/k (e.g. k ~ 150..274 demanding d ~ 6..13).

Equivalently: given parity-check G0 (r = 274-k rows, N = 2187-k cols over F_3), find the **fewest
F_3-linearly-dependent columns**.

Current method (`mindist.py`): meet-in-the-middle (MITM) on columns of G0, cost ~ C(N, d/2)·2^(d/2).
At N~2000, d~10 that is ~10^15 **per candidate** — infeasible across a search.

---

## 1. Ranked comparison

All complexity figures are **per candidate code**, F_3, N~2000, independently recomputed (see
`scratchpad`/inline). "Exact?" distinguishes *certified* (a proof d = value) from *upper-bound*
(a witness of weight = value, minimality unproven).

| Rank | Approach | Exact / certified? | Cost for our params (per candidate) | Feasible at m=7? | Verdict |
|---|---|---|---|---|---|
| **1** | **MacWilliams from small dual** (enumerate all 3^(274-k) words of G0, Krawtchouk-transform to primal weight distribution) | **EXACT + certified; yields full A_d** | 3^(274-k): k=240→1.7e16, k=252→3.1e10, k=254→3.5e9, k=260→4.8e6, k=270→81 | **Yes, but only k ≳ 252** (dual dim ≲ 22) | **winner (in its corner)** |
| 2 | **Structured RM enumeration** (DGM / Kasami-Tokura-Azumi / Leducq affine-flat low-weight codewords of RM_3(9,7), minimize support outside S) | **Mixed**: exact+certified only for tiny k (puncture < ~4th-weight value); else high-confidence **upper bound + candidate A_d** | weight-18 layer = 2.9e9 structured supports → ~5e10 cheap ops/candidate (list fixed across candidates; regenerate on the fly) | Yes as a **fast screen / tight upper bound + A_d for ALL k**; not a certifier for k ≳ 20–30 | promising |
| 3 | **q-ary ISD (Stern/Peters) in syndrome form on G0** + QDistRnd-style sampling | **Heuristic upper bound only** (probabilistic, empirical convergence) | optimized q-ary Stern, k=150: w=6~1e8.5, w=10~1e12.5, w=13~1e15.9, w=18~1e21.6 (2–5 orders below MITM in the gamma<1 band) | Yes as an **upper-bound / candidate-rejection oracle** for d ≲ 13; **not** a certifier | promising (search only) |
| 4 | **Brouwer-Zimmermann / off-the-shelf tools** (Magma, GAP-Guava, Sage; Algorithm 994; quantum-BZ 2408.10743) | EXACT + certified **in principle** | On the high-rate dual (dim 1913): floor(n/K)=1 ⇒ enumerate to weight d−1: C(1913,d−1)·2^(d−1) = 6.8e15 (d=6) … 2.1e46 (d=18) — **worse than MITM** | **No** for the target quantity | does-not-help |

### Why the obvious tool (BZ / Magma) fails here
BZ's only lever is the number of disjoint information sets, t = floor(n/K). The quantum distance
lives on G0^perp with K = 1913 and n = 2187−k ≤ 2187, so 2187 − 1913 = 274 < 1913 ⇒ **t = 1 for
every k**. BZ then degenerates to enumerating all codewords of weight ≤ d−1, i.e. C(1913, d−1)·2^(d−1),
which is ~5 orders of magnitude **worse** than the MITM already in use. BZ cannot transfer a distance
across duality, so running it on the cheap small side G0 (dim 274−k) returns the *wrong* code's weight.
The fastest open BZ codes (Algorithm 994; the quantum/symplectic-BZ of arXiv:2408.10743) are
additionally **GF(2)-only** and would need reimplementation for F_3. GAP-Guava `MinimumWeight` is
deterministic BZ over binary **and ternary** (field match for p=3) but is gated by the same t=1 wall.

---

## 2. Winning recommendation and WHY

**Primary (exact, certified): MacWilliams transform from the enumerated small dual G0.**

- G0 is `[2187-k, 274-k]`. When the dual dimension r = 274−k is small, enumerate **all** 3^r codewords
  of G0, accumulate its weight enumerator W_G0(x,y), then apply the **q-ary MacWilliams identity**
  (Krawtchouk transform, q=3) to obtain the **complete weight distribution of G0^perp = PRM_3(9,7;S)**.
  This yields the exact minimum distance **and every A_d** in one shot — exactly the quantity the search
  needs — with a *certificate* (it is a closed-form linear transform of an exhaustively-computed
  enumerator, not a sample).
- Cost is 3^(274−k), independent of d. Practical workstation reach ~3^22 ≈ 3·10^10. This is the **only**
  method in the survey that is simultaneously (i) exact+certified, (ii) feasible at N~2000 over F_3, and
  (iii) returns A_d for free. It was prototyped end-to-end over F_3 (validated against brute force) in
  two independent verifications.

**Why it beats the alternatives in its regime:** MITM and ISD only ever give an *upper bound* on d
(and not A_d); BZ on the high-rate dual is worse than MITM; structured RM enumeration is uncertified
once the puncture set is larger than the classified weight ceiling. MacWilliams sidesteps all of this
because the *dual is small* — the one piece of exploitable structure that gives an exact answer cheaply.

**Secondary (fast screen + A_d, all k): structured RM enumeration.** Low-weight codewords of RM_3(9,7)
are not arbitrary — DGM/Leducq prove minimum-weight words are evaluations of products of affine forms
= unions of parallel affine flats. d_min = 18 = 2·3^2 (union of 2 of 3 parallel 2-flats inside a
3-flat of AG(7,3)), with A_dmin = (q−1)q^t [m,t]_q [m−t,1]_q C(q,s) = 5,849,021,178 codewords
= 2,924,510,589 distinct supports = 74,987,451 three-flats × 13 directions × 3 coset-pairs (count
independently re-derived; DGM formula checked on RM_3(1,2)→24 weight-6 words via galois). For each
puncture set S, compute |supp \ S| over this fixed structured list to get a **tight upper bound on the
punctured distance and candidate A_d at ~5·10^10 ops/candidate** — ~10^3–10^4× cheaper than MITM and
returning A_d, usable for **every k** to drive and prune the search.

**Tertiary (search acceleration only): q-ary ISD / QDistRnd.** Work in parity-check form on G0 (small
syndrome space r = 274−k makes the Stern collision step prune hard). Use it purely as a fast
upper-bound oracle to **reject** candidates (decide gamma ≥ 1) — punctured-RM duals have large A_d, so
finding *one* low-weight word is much cheaper than the worst case. Do **not** invest in BJMM/MMT:
Canto-Torres–Sendrier prove all ISD variants share the same exponent at sublinear weight (our
w/n ≈ 0.005), so plain q-ary Stern is optimal.

---

## 3. Regime unlocked

| (k, d) regime | What becomes available | By which method |
|---|---|---|
| **k ≳ 252** (dual dim 274−k ≲ 22), any d | **EXACT d AND full A_d, certified**, ~10^10 ops/candidate | MacWilliams-from-dual (winner) |
| k ≳ 240 (dual dim ≲ 34) | exact+certified but ~10^16 — borderline workstation, feasible on a cluster | MacWilliams-from-dual |
| **all k (incl. k~150)**, d ≲ ~18 | **high-confidence upper bound on d + candidate A_d**, ~5·10^10 ops/candidate; exact only for very small k (puncture below the ~4th-weight value) | structured RM enumeration |
| k ~ 150..240, d ≲ 13 | fast **candidate rejection / upper bound**, 10^8–10^16 (2–5 orders below MITM) | q-ary ISD / QDistRnd |
| k ~ 150 (dual dim 124), needing certified d ~ 6..13 | **certified by NO method** — the genuine open gap | — |

The **open-qutrit / high-puncture corner (large k, low d)** — where gamma<1 is reached with heavy
shortening — is exactly the corner MacWilliams makes **exactly and certifiably** reachable, with A_d.
The **hard middle (k ~ 150, dual dim ~124)** remains uncertified by any surveyed method; there the best
available is the structured-enumeration upper bound (very likely tight, since high-weight RM words
rarely concentrate on a generic moderate S) cross-checked against ISD.

---

## 4. Concrete implementation plan

All steps reference only verification-confirmed sources (Section 5).

1. **Build G0 correctly (shared).** From RM_3(4,7) construct the shortened generator G0 = `[2187-k, 274-k]`
   for the chosen puncture/shorten set S, in `galois` over GF(3). **Verify G0 is exactly the shortened
   generator** (so its dual is the intended PRM) before trusting any A_d.

2. **EXACT path — MacWilliams-from-dual (primary, for k ≳ 245–252).**
   - Enumerate all 3^(274−k) codewords of G0 by iterating message vectors; accumulate the weight
     enumerator W_G0 over a length-(2187−k) histogram. Use bit/trit-sliced packing and parallelize over
     message blocks; reach is ~3^22 on a workstation, ~3^34 on a cluster.
   - Apply the q-ary **MacWilliams identity** (Krawtchouk polynomials K_j, q=3) to W_G0 to get the
     full weight distribution {B_w} of G0^perp; the smallest w>0 with B_w>0 is the certified d, and
     B_d = A_d. Validate the transform on a small `[n, dim≤8]` F_3 code against brute force (done in
     verification; reproduce as a unit test). Sanity check: sum B_w = 3^1913 / 3^(274−k) and all B_w
     integers. References: MacWilliams identity standard; see Delsarte-Goethals-MacWilliams [DGM] and any
     coding text for the q-ary Krawtchouk form.

3. **SCREEN path — structured RM enumeration (secondary, all k, returns candidate A_d).**
   - Precompute the **affine-flat generator** for RM_3(9,7) min-weight supports: nested loop over
     3-flats of AG(7,3) (74,987,451 of them) × 13 plane-directions × 3 coset-pairs → 2.92e9 distinct
     weight-18 supports (DGM/Leducq form: union of 2 of 3 parallel 2-flats). Regenerate on the fly
     (storing all is 100s of GB); prune with the geometric condition |T ∩ S| ≤ 18 − d.
   - **Exclude supports T ⊆ S** (they puncture to zero — a real bug surfaced in verification).
   - For each candidate S, the min over structured T (not ⊆ S) of |T \ S| is the upper bound; count
     of minimizers is candidate A_d. Extend with the 2nd/3rd/4th-weight Leducq classes (up to ~2·d_min)
     as a constant-factor refinement. **Honest scope:** only the first ~4 distinct weight *values* are
     classified, so this certifies only when the puncture is below that ceiling; otherwise it is an
     upper bound. References: [DGM]; Leducq arXiv:1001.2554, arXiv:1203.5244, arXiv:1203.4592; Rolland
     arXiv:0902.0058; Kasami-Tokura-Azumi 1976; 4th-weight Inf.Comput. 2023 (DOI 10.1016/j.ic.2023.105110);
     count formula Dang-Ghorpade arXiv:2504.21816.

4. **REJECT path — q-ary ISD oracle (tertiary, search throughput).**
   - Implement Stern/Peters q-ary ISD in **parity-check (syndrome) form on G0** (small r = 274−k
     buckets prune the collision step). Use as a low-cost upper-bound to reject candidates with
     gamma ≥ 1. Optionally call **QDistRnd** (GAP, GF(q), CSS-aware) directly on G0 for the same role.
     Do **not** implement BJMM/MMT (Canto-Torres–Sendrier: no asymptotic gain at sublinear weight).
   - References: Peters ePrint 2009/589; Interlando et al. arXiv:1812.10955 (q-ary cost formulas);
     Canto-Torres–Sendrier PQCrypto 2016; QDistRnd arXiv:2308.15140 / JOSS 7(71):4120.

5. **Orchestration.** Route each candidate by k: k ≳ 245 → MacWilliams (exact+A_d); else structured
   enumeration for screen+candidate-A_d, with ISD/QDistRnd cross-check of the upper bound. Cache the
   fixed structured support list across the whole search. For the hard middle, report the upper bound
   explicitly as **uncertified**.

---

## 5. Verified references (only these were confirmed; author/issue corrections applied)

**Exact-method core**
- M. Grassl, "Searching for linear codes with large minimum distance," *Discovering Mathematics with
  Magma*, Springer 2006, pp. 287–313. DOI:10.1007/978-3-540-37634-7_13. *(Canonical BZ description.)*
- F. Hernando, M. F. Igual, G. Quintana-Ortí, "Algorithm 994: Fast Implementations of the
  Brouwer-Zimmermann Algorithm…," *ACM TOMS* 45(2), 2019. DOI:10.1145/3302389. *(GF(2)-only.)*
- S. Bouyuklieva, I. Bouyukliev, "An Extension of the Brouwer-Zimmermann Algorithm…," *Mathematics*
  9(19):2354, 2021. DOI:10.3390/math9192354. *(Author correction: NOT "Sanvicente"; F_q / info-set
  ordering.)*
- F. Hernando, G. Quintana-Ortí, M. Grassl, "Fast Algorithms and Implementations for Computing the
  Minimum Distance of Quantum Codes," arXiv:2408.10743; *ACM Trans. Quantum Computing*, DOI:10.1145/3795877.
  *(Author correction: NOT "Berent"; speedup over Magma stated as more than one order of magnitude;
  GF(2)/symplectic only.)*
- P. Lisoněk, L. Trummer, "Algorithms for the minimum weight of linear codes," *Adv. Math. Commun.*
  10(1):195–207, 2016. DOI:10.3934/amc.2016.10.195. *(BZ lower bound; info-set count drives cost.)*

**ISD / heuristic upper bound**
- R. Canto Torres, N. Sendrier, "Analysis of Information Set Decoding for a Sub-linear Error Weight,"
  *PQCrypto 2016*, LNCS 9606, pp. 144–161. DOI:10.1007/978-3-319-29360-8_10; HAL hal-01244886.
  *(All ISD variants share the exponent at sublinear weight — don't use BJMM/MMT.)*
- C. Peters, "Information-Set Decoding for Linear Codes over F_q," *PQCrypto 2010*; IACR ePrint
  2009/589; DOI:10.1007/978-3-642-12929-2_7. *(q-ary Stern cost formulas.)*
- J. C. Interlando, K. Khathuria, N. Rohrer, J. Rosenthal, V. Weger, "Generalization of the
  Ball-Collision Algorithm," arXiv:1812.10955. *(q-ary ISD success-prob / per-iter cost.)*
- L. P. Pryadko, V. A. Shabashov, V. K. Kozin, "QDistRnd: A GAP package for computing the distance of
  quantum error-correcting codes," *JOSS* 7(71):4120, 2022. DOI:10.21105/joss.04120; arXiv:2308.15140.
  *(Probabilistic GF(q) CSS upper bound; issue corrected to 7(71).)*

**RM / GRM low-weight structure**
- P. Delsarte, J.-M. Goethals, F. J. MacWilliams, "On generalized Reed-Muller codes and their
  relatives," *Information and Control* 16 (1970) 403–442. DOI:10.1016/S0019-9958(70)90214-7.
- T. Kasami, N. Tokura, S. Azumi, "On the weight enumeration of weights less than 2.5d of Reed-Muller
  codes," *Information and Control* 30 (1976) 380–395. DOI:10.1016/S0019-9958(76)90355-7.
- E. Leducq, "A new proof of Delsarte, Goethals and Mac Williams theorem on minimal weight codewords of
  GRM codes," arXiv:1001.2554.
- R. Rolland, "The second weight of generalized Reed-Muller codes in most cases," arXiv:0902.0058.
  *(Author correction: Rolland, not Leducq.)*
- E. Leducq, "Second weight codewords of generalized Reed-Muller codes," arXiv:1203.5244;
  "Remarks on low weight codewords of generalized affine and projective RM codes," arXiv:1203.4592.
- "On the codewords of generalized Reed-Muller codes reaching the fourth weight," *Information and
  Computation* 293 (2023). DOI:10.1016/j.ic.2023.105110.
- Dang, Ghorpade, "Enumeration of minimum weight codewords of affine Cartesian codes," arXiv:2504.21816.
  *(Closed-form A_dmin count, used and verified here.)*

**Tooling**
- GAP-Guava manual ch.4 (`MinimumWeight` deterministic BZ, **binary and ternary only**;
  `MinimumDistanceLeon` probabilistic, binary only); Magma `MinimumWeight` (BZ, general GF(q),
  closed-source); SageMath `LinearCode.minimum_distance` (wraps Guava/Magma); M. Grassl, codetables.de
  (bounds only, no n~2000 entries). *(Standard, confirmed; same high-rate wall.)*

---

## 6. Honest bottom line

**Partly reachable — the open-qutrit / high-puncture corner opens exactly; the hard middle does not.**

- **m=7, k ≳ 252 (open-qutrit-style heavy shortening):** **SOLVED, exactly and with A_d.**
  MacWilliams-from-the-small-dual costs 3^(274−k) ≲ 3^22 ≈ 3·10^10 — a workstation computation that
  returns the certified minimum distance *and* the full weight distribution. k ≳ 240 (≈10^16) is
  reachable on a cluster. This is a genuine, certified unblock of the regime the search most cares
  about (large k, low d, gamma<1).

- **m=7, k ~ 150 (dual dim ~124, certified d ~ 6..13):** **STILL OUT OF REACH for a certificate.**
  No surveyed method certifies it: BZ degenerates (t=1, worse than MITM), MacWilliams needs 3^124,
  ISD/QDistRnd give only probabilistic upper bounds, and structured RM enumeration is uncertified
  above the ~4th-weight ceiling (puncture ≫ ~30). What *is* newly feasible there is a **fast,
  high-confidence upper bound + candidate A_d** (structured enumeration at ~5·10^10 ops, cross-checked
  by ISD) — enough to **drive and prune the search** and to break the d ≤ 6 cap of the current MITM,
  but **not a proof**.

- **No overclaiming:** the project's stated bottleneck — *exact, certified* distance across the whole
  gamma<1 window — is **not** fully removed. It is removed in the high-k corner and replaced by a
  strong heuristic (with A_d) elsewhere. The honest path forward for the hard middle is either to
  reformulate the search so the target distance sits on a low-dimension code (where BZ/MacWilliams
  become exact), or to push the GRM weight classification (DGM/Leducq) past the 4th weight toward
  ~|S|+d so structured enumeration can certify larger punctures — both are open research, not
  off-the-shelf.
