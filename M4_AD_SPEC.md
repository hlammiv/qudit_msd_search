# Spec: m=4 structured A_d enumerator → optimize the flagship [[519,106,5]]_5

**Status:** proposed 2026-06-29. Goal: compute and minimize A_5 of the paper's γ<1 ququint code
`[[519,106,5]]_5` (p=5, m=4) — the one p=5 code that crosses γ<1 (γ=0.987) — by extending the structured
A_d machinery to m=4, since MacWilliams is infeasible there (dim G0 = 16 ⇒ 5¹⁶ ≈ 1.5e11).

## 0. The single blocker (scoped 2026-06-29)
`qmsd/structured_ad.py` is already **general-(p,m)** (it enumerates codewords by minimal affine-span dim j over
`_flats(m,j,p)`, with a `jmax` cap, and returns `{distance, A_d, low_weight_histogram, dim_histogram, jmax,
exact}`). It works at m=4 for `jmax ≤ 2` but **overflows at j ≥ 3**: `_mitm` (line ~227) calls `_powers(p, r)`
which raises `OverflowError` when `p**rows > int64` (e.g. `5**29`). `rows` = redundancy of the restricted code
`C_F = RM_p(rtilde-(m-j)(p-1), j)` on a j-flat, which is large for j=3 flats at m=4.

**Fix:** make the MITM syndrome encoder overflow-proof — encode each syndrome vector via a **hash (tuple→dict)
or Python-bigint p-adic** instead of the int64 p-adic in `_powers`/`_encode`/`_half_sums`/`_mitm`. The MITM only
needs syndrome **equality** (left_sum + right_sum = 0), so a hash-map match is sufficient (verify on hash hit).

## 1. Constraints
- **Do NOT regress the validated p=3 path.** `structured_ad` reproduces 6/7 paper qutrit A_d exactly. Preferred:
  refactor `_mitm` to a hash encoder (strictly more general) AND re-run the existing structured_ad/qutrit tests
  to prove no regression. Acceptable alternative: a new `qmsd/structured_ad_m4.py` reusing the flat enumeration
  with a hash MITM. Either way ADDITIVE to behavior — keep p=3 results identical.
- **No brute MITM, no `C(p^4,a)` arrays, no `p^dim_G0` full enumeration at the flagship** (5¹⁶ infeasible).
- Bounded discipline: write incrementally, save often, NEVER run a single call > ~2 min unbounded; if slow, cap.

## 2. Key facts (p=5, m=4)
- `r_max(4,5)=5`, `rtilde = 4·4 − 5 − 1 = 10 = 2·4 + 2` (a=2,b=2) ⇒ `d_RM = (5−2)·5^(4−2−1) = 3·5 = 15`.
- The minimum-weight RM₅(10,4) codewords are **plane-supported (j=2, 25-point 2-flats), weight 15**; punctured
  to d=5 needs 10 punctures on the support. Higher-j (j=3,4) codewords may also puncture to 5 (the "flat cap").
- Flagship: n=519, k=106, d=5, dim G0 = 16.

## 3. Tasks (stop at any failure; report)
1. **Fix the encoder** (§0) so `structured_ad(5,4,5,punc, jmax=None)` runs at m=4 without overflow.
2. **Validate (decisive):**
   - **(a)** vs MacWilliams at p=5 m=4 **high-k** (dim G0 ≤ 9 ⇒ 5⁹≈2e6 feasible via `qmsd.weightdist.
     exact_distance_and_Ad`): structured `A_d` (full, jmax=None) MUST equal MacWilliams `B_d` EXACTLY. These
     high-k codes are full-span/high-j dominated, so they genuinely exercise the j≥3 fix.
   - **(b)** consistency: structured_ad restricted to a single 3-flat reproduces the m=3 enumerator
     (`qmsd/structured_m3.py`) on that 3-flat (optional but strong).
   - **(c)** p=3 regression: existing qutrit A_d (572/1104/1128 etc.) still reproduced.
3. **Reach the flagship d=5:** generate puncture sets giving d=5 at p=5 m=4 (extend the flat-spread idea to
   AG(4,5), or search). Determine the **jmax that is provably exact for d=5** (min-weight is j=2; check whether
   j=3,4 contribute via the validated full run on a feasible proxy, or argue from `_min_span_weight`).
4. **Compute + optimize A_5** of `[[519,106,5]]_5`: baseline A_5, then minimize over puncture sets (flat-spread /
   cap-set / local search), as we did for `[[112,13,3]]_5` (478→396) and qutrit `[[206,37,4]]_3` (880→572).

## 4. Deliverable
`qmsd/structured_ad_m4.py` (or a refactored `_mitm` + note) + tests, and a summary: did the encoder fix validate
(structured == MacWilliams at p=5 m=4)? the flagship `[[519,106,5]]_5` baseline A_5 and the optimized A_5 (+
puncture set). Honest scope: if jmax-exactness for d=5 can't be certified, say so and report A_5 as a bound.
