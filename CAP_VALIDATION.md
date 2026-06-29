# Validation of the m=7 cap qutrit code `[[1968,219]]_3`

**Question.** Is the minimum distance `d` of the cap qutrit code `>= 10`?

- `d >= 10`  =>  `gamma = log(1968/219)/log(d) <= log(1968/219)/log(10) = 0.9536 < 1`  (the gamma<1 claim holds)
- `d <  10`  =>  `gamma >= 1`  (the gamma<1 claim is **refuted**)

**Verdict: REFUTED. `d = 1`, established rigorously and exactly. `gamma >= 1`.**

---

## 1. The object under test

Puncture set `S` = the 219 cap points (`cap_qutrit_code.json`, key `puncture_columns_1indexed`).

```python
import qmsd, qmsd.triorthogonal, qmsd.reedmuller, json
S = json.load(open("cap_qutrit_code.json"))["puncture_columns_1indexed"]
built = qmsd.triorthogonal.build_triorthogonal_code(
            3, 7, 4, S, G=qmsd.reedmuller.rm_generator(4, 7, 3))
G0 = built["X_stab"]          # 55 x 1968 over F_3   (55 = 274-219 rows, 1968 = 2187-219 cols)
```

The quantum distance is

```
d = min weight of G0^perp
  = minimum number of F_3-linearly-dependent columns of G0
  = min_{c in RM_3(9,7), c != 0} |supp(c) \ S|        (G0^perp = RM_3(9,7) punctured at S)
```

Exact global methods are infeasible here: MacWilliams needs `3^55`; meet-in-the-middle needs
`~C(1968,5) ~ 2.5e13` stored syndromes (~200 TB). We therefore used two **feasible**,
**qmsd-cross-validated** validators.

---

## 2. The two validators (both validated against `qmsd` ground truth)

### Method 1 — q-ary (F_3) Stern / ISD low-weight-codeword search on `G0^perp`

`cap_validate/stern_isd.py`, `cap_validate/stern_distribute.py`. Low memory (collision
lists, not the full `C(n,w)` table); multithreaded numba hot loop; deep iteration over
random information sets via F_3 Gaussian elimination with a birthday collision on `l`
parity coordinates. It also runs an **exact O(n) zero-/proportional-column scan**
(`trivial_low_weight`) that is *not* probabilistic.

One-sided guarantee: finding a codeword of weight `w < 10` **PROVES** `d < 10` and emits
the witness; finding nothing below 10 after a deep search is strong probabilistic evidence
(not proof) of `d >= 10`.

**qmsd cross-checks (all PASS, in scope = F_3):**
- `[[206,37,4]]_3`: qmsd exact `d=4` (MacWilliams `B=[1,0,0,0,880,...]`); Stern best `=4`, witness verified -> MATCH.
- `[[667,62,4]]_3`: published `d=4`; Stern found a verified weight-4 codeword -> MATCH.
- **Counterexample mechanism caught:** the m=4, r=3 "2-flat-mirror" 18-cap has qmsd-exact `d=3`
  driven by a *higher*-weight RM codeword dropping below the weight-18 structured class bound (=5).
  Stern (run in the correct `l<r` regime) found the weight-3 witness `[28,49,54]` -> MATCH.
  This is the exact "a high-weight codeword punctures low" mechanism that decides the cap target.
- Out of scope: F_5 anchor `[[112,13,3]]_5` — the kernel hardcodes F_3 (`%3`), so non-F_3 input
  is mangled; not usable as an F_5 validator. (Known defect, documented below.)

**Known defects (do not affect this verdict):**
- Degenerate collision window when `l >= r` (`stern_search_multi` caps `l` at `min(l, rows)` but
  never enforces `l < r`); with `l = r` it silently under-reports. The cap target has `r = 55`,
  default `l = 13`, safely in the `l < r` regime.
- No field guard (F_3 only; asserts nothing on non-F_3 input).
- Probabilistic coverage is tied to info-support `= 2p`; thorough runs must sweep `p` and keep `l < r`.

### Method 2 — rigorous structured enumeration of low-weight RM_3(9,7) codewords

`cap_validate/structured/`. Enumerates structured RM_3(9,7) codewords (unions of `j` parallel
`w`-flats — genuine codewords) and computes `min |supp(c) \ S|` exactly over those families,
via the additive reduction "per direction, top-`j` coset occupancy of `S`" in a numba `prange`
loop.

- `w18` = the **complete** minimum-weight class (two parallel 2-flats): 99,463 directions
  x `C(243,2)` = 2,924,510,589 supports. Verified: numba kernel == brute force on sampled
  directions; an independent pure-numpy reimplementation reproduces every number.
- `w27` = single affine 3-flats (a partial weight-27 subfamily).

**Result:** `w18 -> min 10`, `w27 -> min 18`, so the rigorous minimum over the enumerated
families is **10**.

**Scope / honesty (high-severity caveats):**
- This is a **rigorous lower bound on `d` restricted to the enumerated families only**.
  It is **NOT** a proof that `d >= 10`. As framed ("rigorous for RM-weight <= 45") the
  shipped code under-delivers: only `w18` (complete) and `w27` (partial) are enumerated;
  weights 19-26, 28-45, and all non-flat Leducq second/third-weight classes are deferred to
  Method 1. RM weights `> 45` are entirely out of Method-2 scope by design.
- The m=4 analog proves the gap is real: there, the structured class bound (5) sits strictly
  *above* the true distance (4/3). A "structured min = 10" must never be read as `d >= 10`.

---

## 3. Local results (already established, on this machine)

| validator | quantity | value | nature |
|---|---|---|---|
| M2 structured | min `|supp\S|` over enumerated RM-weight<=45 families (`w18`, `w27`) | **10** | rigorous, but only within that window |
| M1 Stern (exact scan) | `trivial_min_weight` (zero/proportional column scan) | **1** | **rigorous, exact** |
| M1 Stern (probabilistic) | lowest weight found in ~12 s, 80 iters, 20 threads | 2 (support `[1003,1241]`) | probabilistic, witness-verified |

**The decisive, rigorous fact.** `G0` (55 x 1968) has **83 all-zero columns**. A zero column
`j` is, by definition, a single linearly-dependent column — i.e. `G0 @ e_j = 0`, so `e_j` is a
genuine **weight-1** codeword of `G0^perp = ` punctured RM_3(9,7). This was independently
re-verified in this session:

```
G0 shape (55, 1968)
num zero columns of G0: 83
syndrome of e_j all zero? True   =>  d = 1
```

Each such `e_j` corresponds to an RM_3(9,7) codeword (weight 77, 103, 112, 145, ...) whose
support lies *almost entirely inside the cap* `S` (e.g. 76 of 77 points in `S`), leaving a
single point outside. These high-weight, cap-concentrated codewords are exactly what
- Method 2's weight `<= 45` window never reaches,
- the prior "cap meets any 2-flat in `<= 4` points" argument (which only bounds the weight-18
  class) never covered, and
- the infeasible `3^55` MacWilliams / `C(1968,5)` MITM could not enumerate.

Because a weight-1 codeword exists, `d = 1 < 10`, exactly and rigorously. **`gamma >= 1`. The
`gamma < 1` claim is refuted.** No probabilistic Stern run is needed to reach this; the deep
distributed pass (below) only re-confirms it and explores the full low-weight spectrum.

---

## 4. What a 'refuted' vs 'supported' outcome means for `gamma < 1`

- **Refuted** (this case): a verified codeword of weight `< 10` exists. This is a *proof*
  that `d < 10`, hence `gamma >= 1`. The witness (here: any of the 83 zero columns, or the
  Stern witness `[1003,1241]`) is checkable in O(rows) by confirming `G0 @ witness = 0`.
  The cap qutrit code does **not** deliver `gamma < 1`.
- **Supported-so-far** (the outcome we would have reported had nothing `< 10` been found):
  Method 2 gives a *rigorous* lower bound on `d` over RM-weight `<= 45`, and a deep Method-1
  sweep finds nothing below 10. That is **high-confidence validation, not a proof** of
  `d >= 10`: Stern's "found nothing" is one-sided probabilistic evidence, and Method 2's
  rigor does not extend past its enumerated families. A genuine proof of `d >= 10` would
  require either completing the rigorous structured enumeration through *all* RM weights or
  an exact global method (both infeasible at these sizes).

For this code we are in the first bucket, and the evidence is the strongest kind: an exact,
re-verifiable weight-1 witness.

---

## 5. Honest bottom line

The validation is **conclusive and rigorous**, and it is **negative**: `d = 1`, far below 10,
so `gamma >= 1` and the `gamma < 1` claim is **refuted**. The refutation does not rely on the
probabilistic Stern search at all — it follows from an exact O(n) structural fact (83 zero
columns of `G0`). The two-machine deep pass described in `CAP_VALIDATE_RUNBOOK.md` is built,
self-tested, and ready; for this code it serves only to (a) reproduce the refutation
independently and (b) map the full low-weight spectrum, **not** to change the verdict.
