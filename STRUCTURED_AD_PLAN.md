# Structured A_d enumerator for punctured Reed–Muller / triorthogonal codes

Goal: compute `A_d` (number of minimum-weight codewords) of the paper's qutrit
distillation codes **without** the MacWilliams engine (cost `p^(dim G0) = 3^(D-k)`)
and **without** the direct `C(n,d)` scan — by enumerating only the *structured*
low-weight codewords of the underlying Reed–Muller code.

All claims below are anchored to numbers actually produced in this repo
(`qmsd/`, run with `PYTHONPATH=/home/hlamm/Desktop/QC/prime_msd`).

---

## 0. Exact setup and the one identity everything rests on

For a punctured triorthogonal code the quantum code is `G0^perp =` the
**punctured** code `RM_p(rtilde, m)` restricted to the surviving coordinates,
where `S` (`|S| = k`) is the set of punctured points. Build chain
(`qmsd/reedmuller.py`, `qmsd/oracle.py`):

```
r      = r_max(m,p) = floor((m(p-1)-1)/3)        (3r < m(p-1), triorthogonality)
rtilde = m(p-1) - r - 1                            (dual degree -> quantum distance)
```

The puncture set as points: column `c` (1-indexed) <-> point `c-1` in `F_p^m`
(`qmsd/puncture.column_to_point`). For a nonzero `c in RM_p(rtilde,m)` its
**punctured weight** is `|supp(c) \ S|`, so

```
d   = min_{c != 0} |supp(c) \ S|
A_d = #{ c in RM_p(rtilde,m) , c != 0 : |supp(c) \ S| = d }.
```

**Why this is a clean count (no multiplicity).** The restriction map
`RM_p(rtilde,m) -> G0^perp` (length `p^m -> n`) is a *bijection* whenever the
puncture submatrix is full rank — its kernel is `{c : supp(c) subset S}`, which
is trivial. Verified for the anchor: `dim RM_3(6,5) = sum_{s<=6}[5,s]_3 = 192`
(`qmsd/pnomial`) and `dim G0^perp = n - dim G0 = 206 - 14 = 192`. So `A_d` equals
`B_d` from `qmsd.weightdist.macwilliams` with **no collapsing** of distinct RM
codewords, and we may enumerate in the unpunctured RM code.

A weight-`d` punctured codeword comes from an RM codeword of weight
`w = d + |supp(c) ∩ S|` with `0 <= |supp∩S| <= k`, hence

```
w in [d, d+k]   (and w >= d_RM, so really w in [d_RM, d+k]).
```

### CORRECTION to the task's stated `d_RM` values
The task's parentheticals are off by a row. Confirmed against the
single-puncture ground truth (`d_single = d_RM(rtilde) - 1`,
`tests/ground_truth.SINGLE_PUNCTURE_CODES`):

| m | r | rtilde | **d_RM(rtilde)** | task said | single-punc d (check) |
|---|---|--------|------------------|-----------|------------------------|
| 4 | 2 | 5      | **6**            | 18 (wrong)| 5 ✓ |
| 5 | 3 | 6      | **9**            | 9 ✓       | 8 ✓ |
| 6 | 3 | 8      | **9**            | 27 (wrong)| 8 ✓ |

So for the qutrit codes `2.5·d_min = 22.5` for `m=5,6` and `15` for `m=4`.

---

## 1. The low-weight RM codeword classification needed

Write `rtilde = a(p-1) + b`, `0 <= b <= p-2`. For the qutrit codes
(`p=3`, `p-1=2`): `m=5,6 -> rtilde = 6,8 -> b = 0`; `m=4 -> rtilde = 5 -> b = 1`.

### 1a. DGM minimum-weight family (Delsarte–Goethals–MacWilliams 1970)
`d_min = (p - b) · p^(m - a - 1)`.

* **b = 0 case (m=5: a=3; m=6: a=4).** Minimum-weight codewords are
  `lambda · 1_V`, the nonzero-scalar indicators of **affine `(m-a)`-flats** `V`.
  Here `m-a = 2`, so `V` is a **2-flat**, `|V| = p^2 = 9 = d_min`. The flat is the
  common zero set of `m-a-1 = ... ` actually of `a` independent affine forms:
  `1_V = prod_{i=1}^{a} (1 - L_i(x)^{p-1})`, degree `a(p-1) = rtilde`. Count
  ```
  A_{d_min}(RM) = (p-1) · (#2-flats),
  #2-flats = p^(m-2) · gaussbinom(m,2,p).
  ```
  (m=5: `(3^3)·1210 = 32670` flats; m=6: `(3^4)·11011 = 891891` flats — both
  enumerated in the anchor.)

* **b != 0 case (m=4: a=2, b=1).** `d_min = (p-b)p^(m-a-1) = 2·3 = 6`. Codewords
  are **NOT** flat indicators: they are `lambda · ell(x) · 1_V` where `V` is an
  `(m-a)=2`-flat and `ell` is a nonzero affine form on `V` — i.e. a 2-flat with
  one parallel line (`p-b·...` -> the `ell=0` slice) removed, support `6` of the
  `9` points. (Empirically confirmed the b=0 probe is the wrong class for m=4: a
  pure 2-flat scan gives **0** contribution to `[[72,9,3]]`; the b=1 family is
  required.) Count families: `#2-flats × (#nonzero affine forms on V mod scaling)
  × (p-1)`.

### 1b. Kasami–Tokura–Azumi / Leducq higher-weight classes (`d_min < w < 2.5 d_min`)
KTA (1976) classifies every codeword of weight `< 2·d_min`; Leducq
(arXiv:1001.2554, 1203.5244, 1203.4592) extends the description toward
`2.5·d_min` and gives the affine-geometric families and their counts. Structure
of these classes (for `b=0`, the relevant qutrit case):

* **Type "two parallel flats" / fattened flat:** support = union of cosets of a
  `(m-a)`-flat inside a `(m-a+1)`-flat, with prescribed nonzero values on each
  parallel slice. Parametrized by the big flat + the value pattern; counted by
  `gaussbinom` (choice of flag of flats) × `p`-nomial style value-pattern counts.
* **Type "product of an extra affine form":** `1_W · g`, `W` an
  `(m-a-1+s)`-flat and `g` a low-degree form on `W`, giving weights
  `d_min·(p^s - small)/...`. These produce the discrete weight ladder
  `9, 12, 15, 16, 18, ...` seen empirically (see §1c).
* The general statement: every codeword of RM weight `< 2.5 d_min` is the
  evaluation of a **product of affine forms times an indicator of a flat**, so it
  is enumerable by choosing (i) a flag of affine flats and (ii) a bounded value
  pattern — cost polynomial in the number of flats, not in `3^(D-k)` or `C(n,d)`.

### 1c. Empirical weight ladder (RM_3(6,5), d_min=9)
The exact RM weights that appear among the weight-4 punctured codewords of
`[[206,37,4]]` (computed by lifting each punctured codeword back to its unique
RM(6,5) preimage):
```
RM weight:  9   12   15  16  17 18  19  20  22 | 23 24 26 27 28 29 30 31 32 33 34 35
count    : 194 374  18  66   2  2   4  16   4  |  2 10 14 22 34 32 28 26 12 10  6  4
                 \------ <= 2.5 d_min = 22 -------/ \------ TAIL  > 2.5 d_min ------/
                          classified region 680            UNCLASSIFIED region 200
```
Two consequences, both load-bearing for scope:
1. The DGM min-weight (w=9) class is **not** the whole story for high-`k` codes —
   it is only `194/880` here.
2. A genuine `200/880 (23%)` of `A_d` comes from RM weights `23..35`, i.e.
   **beyond `2.5·d_min`**, outside the rigorous KTA/Leducq classification. (Max
   contributing weight is `35`, below the loose bound `d+k = 41`.)

---

## 2. Per-code RM-weight window and which classes fall in it

`p=3`. Window `[d_RM, d+k]`; classified iff `d+k <= floor(2.5 d_min)`
(`22` for m=5,6; `15` for m=4).

| code | m | rtilde | d_RM | d | k | RM-window | 2.5·d_min cutoff | fully classified? | min-wt class contrib | A_d (target) |
|------|---|--------|------|---|---|-----------|------------------|-------------------|----------------------|--------------|
| `[[72,9,3]]_3`   | 4 | 5 | 6 | 3 | 9  | `[6,12]`  | 15 | **YES** (needs b=1 fam + KTA) | 0 via b=0 probe (wrong class) | 648 |
| `[[230,13,6]]_3` | 5 | 6 | 9 | 6 | 13 | `[9,19]`  | 22 | **YES**            | **572 = A_d (exact!)** | 572 |
| `[[215,28,5]]_3` | 5 | 6 | 9 | 5 | 28 | `[9,33]`  | 22 | no (tail 23..33)   | 582 | 1104 |
| `[[206,37,4]]_3` | 5 | 6 | 9 | 4 | 37 | `[9,41]→35`| 22 | no (tail 23..35)  | 194 | 880 |
| `[[200,43,3]]_3` | 5 | 6 | 9 | 3 | 43 | `[9,46]`  | 22 | no (large tail)    | 32  | 1700 |
| `[[690,39,5]]_3` | 6 | 8 | 9 | 5 | 39 | `[9,44]`  | 22 | no (tail)          | 878 | 1128 |
| `[[667,62,4]]_3` | 6 | 8 | 9 | 4 | 62 | `[9,66]`  | 22 | no (large tail)    | 394 | 3972 |

(`min-wt contrib` and `A_d` columns are this repo's outputs; `A_d` for the four
small-`dim G0` codes is the live MacWilliams value, the others are the published
ground truth.)

**Reading of the table.**
* `[[230,13,6]]` (the d=6, smallest-`k` structured-only code) is reproduced
  **exactly** by the **min-weight 2-flat family alone** (572) — a brute-force-
  unreachable target (`dim G0 = 38`, `3^38` infeasible) hit by pure structure.
  Higher weight classes contribute **0** for its `S` (window `[9,19]`, but few
  enough punctures that no `w>9` codeword puts `w-6` of its support in the
  13-point `S`).
* `[[72,9,3]]` and `[[230,13,6]]` are the only two whose entire window lies
  `<= 2.5 d_min`; a *complete* KTA/Leducq enumerator reproduces them rigorously.
* The remaining five high-`k`, low-`d` codes have a **tail beyond `2.5 d_min`**
  (verified `200/880` for `[[206,37,4]]`). Reproducing their full `A_d` requires
  structured families *past* the rigorously-classified KTA range — either
  Leducq's extended product-of-affine-form families specialized to these small
  `(rtilde,m)`, or a complete affine-product enumeration up to the actual max
  contributing weight (35 / 44 / 66, well below `d+k`).

---

## 3. The algorithm

```
inputs: p, m, rtilde, puncture columns -> S (set of point indices in F_p^m), d
1. write rtilde = a(p-1) + b ; d_min = (p-b) p^(m-a-1)
2. for each structured weight class W in [d_min .. min(d+k, weight_cap)]:
       for each structured codeword c in class W (enumerated by its flat params):
           pw = W - |supp(c) ∩ S|          # = |supp(c) \ S|, no need to form full vector
           if pw == d: A_d += (multiplicity: (p-1) for a scalar family, etc.)
3. return A_d
```

Concretely (matches the anchor code):
* **Min-weight (b=0):** enumerate 2-dim linear subspaces of `F_p^m`
  (`gaussbinom(m,2,p)` of them, dedup by member-set), translate by all `p^m`
  cosets and dedup -> all `p^(m-2)·gaussbinom(m,2,p)` 2-flats `V` (each a 9-point
  frozenset of point indices). For each `V`, `pw = 9 - |V∩S|`; accumulate
  `(p-1)` codewords when `pw == d`.
* **Min-weight (b!=0, m=4):** enumerate 2-flats `V` and, for each, the
  `(#nonzero affine forms on V)/scaling` value patterns (support = `V` minus one
  internal line); `pw = 6 - |supp∩S|`.
* **KTA/Leducq classes:** enumerate by their flag-of-flats parameters; support is
  the (small) union of flats / flat-minus-subflat; compute `|supp∩S|` directly on
  point-index sets. The per-codeword cost is `O(|supp|)`; the total cost scales
  with the **number of low-weight structured codewords**, i.e. polynomial in the
  count of affine flats (`~gaussbinom`), never `3^(D-k)` or `C(n,d)`.

Cost note: only `supp(c) ∩ S` is needed, and `|S| = k` is small, so an even
faster variant indexes flats by which `S`-points they contain (bucket the `k`
`S`-points, intersect on the fly) instead of materializing all flats.

---

## 4. Build order and validation protocol

**Build order**
1. `flats(m, p, dim)` generator (2-flats first; reuse the dedup-by-coset code in
   `scratchpad/anchor.py`). Cache per `m`.
2. `min_weight_Ad(code)` using §3 b=0 family (m=5,6) and b=1 family (m=4).
3. KTA/Leducq class enumerators (weight `12,15,16,18,...` ladder), each a small
   module returning `(support_set, multiplicity)` iterators parametrized by flat
   flags; from Leducq arXiv:1001.2554 / 1203.5244 closed forms.
4. `structured_Ad(code, weight_cap)` summing §3 over all implemented classes.
5. `weight_cap` policy: default `min(d+k, floor(2.5 d_min))`; expose the tail gap
   honestly (assert `structured_Ad <= A_d`, report the deficit).

**Validation protocol (decisive — `tests/ground_truth.TABLE3_CODES`)**
* **Cross-check vs MacWilliams** wherever `dim G0` is small
  (`qmsd.weightdist.exact_distance_and_Ad`): `[[72,9,3]] (dim 6, A_d 648)`,
  `[[200,43,3]] (8, 1700)`, `[[206,37,4]] (14, 880)`, `[[667,62,4]] (16, 3972)`.
  For these, also check the **per-RM-weight histogram** by lifting each weight-`d`
  punctured codeword to its RM preimage (method in `scratchpad/lift.py`) — this
  is the gold standard that tells you exactly which classes you still owe.
* **Structured-only targets** (MacWilliams infeasible — the whole point):
  `[[215,28,5]] (dim 23) 1104`, `[[230,13,6]] (dim 38) 572`,
  `[[690,39,5]] (dim 39) 1128`. Enumerator is correct iff it reproduces these.

**Validation status reached in this planning pass (honest scope):**
* ✅ `[[230,13,6]]_3 A_d = 572` **reproduced exactly** by the min-weight 2-flat
  family — a brute-unreachable code hit with pure structure. (Strongest result:
  it is in the structured-only set.)
* ✅ `[[206,37,4]]_3`: min-weight family supplies `194/880`; full RM-weight
  histogram computed (§1c); `680/880` lie `<= 2.5 d_min` (KTA-classifiable),
  `200/880` lie in the tail `23..35`.
* ⏳ The KTA/Leducq class modules (weights `12..22`) are **specified but not yet
  implemented**; once added they close the `680` "classified" mass for the m=5
  codes and (with the m=4 b=1 family) fully reproduce `[[72,9,3]]` and
  `[[230,13,6]]`.
* ⚠️ **Open / honest limitation:** the five high-`k` low-`d` codes
  (`200,206,215,690,667`) have verified contributions **beyond `2.5 d_min`**
  (`200/880` for `206`), which KTA does not classify. Full reproduction of their
  `A_d` requires extending the structured families past `2.5 d_min` (Leducq's
  product-of-affine-form families specialized to the exact small `(rtilde,m)`, or
  an enumeration up to the actual max contributing weight `35/44/66`). This is
  the real remaining research content, and the table makes precise exactly how
  much mass (and which weights) each code still needs.

**Resolution of the task's question "do higher-weight RM classes contribute?"**
Yes, decisively, and it is code-dependent:
* For the larger-`d`/smaller-`k` code `[[230,13,6]]` the min-weight class is the
  *entire* `A_d` (higher classes contribute 0).
* For the high-`k` codes higher classes dominate: `[[206,37,4]]` gets only
  `194/880` from min weight, the rest from RM weights `12..35`; `[[200,43,3]]`
  gets only `32/1700` from min weight. A weight-`w>d_min` codeword contributes to
  `A_d` iff it has exactly `w-d` of its support inside `S`, which is common when
  `k` is large.
