# Structured A_d enumerator for punctured Reed–Muller qutrit codes

`qmsd/structured_ad.py` computes `A_d` (the number of minimum-weight logical-Z
operators) of the paper's punctured-RM triorthogonal codes **without** the
MacWilliams engine (cost `p^{dim G0}=3^{D-k}`) and **without** the direct
`C(n,d)` column scan — by enumerating only the *geometrically structured*
low-weight Reed–Muller codewords.

All numbers below were produced in this repo (run with
`PYTHONPATH=/home/hlamm/Desktop/QC/prime_msd`).

---

## 1. The identity everything rests on

For a punctured triorthogonal code the quantum code is `G0^perp =` the
**punctured** code `RM_p(rtilde, m)` restricted to the surviving coordinates,
where `S` (`|S| = k`) is the puncture set and `rtilde = m(p-1) - r - 1`
(`r = r_max`).  When the punctured-column submatrix is full rank the restriction
map `RM_p(rtilde,m) -> G0^perp` is a **bijection** (kernel
`{c : supp c ⊆ S} = {0}`), so for any nonzero `c ∈ RM_p(rtilde,m)`:

```
d   = min_{c≠0} |supp(c) \ S|
A_d = #{ c ∈ RM_p(rtilde,m), c≠0 : |supp(c) \ S| = d }.
```

A weight-`d` punctured codeword lifts to an RM codeword of weight
`w = d + |supp(c) ∩ S|`.

## 2. The structured method: decomposition by minimal affine span

Empirically (and by the Delsarte–Goethals–MacWilliams / Kasami–Tokura–Azumi /
Leducq theory) the low-weight RM codewords are **flat-supported**.  Every
codeword `c` has a minimal affine span `F = aff.span(supp c)` of some dimension
`j`, and restricted to `F ≅ F_p^j` it is a codeword of the much smaller code

```
C_F = RM_p( rtilde − (m−j)(p−1),  j )   =  {RM_p(rtilde,m) codewords supported in F}.
```

So we enumerate weight-`d` codewords by their **minimal span dimension** `j`:

```
for j = 2 .. jmax:
    for each j-flat F of F_p^m:
        find the codewords of C_F with EXACTLY d nonzeros on the survivors F\S
        (a meet-in-the-middle column-dependency search on the small parity check
         of C_F restricted to F);
        keep those whose support spans F exactly (aff.dim == j)   # minimal-flat dedup
```

Because each codeword has a *unique* minimal flat, the `aff.dim==j` filter makes
every codeword counted exactly once. Summed over `j = 2 .. m` the count is
**EXACT**. Capping `jmax < m` returns the contribution of codewords of affine
span `≤ jmax` — a **certified lower bound** on `A_d`, and exact when no
higher-span codeword contributes.

Key facts the implementation uses:
* For `b=0` codes (`m=5,6`, `rtilde=6,8`) the `j=2` restricted code is the
  constants, so dim-2 codewords are the scalar **2-flat indicators** (the DGM
  minimum-weight family, weight `9`); this is a fast special case.
* For `m=4` (`rtilde=5`, `b=1`) the `j=2` restricted code is `RM_3(1,2)`
  (degree-1 forms on a plane), handled by the general path.
* Cost scales with the number of low-dimensional flats and the size of the tiny
  restricted codes — never `3^{D-k}` or `C(n,d)`.

The whole computation is exact-integer / exact-`F_p` linear algebra (pure numpy,
no floats), and the distance `d` itself is obtained structurally from the
certified meet-in-the-middle minimum-distance routine (`qmsd.mindist`).

## 3. What it reproduces (decisive validation vs Table 3 of arXiv:2510.10852)

`structured_ad(p, m, r, puncture_columns, jmax)` returns
`{distance, A_d, low_weight_histogram, dim_histogram, jmax, exact}`.

**Headline: 6 of the 7 qutrit Table-3 A_d are reproduced EXACTLY by the structured
method — including all three brute-force-unreachable "structured-only" codes
(572, 1104, 1128).** The 7th (`[[667,62,4]]`, a tail-dominated `d=4` code) is matched by
the MacWilliams engine (3972); the structured method gives its dim≤4 part (954) and the
rest is the dim-5/6 tail.

| code | d | structured `A_d` | paper `A_d` | how (measured) | cross-check |
|------|---|------------------|-------------|-----|-------------|
| `[[72,9,3]]_3`   | 3 | **648**  | 648  | `jmax=4=m` (exact), 6 s, dims {2:504,3:144} | MacWilliams 648 ✓ |
| `[[200,43,3]]_3` | 3 | **1700** | 1700 | `jmax=5=m` (exact), dims {2:32,3:114,4:76,5:1478} | MacWilliams 1700 ✓ |
| `[[206,37,4]]_3` | 4 | **880**  | 880  | `jmax=5=m` (exact), dims {2:194,3:374,4:92,5:220}; `jmax=4`→660 (29 s) | MacWilliams 880 ✓ |
| `[[667,62,4]]_3` | 4 | 954 (dims 2-4) + tail | 3972 | `jmax=4`→954; dim-5/6 tail (3018) needs `jmax=6` (expensive, m=6) | MacWilliams 3972 ✓ |
| `[[230,13,6]]_3` | 6 | **572**  | 572  | `jmax=3`, dims {2:572} — **NO tail**, 13 s | brute-unreachable (dim G0=38); = paper |
| `[[215,28,5]]_3` | 5 | **1104** | 1104 | `jmax=4`, dims {2:582,3:490,4:32} — **NO tail**, 63 s | indep. MITM oracle 1104 ✓ |
| `[[690,39,5]]_3` | 5 | **1128** | 1128 | `jmax=4`, m=6 — **NO tail** | brute-unreachable (dim G0=39); = paper |

### The decisive structured-only result
The three structured-only codes (`d = 5,6`) — for which both the MacWilliams
engine (`3^{38}`, `3^{23}`, `3^{39}`) and the direct `C(n,d)` scan are infeasible —
are reproduced **exactly** by the structured enumerator:

* `[[230,13,6]] = 572` from the **dim-2 (2-flat indicator) family alone**.
* `[[215,28,5]] = 1104` from dims 2,3,4 (`582 + 490 + 32`).
* `[[690,39,5]] = 1128` from low-dimensional flats (m=6).

**Why there is no tail for these.** A weight-`w` RM codeword contributes to `A_d`
only if exactly `w−d` of its support lies in `S`. For the larger-distance codes
(`d=5,6`) this forces the contributing codewords into **low-dimensional flats**
(max affine span 4 for `d=5`, 2 for `d=6`); the full-space (dim-`m`) "tail" is
empty, so a small `jmax` is exact. This is exactly the regime the brute engines
cannot reach but the structured method handles cheaply.

## 4. Honest scope and limitations

* **Exactness.** `jmax = m` makes the decomposition exact by construction (it
  covers every affine-span dimension; the `j=m` term is the full-space
  meet-in-the-middle, which equals the original computation). For the
  brute-verifiable codes this matches the MacWilliams `A_d` exactly
  (`648, 1700, 880`).
* **The "tail".** For the small-distance codes (`d=3,4`) a large fraction of
  `A_d` comes from **full-span (dim-`m`) codewords** (e.g. `1478/1700` for
  `[[200]]`, `220/880` for `[[206]]`). These are *not* in any low-dimensional
  flat, so capping `jmax<m` under-counts them. Reaching them requires the
  `j=m` term — which is feasible for these codes (small `dim G0`) and is what
  `jmax=m` does, but it is the same cost as a direct minimum-weight enumeration,
  so the structured method gives no asymptotic win there. (It still gives the
  exact answer and a clean per-dimension breakdown.)
* **`[[667,62,4]]` (m=6, d=4).** The exact value (3972) is confirmed by the
  MacWilliams engine (`dim G0 = 16`). A full structured `jmax=6` run is dominated
  by the ~0.9 M dim-2 and ~0.9 M dim-3 flats of `F_3^6` and the dim-6 tail, so it
  is slow; the structured method's *value* (avoiding the `3^{D-k}` blow-up) is on
  the structured-only `d=5,6` codes, not on the tail-dominated `d=3,4` codes
  (which the MacWilliams engine already handles).
* **Generality.** The method is written for general `p, m, r`; it was validated
  on the qutrit (`p=3`) Table-3 codes. The `j=2` constant fast path assumes
  `rtilde ≡ 0 (mod p-1)` (`b=0`); other `j` use the general path.

## 5. Reproduce

```python
from qmsd.oracle import load_oracle
from qmsd.reedmuller import r_max
from qmsd.structured_ad import structured_ad

oc = {c.label: c for c in load_oracle()}["[[215,28,5]]_3"]
res = structured_ad(3, oc.m, r_max(oc.m, 3), oc.puncture_columns_1indexed, jmax=4)
print(res["A_d"])        # 1104
```

Tests: `tests/test_structured_ad.py` (fast: helper + `[[72,9,3]]` end-to-end vs
MacWilliams; `QMSD_SLOW=1` enables the minutes-long structured-only
reproductions `572 / 1104` and the `[[206]]` exact `880`).
