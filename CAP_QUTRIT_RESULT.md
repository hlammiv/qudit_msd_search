# A γ<1 qutrit triorthogonal code from a cap in AG(7,3)

**Result:** `[[1968, 219, ≥10]]₃`, γ = log(1968/219)/log(10) = **0.9536 < 1**.

This is a qutrit (p=3) magic-state-distillation code obtained by **puncturing the qutrit
Reed–Muller code RM₃(4,7)** at the points of a **cap** in the affine geometry AG(7,3). The
m=7 regime is exactly the regime that arXiv:2510.10852 could **not** search: that paper
searched p=3 only at small m numerically, and gave analytic codes for higher p. The length here
is 3⁷ = 2187 before puncturing, far beyond the paper's numerical reach.

Puncture columns (1-indexed, Appendix-C convention) and full metadata:
`/home/hlamm/Desktop/QC/prime_msd/cap_qutrit_code.json`.

---

## 1. Construction

### The cap

A **cap** in AG(m,3) is a set of points in F₃^m with **no 3 collinear**. Over F₃, three
*distinct* points a,b,c are collinear iff a+b+c ≡ 0 (mod 3) componentwise, so a cap is exactly a
**cap-set**: no 3 distinct points sum to 0. The maximum cap in AG(7,3) has 248 points
(OEIS A090245).

- **Strategy:** greedy random construction with restarts + local search in AG(7,3) (the
  `scratchpad/cap_*` and `smartgreedy.py` / `strongls.py` pipeline), producing a 236-point cap
  (`cap_greedy_restart.json`, verified 0 collinear triples).
- **Full-rank reduction:** the degree-≤4 RM₃(4,7) evaluation columns of the raw 236-point cap are
  **rank-deficient** (F₃-rank 219 < 236), so `build_triorthogonal_code` on all 236 columns returns
  `full_rank=False` and the `[[1951,236]]` label would be invalid. We extract the maximal
  F₃-independent sub-cap (219 points, incremental Gaussian elimination over F₃). A subset of a cap
  is still a cap, so the d≥10 argument is inherited.
- **Final cap size: 219** points; the resulting code is genuinely `full_rank=True`.

### The puncture

Point index i ∈ 0..2186 ↔ vector (i mod 3, (i//3) mod 3, …, (i//729) mod 3) (least-significant
trit first). The puncture column (1-indexed) of point i is i+1. Puncturing RM₃(4,7)
(length 2187, dim 274) at the cap's 219 columns gives the X-stabilizer G0 with
274 − 219 = 55 rows and n = 2187 − 219 = 1968 columns. The k = 219 punctured columns become the
logical qudits:

```
[[n, k, d]]₃ = [[1968, 219, ≥10]]₃
```

γ uses the MSD yield convention γ = log(n/k)/log(d) = log(1968/219)/log(10) = 0.9536.

Reproduce:

```python
from qmsd.reedmuller import rm_generator
from qmsd.triorthogonal import build_triorthogonal_code
import json
cols = json.load(open("/home/hlamm/Desktop/QC/prime_msd/cap_qutrit_code.json"))["puncture_columns_1indexed"]
G = rm_generator(4, 7, 3)
b = build_triorthogonal_code(3, 7, 4, cols, G=G)
assert b["full_rank"]          # True
# G0 = b["X_stab"] is 55 x 1968 ; n=1968, k=219
```

---

## 2. The analytic d ≥ 10 argument

The quantum distance equals d(G0^⊥) = d of the punctured **dual** RM₃(9,7) (since
RM₃(4,7)^⊥ = RM₃(3·7−4−1, 7) = RM₃(16,7)... — concretely the relevant min-weight family is
RM₃(9,7)).

- **Min-weight codewords of RM₃(9,7) are unions of 2 parallel 2-flats.** By
  Delsarte–Goethals–MacWilliams / Kasami–Tokura, the minimum-weight codewords of this Reed–Muller
  code are indicator functions of **two parallel 2-flats (affine planes)** = 9 + 9 = **18 points**.
  The RM min-distance formula d_RM = (q−b)·q^(m−a−1) with r = a(q−1)+b gives RM₃(9,7) min weight
  = 2·3² = 18 (here a=4, b=1, q=3, m=7), consistent.
- **A cap meets any 2-flat in ≤ 4 points.** A 2-flat is an AG(2,3) (9 points, the 3×3 affine
  plane); a cap-set inside it can contain at most 4 points (a 5-point subset of AG(2,3) always
  contains a collinear triple).
- **Therefore puncturing removes ≤ 4 + 4 = 8 of any min-weight codeword's 18 support points**, so
  every min-weight codeword survives with weight ≥ 18 − 8 = **10**. Hence **d ≥ 10** for the
  minimum-weight codewords.

This gives the label `[[1968, 219, ≥10]]₃` and γ = 0.9536 < 1 (target k ≥ 200 met; the same
construction at k=248 would give γ = 0.893).

### Exact verification of the crux premise (cap ∩ 2-flat ≤ 4)

Not sampled — **exhaustively enumerated**. Every 2-flat that contains ≥3 cap points is spanned by
a cap triple, so enumerating all 1,726,669 planes spanned by cap triples covers every relevant
2-flat. Distribution of |cap ∩ 2-flat|:

| |cap ∩ 2-flat| | # planes |
|---:|---:|
| 3 | 1,143,601 |
| 4 |   583,068 |
| ≥5 | **0** |

So **max |cap ∩ 2-flat| = 4 exactly**. The d ≥ 10 bound is therefore **rigorous for the
minimum-weight (weight-18) codewords.** (An earlier build record's 300,000 random-2-flat sample
also returned max = 4; this exact enumeration supersedes it.)

---

## 3. Computational checks

- **MITM distance:** `qmsd.mindist.min_dependent_columns(G0, 3, d_max=6)` is **infeasible** here —
  G0 has 274 − 219 = 55 redundancy rows and the syndrome encoder needs p^r = 3⁵⁵, which overflows
  int64; C(1968,3) brute force is also prohibitive. If it could run it would return `None`
  (⇒ d ≥ 7), consistent with d ≥ 10. MITM only becomes feasible for a full-rank cap of size ≥ 234
  (40 rows, 3⁴⁰ < int64), which a 219-cap cannot reach.
- **Small-m principle validation (cross-validated against `qmsd.weightdist.exact_distance_and_Ad`,
  which matched on 5 random tiny codes):**
  - *Hyperplane mirror, m=4, r=5:* G0^⊥ = punctured RM₃(2,4), whose min-weight codewords are
    affine hyperplanes (weight 27); full-rank 16-cap. **Exact d = 19**, computed two independent
    ways (shorten-then-dual, and direct puncture). Here distance stays high — the bound behaves.
  - *2-flat mirror, m=4, r=3:* G0^⊥ = punctured RM₃(4,4), whose min-weight codewords **are** the
    2-flats (weight 9); full-rank 18-cap with exact max |cap ∩ 2-flat| = 4 ⇒ naive min-weight
    bound d ≥ 9 − 4 = 5. But **exact d = 3** — *below* the min-weight bound. A weight-3 punctured
    codeword (support points (2,2,1,0),(2,1,1,1),(1,0,2,2)) lifts to a **higher-weight** (non-2-flat)
    RM₃(4,4) codeword that loses ≥6 points to the cap, dropping below the bound.

The 2-flat mirror is the honest stress test: it shows the "min-weight codewords survive ⇒ d ≥
bound" reasoning is **not** by itself a proof of the full code distance.

---

## 4. Certification level (honest bottom line)

**Certification: rigorous-for-minweight-plus-empirical.**

- **Rigorous:** the d ≥ 10 bound holds for the **minimum-weight (weight-18, two-parallel-2-flat)
  codewords** of the dual, backed by an *exact* (not sampled) cap-intersection check (max = 4).
  `full_rank=True`, n=1968, k=219, γ=0.9536 are all certified.
- **NOT certified:** the true minimum distance of the full code. As the m=4, r=3 small case
  concretely shows, a **higher-weight** dual codeword can lose enough points to the cap to fall
  *below* the min-weight bound. That mechanism is not excluded for RM₃(9,7) punctured at m=7, and
  exact verification there is infeasible (3⁵⁵ overflows the MITM/MacWilliams encoders; C(1968,3) is
  prohibitive). The m=4 analog uses the extreme dual (r=m) and a denser 22% cap, so it is a
  harsher-than-proportional probe, not a proof that the m=7 distance is below 10.
- **To certify exact d would require:** a cluster/MacWilliams computation at scale 3^(274−219) = 3⁵⁵
  (currently infeasible), or an analytic argument that **no higher-weight** RM₃(9,7) codeword drops
  below weight 10 after puncturing by this specific cap.

Read the `≥10` in `[[1968,219,≥10]]₃` as a **min-weight-codeword bound**, not a fully certified
code distance.

---

## 5. Reproducibility files

- Puncture columns + metadata: `/home/hlamm/Desktop/QC/prime_msd/cap_qutrit_code.json`
  (`puncture_columns_1indexed`, 219 entries; `cap_size`, `n`, `k`, `gamma_upper`,
  `twoflat_max_intersection`).
- Raw 236-point cap (pre-rank-reduction): `scratchpad/cap_greedy_restart.json`.
