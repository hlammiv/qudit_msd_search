# prime_msd — Implementation Blueprint

**Goal:** a Python toolkit that *discovers new qudit triorthogonal codes for magic
state distillation with low overhead* (yield `gamma = log(n/k)/log(d) < 1`), grounded
in arXiv:2510.10852 and its lineage.

This document is the **software design**. The **math-to-code mapping** (every formula,
theorem, and pitfall) lives in the project's verified notes on arXiv:2510.10852 (kept
locally), §10 ("Reproduction / Discovery-Tool Checklist") — read that alongside this.

---

## 0. Strategy

The paper hands us two complementary engines and a ready-made correctness oracle:

1. **Analytic engine** (Theorems 4–5). For Manhattan-weight puncturing of maximal-degree
   RM codes, the parameters `[[n,k,d]]` are *pure integer combinatorics* on p-nomial
   coefficients — **no matrices**. This scales to the astronomically large codes in
   Table 2 (e.g. `[[2.9e17, 1.4e13, 21700]]_2`) instantly. Use it to map the whole
   `(p, m, w)` landscape and to reproduce Tables 1–2.

2. **Explicit search engine** (Section 5). Build the actual RM generator over `F_p`,
   puncture an arbitrary column set `S`, and compute `[[n,k,d]]`, `A_d`, `gamma`, `C`
   directly. This is where *new* codes are found (better puncture locations than the
   analytic scheme), but it is bounded by `p^m` (matrix size) and by minimum-distance /
   `A_d` computation cost.

3. **Correctness oracle.** `qmsd/data/puncture_locations.json` holds the
   exact puncture columns for all 10 search codes (already validated: `n+k=p^m`,
   `#punctures=k`). **Before searching for anything new, the explicit engine must
   reconstruct these 10 codes from their puncture sets and reproduce the paper's
   `[[n,k,d]]`, `A_d`, and `gamma` in Table 3.** That is the acceptance test for the
   whole stack.

**Key design consequence:** keep the analytic path (big integers, exact) and the
explicit path (finite-field matrices) as separate layers that *cross-validate* each
other on the Manhattan family.

---

## 1. Dependencies (all present in this environment)

| Lib | Version | Use |
|-----|---------|-----|
| `galois` | 0.4.6 | `GF(p)` arrays + finite-field linear algebra (rank, null_space, row reduce) — the backbone of RM generators, shorten/puncture, and CSS duals. |
| `numpy` | 2.2.6 | array plumbing; integer combinatorics; coordinate enumeration. |
| `sympy` | 1.14.0 | exact rationals/symbols for asymptotics cross-checks; multinomial sanity. |
| `scipy` | 1.16.3 | `brentq` root-find for the saddle equation; special functions. |
| `pytest` (add) | — | the validation harness. |

Python integers are arbitrary precision, so the analytic engine needs no bignum lib.

---

## 2. Package layout

```
prime_msd/
├── IMPLEMENTATION_BLUEPRINT.md        # this file
├── (arXiv paper + notes + tutorial: kept local, not in repo; oracle data is in qmsd/data/)
├── qmsd/                              # the package
│   ├── __init__.py
│   ├── field.py          # mod-p helpers; sum_{x in F_p} x^r identity; GF(p) wrappers
│   ├── pnomial.py        # pnom(m,s,p), pnom_gt, pnom_le; Pascal + multinomial; memoized
│   ├── reedmuller.py     # monomial enumeration; RM_p(r,m) generator; d_RM; r_max; rtilde
│   ├── puncture.py       # base-p column<->point map (Appendix C); Manhattan set; explicit set
│   ├── triorthogonal.py  # Def-1 checks; shorten/puncture; CSS(G0->X, G'^perp->Z); full-rank guard
│   ├── distance.py       # Delta_p (Thm 4, both branches); generic min-distance; A_d enumeration
│   ├── codes.py          # Code dataclass: provenance, params, gamma, A_d, C
│   ├── distillation.py   # delta_out (eqs 38/39), nbar_T, cost C
│   ├── asymptotics.py    # saddle xi(theta), H_p, gamma_0 two-branch, Table-1 reproduction
│   ├── search.py         # Manhattan w-sweep + randomized explicit-set search; scoring; checkpoint
│   └── cli.py            # `python -m qmsd ...` entry points
└── tests/
    ├── test_pnomial.py        # vs binomials at p=2; sum=p^m; Pascal
    ├── test_reedmuller.py     # d_RM vs Section-5 small-code distance tables
    ├── test_analytic.py       # Thm-5 family; reproduce Table 2 rows (n+k=p^m)
    ├── test_oracle.py         # reconstruct the 10 puncture_locations.json codes; match Table 3
    ├── test_distillation.py   # [[519,106,5]]_5: delta_out~8e-18, C~7.4
    └── test_asymptotics.py    # reproduce Table 1 gamma_0(p), t_0(p)
```

---

## 3. Core data model

```python
@dataclass(frozen=True)
class Code:
    p: int                      # prime qudit dimension
    n: int; k: int; d: int      # [[n,k,d]]_p
    # provenance (one of):
    m: int | None = None              # RM variables
    r: int | None = None              # RM degree used to define the triorthogonal SPACE
    w: int | None = None              # Manhattan cutoff (analytic family) — else None
    puncture_columns: tuple[int,...] | None = None  # explicit 1-indexed set (Appendix C)
    A_d: int | None = None
    full_rank: bool | None = None     # surviving-submatrix rank == claimed; else params invalid

    @property
    def gamma(self) -> float:   # log(n/k)/log(d)
        return math.log(self.n/self.k) / math.log(self.d)
```

Codes carry enough provenance to be **regenerated and re-verified** from scratch.

---

## 4. The two engines

### 4a. Analytic (integer) engine — `pnomial.py` + `reedmuller.py` + `distance.py`
- `pnom(m,s,p)`: convolve the length-`p` all-ones vector `m` times (memoize). Reduces to
  `binom` at `p=2`; `sum_s pnom = p^m`.
- Thm-5 family at `m=3*alpha`, `r=alpha(p-1)-1`: `n=pnom_gt(3a,w,p)`, `k=pnom_le(3a,w,p)`,
  `d=pnom_gt(a,w,p)`; sweep `w`, minimize `gamma`. No matrices — scales to Table 2.
- General `Delta_p(m,rtilde,w)` via Theorem 4 (unified sum; assert the `beta=0`/`beta!=0`
  closed branches agree). **Use `m(p-1)`, not `p(m-1)`** (the paper's typo — see NOTES §4).

### 4b. Explicit (finite-field) engine — `reedmuller.py` + `triorthogonal.py` + `distance.py`
- Build RM generator with `galois.GF(p)`: rows = monomials `x_1^{a1}...x_m^{am}`
  (`a_i<=p-1`, `sum a_i <= r`), columns = points of `F_p^m` (`itertools.product`,
  ordered so column `c` ↔ base-p digits of `c-1`, `x_1` least significant — Appendix C).
- Puncture set `S` → `G'` (drop columns) and `G0` (rows vanishing on `S`, restricted).
  Quantum distance `d = d(G0^perp) = d(PRM_p(rtilde,m;S))`; **full-rank guard** on the
  surviving submatrix (else `[[n,k,d]]` is invalid).
- Minimum distance & `A_d`: enumerate **low-weight** codewords of the punctured dual
  (bounded-weight / information-set style) — the real computational bottleneck (§6).

These overlap on the Manhattan family → cross-check analytic vs explicit `[[n,k,d]]`
for small `(p,m,w)`. Disagreement = bug.

---

## 5. Search engine (`search.py`)

```
for (p, m) in targets:
    r = r_max(m, p);  build RM_p(r, m) once
    # Path A — analytic sweep (cheap, exact, large codes)
    for w in 0 .. m(p-1)-r-1:  record analytic [[n,k,d]], gamma   # Thm 5
    # Path B — randomized explicit search (finds NEW codes)
    repeat (seeded; loop-until-dry or budget):
        S = sample puncture columns (size ~ target k, optionally biased by weight)
        build G', G0; full-rank guard
        compute d (Delta_p if S is a Manhattan level set, else low-weight search), A_d
        gamma = log(n/k)/log(d);  C = n / nbar_T(delta_in)
        keep top-K by gamma AND by C (two leaderboards)
    checkpoint leaderboards to disk
```

Scoring: track **both** `gamma` (asymptotic overhead) and `C` (single-round cost,
eqs 38/39) — the paper notes `C` is the better objective when only one round is used.
Parallelize random seeds (independent), each writing to a shared best-of store.

**Open targets the paper leaves on the table (good first discovery goals):**
- a qutrit (`p=3`) code with `gamma < 1` for `n < 729` (paper found none — explicitly open);
- better `(p=5)` codes than `[[519,106,5]]`, and first small `p=7,11,13` `gamma<1` codes;
- optimize `C` directly (specify `delta_in,delta_out`) rather than `gamma`.

---

## 6. The hard part: minimum distance & A_d

General minimum-distance is NP-hard, and these duals can have large dimension, so brute
force over all codewords is infeasible beyond tiny cases. Mitigations, in order:
1. **Manhattan level sets** → closed-form `Delta_p` (Thm 4); no search needed.
2. **Bounded-weight enumeration:** we only need weights up to a small target `d` (≤ ~6 in
   Table 3) and the count `A_d` at that weight — enumerate low-weight codewords via
   information-set / Brouwer–Zimmermann style rather than full enumeration.
3. **Structure:** the dual is a *punctured Reed-Muller* code; exploit affine symmetry and
   the polynomial picture to bound/representative-search minimum-weight words.
4. **Cap & report:** if a code is too large to certify `d`, record it as a *candidate* with
   an upper-bound `d` and flag it — never silently claim an uncertified distance.

This module is the main engineering risk; build it test-first against the 10 oracle codes.

---

## 7. Milestones (test-driven)

| M | Deliverable | Gate (acceptance test) |
|---|-------------|------------------------|
| M0 | `field`, `pnomial`, RM generator, `d_RM`, `r_max`, `rtilde` | `pnom`→binom at p=2, `sum=p^m`; `d_RM` matches Section-5 distance tables (p=3: 2,2,5,8,8,17,26; p=5: 2,3,4,14; p=7: 2,4,6; p=11: 4,7). |
| M1 | analytic Thm-4/5 family | reproduce Table 2 rows exactly (incl. `n+k=p^m`); analytic vs explicit `[[n,k,d]]` agree on small Manhattan cases. |
| M2 | explicit puncture + CSS + distance + `A_d` | **reconstruct all 10 codes in `puncture_locations.json`; match Table 3 `[[n,k,d]]`, `A_d`, `gamma`.** |
| M3 | distillation metrics | `[[519,106,5]]_5` at `delta_in=1e-3` → `delta_out≈8e-18`, `C≈7.4` (eqs 38/39). |
| M4 | search engine | **rediscover** `[[519,106,5]]_5` (or an equal/better `gamma`) from a randomized search over RM_5(r_max,4). |
| M5 | new-code discovery | beat a published benchmark: qutrit `gamma<1` at `n<729`, or new small `gamma<1` codes at `p=7,11,...`; `C`-optimized variants; asymptotics module steers `(p,m,w)`. |
| M6 | asymptotics | reproduce Table 1 `gamma_0(p)`, `t_0(p)`; `gamma_p(1/6)→2.38309/ln p`. |

M2 is the linchpin: once the oracle codes reconstruct exactly, the search engine in M4–M5
is trustworthy.

---

## 8. Conventions & invariants (enforce as asserts)

- All field arithmetic mod `p`; codes over `GF(p)`.
- Triorthogonality threshold and `r_max` use **`m(p-1)`** (NOT `p(m-1)`).
- `rtilde = m(p-1) - r - 1`; quantum distance evaluated at the **dual** degree `rtilde`.
- Column indices are **1-based** (Appendix C); `c-1 = sum_i x_i p^{i-1}`, `x_1` least significant.
- Every `[[n,k,d]]` requires a passed **full-rank** check or it is flagged invalid.
- `nbar_T = (1-(p-1)/p·delta_in)^n · k` (exponent `n`, not `n/k`); `C = n/nbar_T`.
- Never report an uncertified minimum distance as exact.
```
