# qmsd

Discovery of **qudit triorthogonal codes for magic state distillation** with low
overhead — a reimplementation and extension of the search in
[arXiv:2510.10852](https://arxiv.org/abs/2510.10852) (Saha & Prakash, 2025).

Given a prime qudit dimension `p`, `qmsd` builds triorthogonal codes by puncturing
maximal-degree Reed-Muller codes over `F_p` and scores them by the yield parameter
`gamma = log(n/k)/log(d)` (overhead exponent; `gamma < 1` is *sublogarithmic*) and by
single-round distillation cost `C`.

Authoritative math reference: `arXiv:2510.10852`.
Software design: `../IMPLEMENTATION_BLUEPRINT.md`.

## Install / run

Requires Python ≥ 3.10 with `numpy`, `galois`, `scipy`. Run from the project root
(the package imports locally — no install needed):

```bash
# scan a range of m for low-gamma codes at p = 5
python -m qmsd search --p 5 --m 4 --trials 5000

# rebuild a published code from its puncture columns and check it
python -m qmsd reconstruct --label "[[20,5,2]]_5"

# asymptotic optimal yield gamma_0(p)
python -m qmsd asymptotic --p 5
```

To retain search results and make them visible in the Streamlit explorer:

```bash
python -m qmsd search --p 5 --m 4 --trials 5000 --output runs/p5-m4.json
python -m qmsd catalog import runs/p5-m4.json
streamlit run app.py
```

The import is idempotent and preserves distance/full-rank certification metadata.
Imported records live in `qmsd/data/catalog/` and are discovered automatically.

```python
from qmsd.search import search
res = search(5)                       # given p, redo the paper's search
for c in res["best_by_gamma"][:5]:
    print(c.label, round(c.gamma, 3))
```

## Two engines

1. **Analytic (integer) engine** — `codes.code_from_manhattan`, `search.manhattan_sweep`.
   The Manhattan-weight family (Theorems 4–5): `n = [m,>w]_p`, `k = [m,<=w]_p`,
   `d = Delta_p(m, rtilde, w)`, all exact integer arithmetic with **no matrices**, so it
   scales to the astronomically large Table-2 codes (e.g. `[[2.9e17, …, 21700]]_2`).

2. **Explicit (finite-field) engine** — `codes.code_from_puncture`, `search.random_search`.
   Builds the actual `RM_p(r_max,m)` generator over `GF(p)` (via `galois`), punctures an
   arbitrary column set, and computes `[[n,k,d]]`, `A_d`, `gamma` directly. This is how the
   paper's best small codes were found; it is bounded by `p^m` (matrix size) and by the cost
   of certifying minimum distance / `A_d`.

## Sampling & parallelism

`random_search` (and `search`) take:

- **`sampler`** — how each candidate puncture set is drawn:
  - `"uniform"` (default): i.i.d. random sets — stalls at low distance on hard cases.
  - `"capset"`: **cap sets** (no 3 collinear points). RM minimum-weight codewords lie on
    affine lines, so collinear punctures collapse the distance; caps avoid that obstruction.
  - `"capset_climb"`: a cap seed + a cap-preserving distance hill-climb — reaches the rare
    high-distance puncture sets where uniform random fails (it reconstructs the paper's
    `[[72,9,3]]_3` in tens of evaluations vs uniform's 0 in 4000). See
    [`../SAMPLING_INVESTIGATION.md`](../SAMPLING_INVESTIGATION.md). Per-run reliability is
    stochastic — use more `trials` / fresh seeds (or `n_jobs`) for robustness.
- **`n_jobs`** — process-level parallelism (trials are independent): `1` (serial, default),
  or `>1`/`-1` for joblib worker processes. Reproducible for a fixed
  `(seed, sampler, n_jobs, trials)`.

```python
from qmsd.search import random_search
codes = random_search(3, 4, trials=200, target_k=9, sampler="capset_climb", n_jobs=-1)
```

## Module map

| module | contents |
|--------|----------|
| `field` | `F_p` arithmetic helpers, power-sum identity |
| `pnomial` | p-nomial coefficients `[m,s]_p`, tail/head sums, multinomial cross-check |
| `reedmuller` | `r_max`, `r_tilde`, `d_rm` (Schwartz–Zippel), monomials, points, `rm_generator` |
| `puncture` | base-p column↔point map (Appendix C), Manhattan puncture set |
| `triorthogonal` | triorthogonality check, puncture/shorten/dual, `build_triorthogonal_code` |
| `distance` | analytic `delta_p` (Thm 4), exact `min_distance`, `A_d_logical_Z` |
| `codes` | `Code` dataclass, `code_from_puncture`, `code_from_manhattan` |
| `distillation` | `nbar_T`, `cost`, `delta_out_*` (eqs 38/39) |
| `asymptotics` | saddle `xi`, `H_p`, `gamma0`, `optimal_gamma` (Table 1) |
| `search` | `manhattan_sweep`, `random_search` (samplers + `n_jobs`), top-level `search(p)` |
| `sampling` | cap-set geometry (no 3 collinear points) for structure-aware puncturing |
| `cli` | `python -m qmsd search/reconstruct/asymptotic` |
| `oracle` | loads the 10 validated paper codes (the correctness oracle) |

## Correctness oracle

The package is anchored to the paper's published results. `tests/` reproduces:

- **Table 2** (9/9 sublogarithmic codes) via the analytic engine — exact, including big-ints.
- **Table 1** (9/9 asymptotic `gamma_0(p)`, `t_0(p)`).
- **All 10 search codes** rebuilt from their puncture columns to the correct `[[n,k,d]]` and
  `full_rank` — distances **certified for all 10** via the meet-in-the-middle routine
  (`mindist.py`), including `[[519,106,5]]_5`. `A_d` is certified for the small codes; exact
  large-`A_d` counting raises rather than guesses.
- The single-round distillation example (`[[519,106,5]]_5`: `delta_out ≈ 8e-18`, `C ≈ 7.4`).

Run the suite from the project root: `python -m pytest -q`.

### Key conventions (enforced)

- Triorthogonality / `r_max` use **`m(p-1)`**, not the paper's misprinted `p(m-1)`:
  `r_max = floor((m(p-1)-1)/3)`, `r_tilde = m(p-1) - r - 1`.
- Puncture columns are **1-indexed**; column `c` ↔ point `x` via `c-1 = x_1 + x_2 p + …`
  (`x_1` least significant).
- A minimum distance is **never reported as exact unless certified**; uncertified codes carry
  `d_certified = False`.
