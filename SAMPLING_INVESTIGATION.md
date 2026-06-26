# Smarter-Sampling Investigation for Low-Overhead Qudit Triorthogonal Codes

## Problem

We search for low-overhead qudit triorthogonal codes by puncturing a Reed-Muller
generator at a set of `target_k` columns. Each candidate is scored by its certified
minimum distance `d` and the overhead yield `gamma = log(n/k)/log(d)` (lower is better).

The default explorer, **uniform-random puncturing**, is weak: on the primary benchmark
it never escapes `d = 2`. This investigation evaluates five "smarter" sampling strategies
against a fixed benchmark to find one that reliably reaches the higher-distance puncture
sets that are known to exist (e.g. the paper's `[[72,9,3]]_3`).

## Benchmark (identical for every strategy)

- **PRIMARY:** `p=3, m=4, target_k=9`. Budget: <= 2000 `code_from_puncture` evals OR 90 s.
  **SUCCESS** = a full-rank candidate with certified `d >= 3`.
- **SECONDARY:** `p=5, m=3, target_k=13` (paper `[[112,13,3]]_5`).
- Prebuilt `r = r_max(m,p)`, `G = rm_generator(r,m,p)`; `max_distance=6`. RNG seeded.

### Uniform baseline (the bar to beat)

**Uniform random: best `d = 2`, 0 successes at `d >= 3` in 4000 trials.**
Confirmed not to be a harness artifact: the paper's oracle puncture set
`(12,29,34,36,53,57,63,67,75)` certifies `d=3` / `[[72,9,3]]_3` / `gamma=1.893`
through this exact API, and 1000 uniform full-rank draws give `d in {1,2}` only
(338x d=1, 657x d=2). The `d>=3` sets are genuinely rare needles (<< 1/4000).

## Results

| Strategy | Primary best d | Success | Evals-to-first-success | Best gamma | Beats uniform | Verified? |
|---|---|---|---|---|---|---|
| **affine_structured** (cap-set) | **3** | **yes** | **55** | **1.8928** | **yes** | yes (d=3, full-rank) |
| hill_climb (+A_d gradient) | 3 | yes | 235 | 1.8928 | yes | yes (d=3, full-rank) |
| simulated_annealing (+A_d) | 3 | yes | 481 | 1.8928 | yes | yes (d=3, full-rank) |
| seed_perturb_manhattan | 3 | yes | 594 | 1.8928 | yes | yes (d=3, full-rank) |
| genetic (+A_d, Manhattan prior) | 3 | yes | 1865 | 1.8928 | yes | yes (d=3, full-rank) |
| manhattan_biased | 2 | no | — (never) | 3.0 | **no** | yes (d=2, full-rank) |
| *uniform (baseline)* | *2* | *no* | *— (0/4000)* | *3.0* | *—* | — |

All five non-trivial strategies were independently reproduced (re-run, distance and
full-rank re-certified through `code_from_puncture`); no discrepancies were found.

### Which strategies beat uniform

Five of six beat the baseline decisively: **affine_structured, hill_climb,
simulated_annealing, seed_perturb_manhattan, genetic** all reach `d=3` and reconstruct a
`[[72,9,3]]_3` code. Only **manhattan_biased FAILS** — an i.i.d. Manhattan-weight marginal
bias (swept over beta in `[-3, 2]`, deterministic low/high-weight, and weight-window
pools) never exceeds `d=2`, because the paper's oracle puncture points have *mixed/high*
Manhattan weights `(2,3,3,4,4,5,6,6,6)`, not a low-weight shell. `d>=3` is essentially
uncorrelated with the marginal weight of punctured points, so no i.i.d. weight-biased
sampler can concentrate on the valid sets.

## Winner: `affine_structured` (cap-set-biased sampling + cap-preserving hill climb)

File: `/tmp/claude-1000/-home-hlamm-Desktop-QC-prime-msd/2ffd6aea-a00e-403c-8252-73440ed02d50/scratchpad/sampling/affine_structured.py`

**Why it wins.** It is the only strategy grounded in *why* `d` is low, not just a generic
metaheuristic. RM minimum-weight codewords live on affine subspaces, so **collinear**
puncture points create short dual codewords that kill the distance. The key empirical
discovery is that the paper's own `[[72,9,3]]_3` puncture points form a **9-cap** (a set
in `F_3^4` with no 3 collinear points — 0 of 84 triples sum to 0 mod 3). The strategy
biases sampling toward cap sets via `column_to_point`, then runs a cap-preserving,
full-rank-preserving swap climb with restarts.

Concrete wins over the field:
- **Fastest first hit:** `d>=3` at **55 evals / 0.75 s** (primary, seed 7), vs 235 / 481 /
  594 / 1865 for the others.
- **Best structural insight / generalization:** the cap property transfers to the
  secondary benchmark essentially for free — `p=5,m=3` succeeds **5/5 seeds, on the first
  cap seed (1-3 evals, ~0.1 s)**, producing `[[112,13,3]]_5`, gamma=1.960.
- Reconstructs the paper code exactly: primary best `[[72,9,3]]_3`, gamma=**1.8928**,
  cols `[16,19,29,31,32,45,50,67,77]`.

**Runner-up nuance.** If the single most important axis were *single-run reliability*
rather than speed/insight, `hill_climb` (7/8 seeds) and `simulated_annealing` (9/10 seeds)
are more robust than `affine_structured` (6/10 seeds) — all three share the same enabling
trick (an `A_d` minimum-weight-multiplicity gradient; see caveats). The cap-set winner is
chosen for its fastest convergence, its explanatory power, and its near-instant secondary
performance, with the understanding that its reliability gap is easily closed by more
restarts.

## Comparison vs the paper (Table 3)

Paper best at `p=3, m=4`: `[[72,9,3]]`, gamma = 1.89.

Every successful strategy **matches** this entry: each finds a *full-rank* `[[72,9,3]]_3`
code with gamma = 1.8928 (the same overhead). The winner finds it in 55 evals via a
*distinct* puncture set from the oracle. **No strategy strictly beats the paper** — none
found a code with gamma < 1.89 or with `d > 3` at the target `k`. This is expected and
consistent: `d` certification is capped at `max_distance=6`, the paper itself reports
`d=3` for both targets, and at fixed `(n,k)` gamma is a pure decreasing function of `d`, so
beating the paper would require `d>=4` puncture sets, which do not appear in the explored
neighborhoods. Secondary likewise matches the paper's `[[112,13,3]]_5`, gamma=1.960.

## Honest caveats

1. **The `d>=3` points are isolated needles.** All 648 one-swap neighbors of the oracle
   `[[72,9,3]]` set collapse to `d=2`. The integer certified distance therefore has **zero
   local gradient** — this is exactly why uniform random and naive distance-climbing fail.
2. **`A_d` is the load-bearing trick for three strategies.** `hill_climb`,
   `simulated_annealing`, and `genetic` only work after switching fitness to use `A_d`
   (the multiplicity of minimum-weight codewords, `compute_A_d=True`) as a smooth
   surrogate: among `d=2` sets, low `A_2` basins (A_2 ~ 10) sit adjacent to the rare `d=3`
   spikes. Without it these degenerate to random search. `compute_A_d=True` costs ~50%
   more per eval (~0.015 s vs ~0.010 s) and is a deliberate deviation from the suggested
   `compute_A_d=False`. The winner, `affine_structured`, does **not** need `A_d` — the cap
   constraint supplies the structural gradient instead.
3. **Stochastic reliability (per single run, primary, 2000-eval budget):**
   `simulated_annealing` ~9/10, `hill_climb` ~7/8, `affine_structured` ~6/10,
   `seed_perturb_manhattan` succeeds on all seeds tried, `genetic` ~1/5 (weakest;
   first hit at 1865 evals leaves almost no margin). A restart with a fresh seed recovers
   failures in every case.
4. **Budget binding:** for primary the **eval count** binds; for secondary (`p=5`, ~34-40
   ms/eval) **wall-clock** is the tighter constraint for the metaheuristics (SA hit 80.7 s
   of the 90 s cap). The winner's secondary is budget-trivial (first cap seed).
5. **`max_distance=6` cap** means any `d>6` would be invisible; observed best `d` stayed at
   3 throughout, matching the paper.

## Recommendation: integrate the winner into `qmsd.search.random_search`

Do **not** modify `qmsd`. The proposed change (for a maintainer) is additive and
backward-compatible. Current signature:

```python
random_search(p, m, trials, seed=0, target_k=None, max_distance=6, n_jobs=1) -> list
```

Add an optional `sampler` argument that selects how each candidate puncture set is drawn,
defaulting to today's behavior so nothing changes for existing callers:

```python
random_search(p, m, trials, seed=0, target_k=None, max_distance=6, n_jobs=1,
              sampler="uniform")   # NEW: "uniform" | "capset" | "capset_climb"
```

- `sampler="uniform"` (default): unchanged — preserves reproducibility and the documented
  baseline.
- `sampler="capset"`: replace the per-trial uniform `target_k`-column draw with a greedy
  cap-set draw — map columns to points with `column_to_point`, accept a candidate column
  only if it forms no collinear triple (`a+b+c == 0 mod p` for `p=3`; general parallel
  test `b-a parallel to c-a` for `p>3`) with the already-chosen points, restart the greedy
  build on a dead end. This is a drop-in replacement for the sampling step inside the
  existing trial loop; the dedup/certify/sort/keep machinery is untouched.
- `sampler="capset_climb"`: after each cap seed is certified, run the cap-preserving,
  full-rank-preserving swap climb with restarts (accept non-decreasing `d`) until that
  trial's local eval sub-budget is exhausted, then record the best. This is the full winner
  and is what delivers the 55-eval first hits.

Lift the cap-set seed builder, the collinearity predicate, and the swap-climb loop directly
from `affine_structured.py` (referenced above). Keep `n_jobs` parallelism intact (cap seeds
are independent per worker, derive each worker RNG from `seed` as today). Suggested doc note:
`"capset"`/`"capset_climb"` are structure-aware modes that dramatically raise the hit rate
for `d>=3` puncture sets where uniform sampling stalls at `d=2`; they change which sets are
sampled, so results are reproducible for a fixed `(seed, sampler, n_jobs, trials)` but are
not comparable to the uniform-mode stream. For robustness, expose the restart count and
recommend re-running with a fresh seed on a failed single run (the documented ~6/10 ->
near-certain over a few restarts).
