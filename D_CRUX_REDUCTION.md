# The corrected m=2 crux, reduced to a classical arc bound

2026-06-30 — "the real run" at the corrected crux from `D_PROOF_MAP.md` §4. Outcome: **a clean
reduction of the scoped m=2 no-go (p=11, p=13) to a single classical finite-geometry number**, with
the surrounding structure proved and the proof-map's crux corrected — but the final arc bound itself
left open (resists HiGHS/PicoSAT in available time). Companion to `D_NEGATIVE_RESULT.md` §4.1.

## 1. Clean restatement

For p≤13, m=2: let `G0 = shorten(RM_p(r_max,2), S)`, `3r<2(p-1)` (so `G0` is automatically
triorthogonal — the moment hypothesis is inert, `D_PROOF_MAP` L7), `|S|=k`, `n_c = p²−k`, `d` = Z-distance.
The scoped no-go ⟺ **`d ≤ n_c/k` (γ≥1) ⟺ `k(d+1) ≤ p²`.** Tight: p=13, k=26, d=5 ⟹ `156 ≤ 169`, but
`d=6` would need `182 > 169`.

## 2. The reduction (two ingredients) — and the correction to the proof map

> **`d ≤ d_RM(dual) − max_line(S)` (line bound, L6 — PROVED)** with `d_RM` = 8 (p=11) / 9 (p=13),
> and **`max_line(S) ≥ d_RM − n_c/k` (arc bound)** ⟹ `d ≤ n_c/k`, i.e. **γ≥1**.

The arc bound follows from the **(k;s)-arc maxima** `m_s(2,p)` (max size of a set in AG(2,p) with ≤s
points on every line): `k > m_s ⟹ max_line ≥ s+1`. The relevant instances:

| p | k_opt | needed `max_line ≥` | reduces to |
|---|---|---|---|
| 11 | 21 | `⌈8 − 100/21⌉ = 4` | `m_3(2,11) ≤ 20` |
| 13 | 26 | `⌈9 − 143/26⌉ = 4` | `m_3(2,13) ≤ 25` |

**Correction to `D_PROOF_MAP` (which proposed a "rank + 2D full-span codeword" crux).** Computationally
verified that at the *optimal* k the **line bound is already tight** (`true d = d_lines`) and that
search **cannot push `max_line` below 4**; the 2D codeword only binds at *higher* k (k≥22 / 27), where
γ>1 already. So the real obstruction is the **arc bound on `max_line`**, not a 2D codeword.

## 3. What is rigorously established

- **Line bound L6** (workflow, 3× SOUND): `d ≤ d_RM − max_line`, min-weight dual words line-supported
  (the `a=1` GRM regime, `rtilde`∈[10,20)).
- **Two elementary `max_line` lower bounds (proved here):** from `Σ_ℓ|S∩ℓ| = k(p+1)` and
  `Σ_ℓ|S∩ℓ|² = k(p+k)`, the second moment gives `max_line ≥ (p+k)/(p+1)`; the pencil through a point
  gives `max_line ≥ ⌈(k−1)/(p+1)⌉ + 1`. **Both give `max_line ≥ 3` at p=11,k=21** ⟹ `d ≤ 5` — exactly
  **one short** of the needed `d ≤ 4`.
- **Strong counterexample-negative:** simulated annealing, greedy, ILP (HiGHS) and SAT (PicoSAT) all cap
  the ≤3-per-line construction at **20** (p=11) / 22 (p=13), and **no γ<1 crossing** appears in any
  search (best true `d = 4` / `5`, matching `[[100,21,4]]₁₁`, `[[143,26,5]]₁₃`).

## 4. The single remaining gap

**Lemma (open): `m_3(2,11) ≤ 20`** (and `m_3(2,13) ≤ 25`, which has wide margin). Status:
- **Lower bound `≥ 20` is constructed** (HiGHS + greedy both exhibit a 20-set; no 21-set ever found).
- **Upper bound resists computation:** HiGHS could not prove optimality in 230 s; PicoSAT (even with the
  affine-group symmetry breaking — two points fixed by transitivity) did not close the boundary UNSAT in
  1500 s. The over-constrained UNSAT is the hard part.
- **Literature brackets it:** Barlotti `m_n(2,q) ≤ (n−1)q+n` gives `≤ 25`; Ball–Blokhuis–Mazzocca (q odd)
  rules out maximal arcs, so `< 25`. The elementary bounds give only `max_line ≥ 3` (the +1 gap).

**Residual uncertainty (small but honest):** if `m_3(2,11) = 21` (contra all evidence), a `max_line=3`
21-set would have `d_lines = 5`, and one must then check its *true* d — `d = 5` would be a new
`[[100,21,5]]₁₁`, γ=0.97 (a crossing); `d = 4` (a 2D codeword capping it) would keep γ≥1. Every search
indicates `m_3 = 20` and `d ≤ 4`, so this residual is very unlikely, but it is not rigorously closed.

## 5. How to finish

1. **Settle `m_3(2,11) ≤ 20`** — by a literature value (PG(2,11)/AG(2,11) (k;3)-arc tables), a tuned
   exact solver (a leaner at-least-k cardinality encoding; a portfolio/parallel SAT; or hours on lenore),
   or a dedicated combinatorial argument sharpening `max_line ≥ 3` to `≥ 4` for k=21.
2. With that lemma, the **m=2 p=11/13 no-go is a complete proof** (line bound + arc bound). The p=5,7 m=2
   cases are already closed (Singleton-infeasible); m=3 / p=7-m=4 are separate (see `D_P7M4_WINDOW.md`).

## Bottom line

The corrected crux is no longer vague: it is the **single classical inequality `m_3(2,11) ≤ 20`**, with
everything around it proved and the empirics overwhelmingly consistent. This is a clean, well-posed
milestone — short of a complete proof only by one stubborn (but standard) finite-geometry bound.
