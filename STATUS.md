# qmsd — build status

**Status: complete and green.** All modules implemented; full test suite **243 passed,
0 failed, 2 skipped** (`python -m pytest -q` from the project root; the 2 skips are the
slow `[[667,62,4]]₃` 43.0M-message enumerations, gated behind `QMSD_RUN_SLOW=1` — both pass
when run). Built with adversarial
verification at every module, and the meet-in-the-middle minimum-distance routine was
independently adversarially verified (2,937 fuzz cases vs. null-space ground truth, plus
directional bug-hunts and edge/guard checks — no counterexample).

## What is verified (trustworthy)

The **deterministic math core** reproduces the paper's published results:

| Target (from arXiv:2510.10852) | Status |
|---|---|
| p-nomial identities (binomial at p=2, Σ=pᵐ, Pascal, multinomial) | ✅ exact |
| `d_RM` (Schwartz–Zippel) vs the 16 single-puncture codes | ✅ exact |
| `delta_p` (Theorem 4 distance), both β branches | ✅ exact |
| **Table 2** — 9/9 smallest sublogarithmic codes via the analytic engine | ✅ exact (incl. 18-digit big-ints) |
| **Table 1** — 9/9 asymptotic `gamma_0(p)`, `t_0(p)` | ✅ within 2e-3 |
| **All 10 oracle codes** rebuilt from puncture columns → `n, k, full_rank` | ✅ exact (10/10) |
| **Oracle distance `d`** (meet-in-the-middle, certified) | ✅ **10/10** |
| Oracle **A_d** | ✅ 3/3 attempted (small codes); large ones raise rather than guess |
| Distillation example `[[519,106,5]]₅`: `δ_out≈8e-18`, `C≈7.4` | ✅ reproduced |
| CLI `search` / `reconstruct` / `asymptotic` | ✅ smoke-tested |

**The minimum-distance upgrade.** `qmsd/mindist.py` implements an exact, certified
meet-in-the-middle column-dependency search (syndrome collision on the small parity check,
int64-encoded). It certifies **all 10 oracle distances**, including the 5 the original naive
scan could not — notably `[[519,106,5]]₅` (d=5) in ~10 s and `[[215,28,5]]₃` in 0.7 s (the
authors' own search could not finish the latter in 200 s). Independently adversarially
verified: never under-reports (false collision) nor over-reports (missed codeword); guards
on `d_max`, the distance>6 cap, and the int64 overflow boundary all correct.

Key conventions enforced: triorthogonality / `r_max` use **`m(p-1)`** (not the paper's
misprinted `p(m-1)`); puncture columns are 1-indexed with `x_1` least significant; a distance
is never reported exact unless certified (`d_certified`).

**Search capabilities (parallel + structure-aware).** The explicit `random_search` runs across
processes (`n_jobs`, joblib — ~linear speedup) and supports a `sampler` mode: `uniform`
(default), `capset` (no-3-collinear point sets), and `capset_climb` (cap seed + cap-preserving
distance climb). The cap-set climb reaches the rare high-distance puncture sets where uniform
random stalls — it reconstructs the paper's `[[72,9,3]]₃` in tens of evaluations vs uniform's
0 in 4000. The cap structure was identified by a multi-agent investigation
(`SAMPLING_INVESTIGATION.md`): RM minimum-weight codewords lie on affine lines, so cap sets
(no 3 collinear) avoid the short dual codewords that collapse the distance.

## The exact MacWilliams engine (small-dual regime) — distance > 6 and exact A_d

`qmsd/weightdist.py` is an exact, certified weight-distribution engine for the **small-dual /
high-puncture** regime, where the shortened generator `G0` (= `X_stab`) has small dimension
(`dim(G0) = G0.shape[0]` is tiny when the puncture count `k` is large). It enumerates the
`q**dim(G0)` codewords of `G0` exactly (chunked numpy, big-int histogram) and applies the q-ary
**MacWilliams identity** to get the FULL weight distribution `B` of `G0^perp` — hence the
certified minimum distance (`= min w>0 with B_w>0`) and `A_d = B_d`. All arithmetic is exact
python int (Krawtchouk values, the transform, the divisibility/`B_0=1`/`sum B = q^(n-dim)`
invariants are asserted). It is wired into `codes.py` via `code_certify` (and an opt-in
`exact_budget` on `code_from_puncture`, default 0 = OFF so the search and all existing tests are
unchanged). This lifts two prior limitations:

- **Distance > 6 is now certifiable.** Verified on an MDS witness — a Reed-Solomon `[10,6]`
  generator over `F_11` whose dual is `[10,4,7]`: the engine certifies distance **7** (with
  `A_7 = 1200`), where `min_dependent_columns` provably raises (its d≤6 cap).
- **Exact `A_d` for large codes** where `A_d_logical_Z` refuses. Verified on `[[667,62,4]]₃`
  (`C(667,4) ≈ 8.2e9` blows the subset-scan budget, so the reference raises): the engine returns
  the paper's `A_d = 3972` exactly by enumerating `3**16 = 43.0M` messages.
- All six small-dual oracle codes (`[[20,5,2]]₅`, `[[72,9,3]]₃`, `[[112,13,3]]₅`, `[[200,43,3]]₃`,
  `[[206,37,4]]₃`, `[[667,62,4]]₃`) reproduce the published `A_d` exactly via `B_d`; where cheap
  this is cross-checked against the slow logical reference `A_d_logical_Z` (the minimum-weight
  dual codewords are all logical — no weight-d stabilizers — so `A_d_logical == B_d`).

## Remaining limitations (honest)

1. **Distance certification is capped at d ≤ 6 *only for the meet-in-the-middle path*** (the MITM
   splits into halves of size ≤ 3). Codes with distance > 6 whose dual is **not** small-dual
   (so the exact MacWilliams engine cannot enumerate `G0`) are still left uncertified rather than
   mis-reported. Lifting the MITM cap itself needs size-4 halves (straightforward extension).
2. **Exact `A_d` outside the small-dual regime.** For codes where `dim(G0)` is too large to
   enumerate (`q**dim(G0)` over budget) AND `C(n,d)` is over the `A_d_logical_Z` budget, exact
   `A_d` still raises rather than guess. The exact engine closes this gap in the small-dual /
   high-puncture regime (e.g. all six small-dual oracle codes, incl. `[[667,62,4]]₃`).
3. The explicit search builds `GF(p)` matrices, so it is bounded by `p^m`
   (`search.EXPLICIT_MAX_BLOCK = 750`, comfortably covering the paper's regime). The analytic
   Manhattan engine has no size limit.

## Recommended next steps

- **Extend `A_d` counting** with the MITM idea (the natural next deliverable now that distance
  is fast) — this completes single-round distillation scoring for the large codes.
- Lift the distance cap to d ≤ 8 (size-4 halves) if targeting higher-distance codes.
- Then pursue the paper's open targets — now reachable with the cap-set sampler + parallel
  search: a qutrit `gamma<1` code at `n<729` (the authors found none and were compute-limited),
  and cost-`C`-optimized (vs `gamma`-optimized) search.

## Layout

`qmsd/` (package, incl. `mindist.py` and the validated `data/puncture_locations.json`),
`tests/` (165 tests + `ground_truth.py` + oracle loader), `IMPLEMENTATION_BLUEPRINT.md`
(design), `qmsd/README.md` (usage). The arXiv:2510.10852 paper, verified notes, and typeset
tutorial are maintained locally (not in this repo).
