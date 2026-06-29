# M1.5 + M2 findings — ring linear algebra, and the search for a single-ququart distillation code

*Status 2026-06-29. Code: `ringlinalg.py` (M1.5), `ring_css.py`, `ring_transversal.py`, `m2_search.py`;
37 regression tests pass (`python -m pytest primepower_msd/tests/`).*

## Reproduce

```bash
python -m primepower_msd.m2_search           # the search (verdict at the end)
python -m pytest primepower_msd/tests/ -q
```

## M1.5 — Howell-normal-form ring linear algebra (`ringlinalg.py`)

`galois` represents the *field* GF(2^k), not the *ring* Z_{2^k}.  `ringlinalg.py` supplies the
canonical replacement — the **Howell normal form** over Z_n — plus span membership, module size,
and the (vectorized) right kernel / dual module.  **Validated against brute-force enumeration**:
1000+ random matrices over Z_4/Z_8/Z_16, exact agreement on module equality, size, and membership
(`test_m1_5_m2.py`).  This is the linear-algebra substrate for all ring-code work (M2 here, distance
in M3).

## M2 — the key reduction (why this is tractable)

Because the magic gate is **diagonal**, transversality has a sharp finite form.  With the X-basis
logical states |psi_c> ∝ sum_{a in M_X} |c+a>, the transversal gate U = D^{⊗n} (D = diag(zeta_N^phi))
never mixes M_X-cosets, so

> **U preserves the codespace  ⟺  the total phase Phi(y) = sum_i phi(y_i) is constant (mod N) on
> every M_X-coset**,  and then the induced logical gate is diag over the d cosets.

No d^n state vectors are needed — M2 is finite phase-arithmetic over cosets.  The shortcut is
**cross-checked by direct simulation** in C^{d^n} (`brute_force_transversal_check`): on the witness
below the codespace leakage under U^{⊗n} is exactly 0.

## M2 — results

### (a) Transversal level-3 single-ququart gates DO exist on Z_4 CSS codes

The **non-free loophole** the literature flagged (Gluesing-Luerssen–Pllaha: free ring codes are
distance-bounded, non-free codes are the open case) is decisive: requiring *free* modules gives
**nothing** (no k=1 cyclic free weakly-self-dual code even exists for n ≤ 5).  Allowing **non-free**
codes immediately yields transversal hits, e.g.

```
[[3,1]]_4 CSS,  X-stabilizers (0,0,2),(0,2,0),  uniform transversal antidiff(S^2) = diag(1, ζ8^-1, 1, ζ8^-1)
  -> induced logical gate = a strict level-3 single-ququart gate;  brute-force codespace leakage = 0.
```

So a transversal, single-system, Clifford-hierarchy-level-3 gate on a *cyclic* Z_4 CSS code is real
(not a GF(4)/two-qubit artifact — it passes the anti-collapse certificate).

### (b) ...but only on TRIVIAL (distance-1) codes — the distillation question is a NULL

The witness above has the X-stabilizers acting only on qudits 2–3, so qudit 1 is the **unencoded**
logical qudit: **distance 1**, no error protection, cannot distill.  Adding a distance filter:

| n | k=1 cyclic codes searched | with distance ≥ 2 | transversal level-3 hits |
|---|---|---|---|
| 3 | 4 (exhaustive, all ranks) | 0 | distance-1 only |
| 5 | 89 (exhaustive small-rank + sampled) | 0 | distance-1 only |
| 7 | 350 (exhaustive small-rank + sampled) | **346** | **none on the 346 protected codes** |

**Verdict.** Across **443** weakly-self-dual k=1 cyclic Z_4 codes × the **complete 12-gate** level-3
single-qudit diagonal family, under **uniform AND {0,1}-addressable** transversal application, the
only transversal level-3 hits (58 of them) are **distance-1**.  **No distance-≥2 single-ququart
distillation code was found.**  There is a genuine tension between *transversal level-3* and *error
protection* in this family — exactly the ring-precision-vs-orthogonality obstruction the literature
predicts.

## What this does and does NOT establish

- **Does:** single-qudit level-3 magic gates exist *and* act transversally on cyclic-Z_4 CSS codes
  (K1 + a transversal-gate existence proof), but in the natural weakly-self-dual family they coincide
  only with unprotected (distance-1) codes.  The non-free structure is essential; free codes give
  nothing.
- **Does NOT:** prove a no-go.  The search is bounded to: weakly-self-dual codes (M_X = M_Z),
  n ≤ 7, uniform/{0,1}-addressable gates, d = 4.  It says the *simplest* family does not yield a
  distillation code and pinpoints *why* (the protected codes reject every transversal level-3 gate).

## M2b — general CSS (M_X != M_Z) + a sharp structural fact (`m2b_search.py`)

Decoupling X- and Z-stabilizers (general CSS) was the most promising broadening, via the gate's
**translation stabilizer** Stab(phi) = {a : Phi(.+a) == Phi everywhere} (a module).  A diagonal gate is
(globally) transversal exactly when the X-stabilizers M_X subseteq Stab(phi).  This exposes the crux:

> **Of the 12 strict-level-3 single-ququart gates, only 2 have a nontrivial Stab(phi)** (stab1 = {0,2}),
> and those 2 are the **qubit-like** gates whose phase depends only on x mod 2 (e.g. diag(1, ζ8⁻¹, 1, ζ8⁻¹)).
> The other **10 — the genuinely-quadratic gates** (e.g. diag(1, ζ16, i, ζ16¹³)) — have **trivial
> stab1 = {0}** and **cannot be globally transversal on any nontrivial code.**

A Stab(phi)-guided general-CSS search (537,600 k=1 cyclic codes) found that the 2 transversable
(qubit-like) gates still meet **only distance-1 codes** — no distance-≥2 witness.

**Honest caveat (important).** "Cannot be transversal" above is the *global* condition.  The true
per-coset transversality (Phi constant on M_X-cosets within M_Z^perp only) is weaker, so a
genuinely-quadratic gate could still be transversal on a special code where the constraint holds on
M_Z^perp — this is precisely the **ring-triorthogonality (R1)** condition (the c-linear part
2·sum c_i a_i == 0 mod N for c in M_Z^perp).  Random search will not find these measure-zero codes;
they require the analytic co-design.  **So the genuinely-quadratic transversal question is open, = R1.**

## Performance — parallel + optimized (`parallel.py`)

The searches are embarrassingly parallel (independent trials → hit set).  Parallelized with joblib/loky
(`parallel.py`): single-threaded BLAS per worker (no oversubscription), conservative default of
`cpu − ⌈cpu/3⌉` workers (20 → 14, override via `QMSD_JOBS`), bounded n so the d^n dual grid stays small.
Algorithmic wins (~21×/core hot loop): a **k=1 size-prune** (`|M_X|·|M_Z| == d^{n-1}`) skips the d^n dual
on codes that cannot be k=1; the coset-phase check is **vectorized**; `dual(A)` is **reused** across
Z-stabilizer samples.  Result: `m2_search` 53 s/14 cores (was timing out serially); `m2b_search`
537,600 codes in 151 s/14 cores (was 216,000 in 352 s/18 cores). 37 tests pass.

## Next (M3 / R1)

1. **Per-coset (ring-triorthogonality, R1) search for the quadratic gates** — the real open question
   M2b leaves: co-design M_X, M_Z so 2·sum_{i} c_i a_i == 0 (mod N) on M_Z^perp.  Analytic, not random.
2. **M3 — Howell-scaled ring distance/A_d** to reach larger / punctured n (prime codes live at
   [[8,1,2]], [[80,1,5]]); the brute d^n dual must become a Howell kernel.
3. **Escalate toward a no-go** for the weakly-self-dual family (distance ≥ 2 + transversal level-3
   impossible) — a publishable negative.
4. **d = 8, 16, 32** — repeat with the (already parametric) substrate; the d=8 Hoggar structure may
   differ.
