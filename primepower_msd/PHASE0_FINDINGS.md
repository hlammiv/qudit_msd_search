# Phase 0 findings — cyclic-clock substrate, M0 calibration, M1 census, K1 existence

*Status 2026-06-28. Code: `primepower_msd/{ring,weyl,clifford_ring,single_qudit_gate,census,hierarchy_search}.py`;
30 regression tests pass (`python -m pytest primepower_msd/tests/test_phase0.py`). Everything is
parametric in d = 2^k and run for **d = 2, 4, 8, 16, 32**.*

## Reproduce

```bash
python -m primepower_msd.census            # M0 calibration + M1 naive-gate census
python -m primepower_msd.hierarchy_search  # K1 existence search (the decisive result)
python -m pytest primepower_msd/tests/test_phase0.py -q
```

## M0 — the level oracle is trustworthy

The Clifford-hierarchy level oracle (`level_of`, built from exact d×d matrices and the operational
tests `is_pauli`/`is_clifford`/`is_level3`) reproduces every textbook qubit (d=2) fact — **S is
Clifford, T = diag(1, e^{iπ/4}) is level 3, H is Clifford** — and the three verified algebraic facts:

- `diag(ζ_{2d}^{x²})` is the **Clifford S** for every d (so on Z_4 the 8th-root quadratic `diag(ζ_8^{x²})`
  is Clifford, **not** magic — the corrected precision fact).
- `x→x²` and `x→x³` are **non-additive** over Z_{2^k} (k≥2) — the anti-collapse lever vs Frobenius.
- the quadratic Gauss sum `Σ_x e^{2πi x²/d} = (1+i)√d` for d ≡ 0 (mod 4).

## M1 — the *naive* single-qudit magic candidates FAIL for d = 2^k (k ≥ 2)

Scanning linear / quadratic / cubic monomial phase gates `diag(ζ_N^{x^j})` at precisions N = 2d, 4d, 8d:

| candidate | d=2 (qubit) | d = 4, 8, 16, 32 |
|---|---|---|
| `diag(ζ_{2d}^{x²})` = S | Clifford | Clifford |
| `diag(ζ_{4d}^{x²})` (lit-map's claimed magic gate) | **level 3** | **level ≥ 4 — NOT in C₃** |
| `diag(ζ_{2d}^{x³})` cubic | Clifford | level ≥ 4 |
| `diag(ζ_{4d}^{x³})` cubic | level 3 | level ≥ 4 |
| BRG `T_s = diag(e^{2πi x/4d})` | level 3 | level ≥ 4 |

> **Correction to LITERATURE_MAP.md.** It asserted the d=4 single-ququart magic gate is
> `diag(ζ_16^{x²})`. The oracle **refutes** this: that gate is level ≥ 4, not 3. The cause is the
> **cyclic wraparound defect** — the finite difference of `x²` is not consistent across the
> x=d−1 → 0 boundary, which exists for d=2^k (k≥2) but not for d=2. BRG's `T_s` gives *universality*
> but is likewise not a Clifford-hierarchy level-3 gate, so it is **not** the right MSD target either:
> universality ≠ level-3. **MSD-by-triorthogonality needs a C₃ gate.**

## K1 — genuine single-qudit level-3 magic gates DO exist (kill-criterion cleared)

A diagonal U whose cyclic phase antidifference satisfies `U X U† = X·diag(v)` lies in C₃ **iff**
`diag(v)` is Clifford. Integrating (antidifferencing) every diagonal Clifford `S^a Z^b` and testing
the result gives an exact, exhaustive existence search:

| d | strict level-3 single-qudit gates (= C₃\C₂) | count | anti-collapse |
|---|---|---|---|
| 2 | qubit T family | 2 | **fails** (shift order 2 → ordinary qubit T, not new) |
| 4 | wraparound-corrected quadratics | **12** | passes → genuine single-ququart magic |
| 8 | | **56** | passes |
| 16 | | **240** | passes |
| 32 | | **992** | passes |

The count is exactly **d(d−1)**. For d = 2^k (k≥2) every one passes the anti-collapse certificate
(the cyclic Pauli group has an order-d shift, so it is not isomorphic to any tensor product of
smaller-qudit Pauli groups) — i.e. these are genuine single-ring-qudit resources, **not** "qubit MSD
in GF(4) disguise". At d=2 the same construction returns the ordinary qubit T and is correctly **not**
flagged as a new resource.

### The actual target gate (d=4)

The single-ququart magic gate is the **wraparound-corrected quadratic**

```
T_4 = diag(ζ_16^{0}, ζ_16^{1}, ζ_16^{4}, ζ_16^{13}) = diag(1, ζ_16, i, ζ_16^{13})
```

i.e. `diag(ζ_16^{x²})` with a **+4 correction at x=3** so the phase closes around the cyclic clock
(naive exps (0,1,4,9) → corrected (0,1,4,13)). This finite, enumerable family — `antidiff(S^a Z^b)` —
is exactly what the ring code must implement **transversally** in M2.

## Consequences for the program

1. **K1 is cleared for d = 4, 8, 16, 32.** Single-qudit prime-power magic is not vacuous — the
   target resource exists and is explicitly enumerated (d(d−1) gates per dimension).
2. **The target gate is sharper than the literature.** It is the wraparound-corrected quadratic
   (antidifference of a diagonal Clifford), not the naive `diag(ζ_{4d}^{x²})` and not BRG's `T_s`.
   The M2 transversality search should target *this* family.
3. **Next: M1.5 + M2.** Build Howell-normal-form ring linear algebra (`galois` cannot represent
   Z_{2^k}), then search small cyclic-Z_4 stabilizer/CSS codes for one on which `T_4^{⊗n}` (or a
   member of the target family) acts transversally and induces a logical T_4. This is the real
   existence question for a *distillation code*; K1 only established the *gate* exists.
4. **Open guard.** The anti-collapse certificate here is the substrate-level Pauli-order argument;
   the full gate-level anti-collapse (no single-qudit Clifford basis change factorizes the *logical*
   gate into a 2-qubit gate) must still be checked on any code found in M2.
