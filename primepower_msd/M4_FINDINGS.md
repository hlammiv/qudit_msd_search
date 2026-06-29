# M4 — Low-overhead search: the γ=3 floor and the Z₄ zero-divisor barrier

*Status 2026-06-29.  Multi-strategy workflow (16 agents) + independent re-verification.  Minimal code
frozen in `tests/test_m3_witness.py::test_minimal_8_1_2_code` (43 tests pass).*

## Headline

> **γ < 1 is NOT reachable for the single-ququart transversal-diagonal construction over cyclic Z₄ —
> a hard structural barrier, not an under-searched region.**  The lowest achievable overhead is
> **γ = 3.0**, realized by the minimal code **[[8,1,2]]₄** (and its direct sums [[16,2,2]], [[24,3,2]]).

For context γ = log(n/k)/log(d); the prime baselines are [[15,1,3]]₂ at 2.47 and sublogarithmic
families at γ<1.  Here d (distance) is frozen at 2, so γ ≥ log(n/k)/log 2 ≥ 3.

## The minimal verified code (independently confirmed)

**[[8,1,2]]₄**, transversal **antidiff(S³Z¹)** = diag(ζ₁₆^[0,1,12,13]), inducing logical **diag(1,−i,1,i)**:
- X-stabilizer M_X = ⟨(1,3,3,1,1,1,1,1)⟩;  logical-X gen g = (1,0,0,1,1,2,0,3).
- This is the verified [[9,1,2]]₄ with its free |0⟩ spectator qudit removed (D acts trivially on |0⟩).
- **Dense C^(4⁸) simulation: leakage 1.3×10⁻¹⁵**, pure-diagonal logical, re-certified strict level-3.
- γ = log(8)/log(2) = **3.000**.

**k≥2 codes** ([[16,2,2]], [[24,3,2]], …) are **direct sums** of [[8,1,2]] — they raise k but keep n/k=8
and d=2, so **γ stays 3.0** (the frontier is flat; adding logical ququarts does not lower overhead).

## The barrier (why γ<1 fails) — convergent across 6 independent strategies

Root cause: the **zero divisor 2 ∈ Z₄** (2·2=0).  Two independent ceilings, each alone blocking γ<1:

1. **Distance capped at d = 2.**  The self-condition Φ=0 on M_X is a *ring-triorthogonality* satisfiable
   only at very low M_X rank (≈1); higher-rank self-orthogonal modules are vanishingly rare (0/30k).
   Low-rank M_X ⇒ M_Z = L⊥ is large, and any all-even column carries the **weight-1 Z-logical 2·eₚ**
   (2·eₚ ⊥ every all-even column).  So d=1 generically; the escapes reach d=2; **d=3 was never observed**
   — across all 12 level-3 gates, n ≤ 23, ~470k+ constructions, even after dropping separability and
   transversality and pushing k to 6.
2. **Rate capped at n/k ≥ 8** for a genuinely single-qudit (non-entangling) resource.  Separability ≡
   *disjoint logical support* (overlapping generators induce a logical CCZ — the rejected field route —
   in ~96% of cases, and reintroduce a weight-1 logical otherwise, confirmed by exact V-enumeration at
   n=9,11,13).  So the only separable way to grow k is disjoint direct sums, which freeze n/k and d.
   The minimal block is provably [[8,1,2]] (exhaustive n=6,7).

γ<1 needs n/k < d, i.e. d ≥ 3 at high rate — both ceilings forbid it.  This is a **closed-form catch-22**
(ring algebra), corroborated by 388k+ null searches, not inferred from nulls alone.  Contrast the prime
case: GF(p) has **no zero divisors** and **linear, closed** triorthogonality, which is exactly why
punctured-Reed-Muller over a field reaches high k and high d with γ<1.

## What this establishes

- **Existence (settled, M3):** a genuine single-ququart MSD code exists — minimal **[[8,1,2]]₄**.
- **Low overhead (settled here, negative):** γ<1 is **structurally out of reach** for the single-ququart
  transversal-diagonal construction over cyclic Z₄; the floor is γ=3.0.  This is a publishable no-go for
  *this construction class*, with the mechanism (zero divisor) pinned down.

## Recommended next direction (from the synthesis)

**Attack the exact obstruction — the zero divisor.**  Pursue the prime punctured-Reed-Muller mechanism on
a **GF(4) field** code (no zero divisors, closed linear triorthogonality, where high-rate high-distance
γ<1 codes provably exist) under a **forced single-qudit (separable, non-CCZ) restriction**:

> Does a GF(4) punctured-RM family (k≥2, d≥3, already γ<1) admit a transversal level-3 diagonal gate whose
> induced logical action is **separable** (single-qudit T^⊗ᵏ), not the entangling CCZ the field route
> usually gives?

Either a separable+transversal gate survives ⇒ **γ<1 single-qudit prime-power MSD is reached**, or all
surviving transversal gates are entangling ⇒ a **clean no-go theorem** ("single-qudit MSD with γ<1 forces
a ring with zero divisors, which forces d=2").  Publishable either way; it is the one direction every Z₄
strategy independently pointed to.
