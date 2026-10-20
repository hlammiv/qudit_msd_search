# Universality check for ring-qudit (Z₄) magic gates — what is actually distilled

*2026-10-19.  Modules: `a2.py` (sector-escape test), `phase_tables.py` (table hierarchy oracle,
catalog), `rgt.py` (hosting + induced logical).  Tests: `tests/test_a2.py`, `tests/test_phase_tables.py`,
`tests/test_rgt.py`.  Gate catalogs: `runs/diverse/gate_catalog_m{1,2,3}.json`.*

## 1. The criterion

For d = 2^k the Clifford adjoint action on the Pauli algebra is reducible (Borda–Rincón–Galindo,
arXiv:2512.20787): the Paulis X^uZ^v split into sectors W_j by gcd(u,v,2^k) = 2^j.  Define

    A_j = span{ X^a Z^b : all a_i, b_i ∈ 2^j Z_d },   j = 1..k−1.

At d = 4 the single sector algebra **A₂ = span{XᵃZᵇ : a, b all even}** is a *maximal abelian*
subalgebra (X², Z² commute): one fixed product basis {(|0⟩±|2⟩), (|1⟩±|3⟩)}.  Its normalizer N(A₂) is
the group of monomial unitaries in that basis — a closed proper subgroup of U(4ⁿ).  Hence

> **Clifford + G is universal on ququarts only if G escapes A₂ under conjugation.**

Exact test for a diagonal gate D = diag(ζ_{4d}^{P(x)}) (`a2.preserves_sector`):

    D normalizes A_j  ⇔  Δ_{2^j e_i} Δ_{(d/2^j) e_l} P ≡ 0 (mod 4d)  for all i, l;

at d = 4: `q(x+2) − q(x) ∈ {0, 8}` for every x.  Validated against a dense Pauli-expansion test on all
4096 single-ququart tables, random 2-ququart tables, and d = 8 (both sectors).

## 2. Results

| gate | level | escapes A₂ | odd 3rd diff. |
|---|---|---|---|
| Fourier, S = ζ₈^{x²}, CSUM, CZ = i^{xy} (Cliffords) | ≤2 | no | – |
| **ring-CCZ i^{xyz}** (Phase-5 target) | 3 | **no** | yes |
| i^{x³} = diag(1,i,1,−i) (the [[8,1,2]] gate) | 3 | no | no |
| every diag(i^F), F polynomial mod 4 (40 random, 1–3 qudits) | ≤3 | no | – |
| **carrier ζ₁₆^{[0,1,4,13]}** (the [[8,2,2]] gate) | 3 | **yes** | yes |
| ζ₈^{x³} = [0,2,0,6], ζ₁₆^{x} | 4 | yes | – |

Consequences:
* **The γ<1 Phase-5 codes `[[10,4,3]]₄`, `[[12,5,3]]₄` distill a non-universalizing gate.**  Their
  logical ring-CCZ is genuine order-4 level-3 magic, but Clifford + ring-CCZ ⊂ N(A₂).  (Reproduced
  with `rgt`: hosted, k = 4, induced logical level 3, `escapes_A2 = False`, odd third difference True.)
* Of the 48 strict level-3 single-ququart tables (with P(0)=0), **32 escape** and 16 do not; the 12
  `hierarchy_search` carriers are among the 48.
* **Transmission matters.**  On `[[8,2,2]]₄` every one of the 48 tables is hosted, but **every induced
  logical gate is non-escaping and has no odd third difference** (carrier → F = [[0,4,0,4],[12,4,12,4],
  [0,12,0,12],[4,4,4,4]], a mod-4 phase).  For a uniform carrier-type layer (Δ³q = 4):
  `Δ_{2e_i}Δ_{2e_l}F = 8·#{t : a_{i,t}, a_{l,t} both odd}`, so the logical escapes iff some logical-X axis
  has an odd number of odd entries or two axes have an odd odd-overlap (`rgt.axis_parity_report`).
  `[[8,2,2]]`'s axes have odd-weights 4, 4 and overlap 2 ⇒ no transmission.  **So at present no code in
  the program distills a universality-enabling logical gate.**  The Phase-1 "genuine ω=i logical"
  claim for `[[8,2,2]]` tested level-3-ness and the code's logical commutator, not the gate's sector escape.

## 3. Entangling gates: the one-body "lemma" is false, but nothing new hides in them

A first derivation (polynomial phases with translation-consistent formulas) suggested A₂-escape is a
purely one-body phenomenon.  On genuine function tables it is **false**: the exhaustive m = 2 catalog
(131 072 level-≤3 tables with P(0)=0; 130 048 strict) contains **65 536 escapers whose escape is carried
by the entangling part**, e.g. P = [[0,0,0,0],[0,2,0,10],[0,0,8,8],[0,10,8,10]] (ζ₈ phases, level 3
dense-verified, order 8, zero one-body marginals, not single-frame local).

However, with W_m := span{t(l·x) : l ∈ Z₄^m a linear form, t a level-≤3 single-qudit table} + Clifford
tables, **W_m equals the entire level-3 module for m = 2 (2²¹) and m = 3 (2³⁹)**
(`phase_tables.rank1_module`).  Every level-3 diagonal gate on 2–3 ququarts — ring-CS, ring-CCZ, the
entangling escapers — is therefore an exact product of Clifford-frame-rotated single-qudit level-3 gates
and Cliffords, exactly as qubit CS = T_x T_y T†_{x+y} and CCZ = 7 T's.  Distilling entangling gates is an
*efficiency* choice (as for qubit CS/CCZ), not a new resource; universality is certified by sector escape,
which entangling (e.g. across-block CS-like) gates can carry.

## 4. The ring triorthogonality theorem for 0/1 lifts, and the first transversal T-like ring codes

For a binary CSS code lifted by `lifts.lift_css_from_x` (M_X = 0/1 lift of the X-stabilizers, logical
X's = 0/1 lifts, M_Z = (M_X + L)^⊥ over Z₄) and the carrier q = [0,1,4,13] (Δq = [1,3,9,3], Δ₂q = [4,12,12,4],
Δ³q = 4), the hosting conditions `Δ_aΦ ≡ 0 on supp` reduce to binary overlap conditions:

    Q1  |s| ≡ 0 (mod 16)            for X-stabilizers s
    Q2  |s ∧ u| ≡ 0 (mod 8)         for u a stabilizer or a logical representative
    Q3  |s ∧ u ∧ u'| ≡ 0 (mod 4)

Q1–Q3 on stabilizers alone ⇔ **the binary X-stabilizer code is quadruply even (all weights ≡ 0 mod 16)**
(|s∧s'| = (|s|+|s'|−|s+s'|)/2, |s∧s'∧s''| = (Σ|·| − Σ|·+·| + |s+s'+s''|)/4).  With logicals, Q2 says the
logical weight is constant mod 16 on each stabilizer coset; an ODD coset weight transmits escape
(`rgt.axis_parity_report`).  This is the Z₄ analog of Bravyi–Haah triorthogonality (mod 8 → mod 16).
By Ax's theorem RM(r,m) is quadruply even iff 4r < m, so 0/1 lifts of RM(2,m) need m ≥ 9 (n = 512) and
RM(1,m) needs m ≥ 5.

Results (`lifts.punctured_rm_triorthogonal` + `lift_css_from_x`, distances by `z4_distance`):

| code | hosts (of 48) | escaping logical | d | γ |
|---|---|---|---|---|
| [[15,1,3]]₄ (RM(1,4)*) | 16 (non-escaping only; weight 8) | 0 | 3 | – |
| **[[31,1,3]]₄ (RM(1,5)*)** | **48** | **32** (carrier → logical ζ₁₆^{[0,15,12,3]}) | 3 | 3.13 |
| **[[63,1,3]]₄ (RM(1,6)*)** | 48 | 32 | 3 | 3.77 |
| cube [[8,3,2]]₄ | 16 | 0 (even logical weights) | 2 | – |

These are the first ring-qudit codes with a transversal, universality-enabling (sector-escaping) logical
gate; hosting was cross-checked by 30 000 random Δ_aΦ(x) evaluations on [[31,1,3]].  Ring punctured
RM(1,m) with k ≥ 2 punctures hosts and escapes automatically but the affine span of the punctures kills
the distance (k=2 → d=2, k≥3 → d=1); RM(2,m) lifts fail Q2/Q3 for m ≤ 7 (as predicted).  The γ<1 front
is therefore the quadruply-even families with many logicals (ring punctured RM(2,9), GR(4,m) non-CSS
codes) — see `diverse_search.py` and `runs/diverse/`.

### Experiments E3 / E4 (negative, as the theorem predicts)
* **E3** — 400 GR(4,3) non-CSS graph codes (N = 12, k = 6): none hosts any single-qudit level-3 table; only 29
  pass S1 for some table (`runs/diverse/catalog_E3.json`).
* **E4** — the constrained sampler (X-generators drawn from the q-weight-0 set, S2/S3 enforced) stalls at 1–2
  generators at n = 12: random vectors essentially never satisfy the pairwise mod-16 conditions.  Hosting
  families must be CONSTRUCTED (quadruply-even binary codes), exactly as triply-even codes for qubit T.

## 5. Loophole probe (E6): entangling-only universality via non-semi-Clifford gates

Pauli-controlled Cliffords Σ_x |x⟩⟨x| ⊗ C_x are level-3 only if C_{x+1}C_x† is a Pauli for all x, i.e.
C_x = P_x C_0 with the controlled-Pauli Σ|x⟩⟨x|⊗P_{x+1}P_x† Clifford (linear control); such gates are
Clifford-conjugates of diagonal 2-ququart tables, hence covered by §3 (W_2 = level-3 module).  A genuinely
non-semi-Clifford level-3 ququart gate (Gottesman–Mochon type needs ≥3 qudits for qubits) remains the only
unexplored route; for 1–2 ququarts the diagonal catalog is complete.

## 6. Status of the "T-like" claim

What the program has: existence of genuine level-3 ring gates hosted by codes (`[[8,2,2]]₄`), and γ<1
codes hosting a (non-transversal, non-escaping) ring-CCZ.  What it does not have: any code whose hosted
transversal gate induces an **escaping** logical gate.  The search objective in `diverse_search.py` is
therefore: hosting ∧ transversal/across-block ∧ `induced_logical` escapes ∧ d ≥ 2 (then maximize k/n).

Reproduce: `python -m pytest primepower_msd/tests/test_a2.py primepower_msd/tests/test_phase_tables.py -q`;
`python -c "from primepower_msd.phase_tables import catalog; print(catalog(2))"`.


## 7. Two-tier ("CCZ most of the time, T for universality") — the exponent cannot be shared

**Theorem (finite cheap group).**  In the A₂ eigenbasis {(|0⟩±|2⟩)/√2, (|1⟩±|3⟩)/√2}^{⊗n} every Clifford
generator (F, S, CSUM) and the ring-CCZ i^{xyz} is a MONOMIAL matrix whose phases are 8th roots of unity
(checked exactly: phase exponents in {0,2,4,8}/16).  Hence H = ⟨Clifford, ring-CCZ⟩ ⊂ μ₈ ≀ S_{4ⁿ} is FINITE
(the large element orders observed — up to 9840 — are permutation orders, not evidence of infinitude).  The
same holds for any A₂-preserving gate with root-of-unity phases.

**Consequence.**  With t universality-enabling (escaping) gates and arbitrary H-padding at most |H|^{t+1}
unitaries are reachable, so ε-approximating a generic unitary needs t = Θ(log 1/ε) escaping gates.  The
asymptotic overhead exponent of any protocol is therefore γ of the escaping (T-like) resource; cheap ring-CCZ
factories improve constants (Toffoli-heavy algorithms) but not the exponent.  (A cheap gate with irrational
phases would break this, but such gates are not Clifford-hierarchy gates and are not distillable by MSD.)

## 8. Fault-tolerant ring-CCZ codes (step 2, `ccz_search.py`)

Across-block ring-CCZ on three copies of a lifted CSS code hosts iff
    (C4)  |a ∧ u ∧ v| ≡ 0 (mod 4)   for all X-stabilizers a and u, v ∈ stabilizers ∪ logical reps,
and induces a genuine order-4 logical ring-CCZ iff some logical triple overlap |l∧l'∧l''| is odd.  (C4) is
the mod-4 strengthening of Haah–Hastings' CCZ-triple condition.  Results: punctured RM(1,5) gives the FT
ring-CCZ codes [[31,1,3]]₄ (γ 3.13) and [[30,2,2]]₄ (γ 3.91); CSS Hamming codes and RM(2,m≤7) puncture codes
violate (C4) even after dropping to one logical (triple intersections 2 mod 4).  By Ax, RM(1,m) satisfies
(C4) for m ≥ 5 and RM(2,m) needs m ≥ 13 (n = 8192) — the same "Haah–Hastings regime" as the T-like family.
**The Phase-5 γ<1 values are existence statements about codes, not protocols**: no FT ring-CCZ code with
γ<1 is known, and the residue bound ties lifted ring codes to binary codes with 4-divisible triple
intersections, which at n ≤ 10⁴ reach γ ≈ 1.2–1.3 at best.


## 9. Step 1 — Galois-ring evaluation codes (full alphabet): negative at m ≤ 3

`galoisring/evalcodes.py`: Reed–Solomon-type CSS codes over GR(4,2), GR(4,3) (Teichmüller points, all
stabilizer/logical degree bands, optional 2-torsion X layers, per-qudit sign choice {q, q∘neg} optimized
on the S1–S3 violation count, all 32 escaping single-ququart tables).  54 codes, **0 hosting hits**.  The
component-wise carrier conditions are not ring-algebraic (q acts on Z₄-components, not on GR elements),
so the Teichmüller structure that gives RS triorthogonality over F_p does not transfer.  Larger-alphabet
ring codes would need a q-compatible structure (e.g. codes whose symmetrized weight enumerators satisfy
the weighted mod-16 conditions — Type-II-like Z₄ codes) rather than Galois-ring evaluation.
