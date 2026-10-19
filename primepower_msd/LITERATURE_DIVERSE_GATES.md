# Literature review — what is necessary / known for low-overhead, diverse-gate MSD on Z_{2^k} qudits

*2026-10-19.  Complements `LITERATURE_MAP.md` (2026-06).  Scope: single-qudit universality for
prime-power dimension, the Clifford hierarchy for even d, the triorthogonality templates for T / CS / CCZ
distillation, ring (Z_4) stabilizer codes, and the field-collapse "2^k" MSD literature.  Our own
results that fix the framing are in `UNIVERSALITY_CHECK.md`.*

## 1. What "universal" requires in dimension d = 2^k

**Borda–Rincón–Galindo, "Quantum Universality in Composite Systems: A Trichotomy of Clifford
Resources", arXiv:2512.20787 (v3, Dec 2025).**  The key reference.
- For prime d every non-Clifford gate is universal with the Clifford group (primality is necessary).
- For d = p^m, m ≥ 2, the Clifford adjoint action on sl(d) is *reducible*: sl(d) = ⊕_{j} W_j with W_j
  spanned by Paulis P_(u,v), gcd(u,v,p^m) = p^j (Thm 30).  Universality ⇔ irreducibility + infiniteness
  (Thm 4/17); a diagonal gate T_ζ is "orbit-mixing" iff for each proper divisor u the Fourier transform of
  its difference data is nonzero at some unit (Def 31/Thm 32).  T_s = diag(e^{2πik/s}) is non-Clifford iff
  s ∤ 2d (even d), irreducible iff s ∤ d (Thm 33; the sufficient infiniteness bound is loose, s ≳ 19 at d=4).
- Permutations: the transposition (0 1) is universal for every p^m, m ≥ 2 (Thm 39; d = 4 checked in Magma).
- Coprime composite d: intra-qudit CNOTs between the factors suffice — the arithmetic (GF(2^k) = (Z_2)^k)
  collapse we reject.
- **Open problem stated by the authors: whether the trichotomy extends to the multi-qudit setting.**
  Our sector test `a2.py` is the n-qudit version of their orbit-mixing criterion for diagonal gates, and the
  m = 2, 3 catalogs show entangling diagonal escapers exist but are frame-products of single-qudit ones.

**Hostens–Dehaene–De Moor, quant-ph/0408190.**  Pauli/Clifford bookkeeping for arbitrary d: phases are
powers of ζ = √ω (mod 2d) — necessary exactly when d is even and aᵀUa is odd (order-2d Paulis); Clifford
= (C ∈ Sp(2n, Z_d), h ∈ Z_{2d}^{2n}) with the parity constraint linking h to C.  For d = 4 they note
non-graph stabilizer states such as (|0⟩+|2⟩)/√2 with m > n generators — precisely the joint eigenbasis of
X², Z², i.e. the A₂ basis whose normalizer obstructs universality.

## 2. Clifford hierarchy for even / composite d (what is classified, what is open)

- Cui–Gottesman–Krishna 1608.06596: complete classification of *diagonal* hierarchy gates for prime d;
  de Silva–Lautsch 2501.07939: one qubit or qudit, odd prime.  "Groups of diagonal gates in the Clifford
  hierarchy" (Theor. Math. Phys. 2024): level-2 diagonal groups in arbitrary d, primes 2 and 3 special.
- 2502.01864 (addressable transversal non-Clifford gates), appendix: the full hierarchy for 2^t-dimensional
  qudits is **open** (their gates are GF(2^t)-field constructions, i.e. qubit collapses).
- Semi-Clifford structure: Beigi–Shor (every level-3 gate on 1–2 *prime* qudits is semi-Clifford);
  Gottesman–Mochon (3-qubit non-semi-Clifford level-3 gate).  Composite d: **open** — the only loophole for
  an entangling-only route to universality (probe E6).
- What we add: the finite-difference oracle for diagonal tables (`phase_tables.level`), exhaustive level-3
  catalogs on 1–2 ququarts (48 strict single-qudit tables; 130 048 strict 2-qudit tables), and
  W_m = level-3 module for m = 2, 3 (all diagonal level-3 gates are products of frame-rotated single-qudit
  gates and Cliffords).

## 3. Distillation templates (qubits and prime qudits)

- **Bravyi–Haah (2012)** triorthogonal matrices: |G_a ∧ G_b| even, |G_a ∧ G_b ∧ G_c| even for distinct
  rows, odd self-weights for logical rows; transversal T induces logical T.  [[3k+8, k, 2]] family, γ → log₂3.
- **Haah–Hastings, "Codes and Protocols for Distilling T, controlled-S, and Toffoli Gates", 1709.02832
  (Quantum 2018).**  *Generalized triorthogonal matrix* (Def. 1, Eq. 2.7): blocks G_T (T rows, odd weight,
  self-triple-overlap 1), G_CS (pairs with odd mutual overlap), G_CCZ (triples with odd 3-way overlap), G_0
  (even-weight X-stabilizers); all other triple overlaps even.  Transversal T then induces T / CS / CCZ on the
  respective logical blocks up to a Clifford.  RM(r,m), r = m/3 − 1 gives CCZ codes: m=3 [[8,·,2]] (28p²),
  m=6 [[64,2,4]], m=9 [[512,30,8]] (10 CCZ, 512 T's); randomized construction γ → 1.  This is the template
  our `rgt.py` generalizes: X-generators s with q-weight / pair / triple conditions S1–S3 over Z_4, with
  logical blocks read off from the induced table F(p).
- **Krishna–Tillich 1811.08461 (prime d).**  Tri-orthogonal matrix over F_p: Σ_i h_{ai}h_{bi} = 0 (a≠b),
  triply-even V∗V ⊆ V^⊥, and the normalization Σ_i h_{ai}³ = 1, Σ_i h_{ai}² = −1 (mod p) for logical rows —
  conditions that use inverses in F_p and do not transfer to Z_4 (zero divisors), which is why the ring case
  needs the table/finite-difference formulation.  Reed–Solomon: γ = O(1/log p).
- Campbell–Anwar–Browne PRX 2 041021 (odd prime RM protocols); Prakash–Saha 2510.10852 (Hastings–Haah
  γ ≈ 0.678 for all primes; the `qmsd` search substrate); low-overhead qutrit 2403.06228.
- **Watson–Campbell–Anwar–Browne 1503.08800, qudit colour codes:** modular-qudit colour codes in μ spatial
  dimensions carry a transversal level-μ gate when μ! is a multiple of... (the construction needs the lattice
  factorial condition; for d = 4 the 3D code's level-3 gate is obstructed and 4D is needed).  The one
  published family of *ring*-flavoured codes with transversal non-Clifford gates — a cultivation substrate.

## 4. Stabilizer codes over rings with zero divisors

- Gluesing-Luerssen–Pllaha 1710.09884 (FFA 2019): stabilizer codes over local Frobenius rings; a FREE code
  over the ring has distance ≤ that of its residue-field reduction in all studied cases (conjectured
  general).  ⇒ Z_4 lifts of binary codes inherit at most the binary distance (`lifts.residue_distance_bound`).
- Nadella–Klappenecker (Frobenius rings, ISIT 2012); Gunderman 2303.17000 / 2501.04888 (exotic and
  beyond-integral-domain local dimensions); EC-Zoo "modular-qudit stabilizer code": composite-d codes contain a
  prime-qudit subcode per prime factor bounding the distance, and need not encode an integer number of qudits
  (our half-integer k / "soft" directions).
- Hammons–Kumar–Calderbank–Sloane–Solé (1994): Z_4-linearity of Kerdock/Preparata, ZRM(r,m) codes — the
  classical substrate for `lifts.zrm`.

## 5. "2^k-dimensional" MSD in the literature is the field collapse

2512.21874 (Lamm et al.: F_{2^{2m}} AG triorthogonal, transversal CCZ, q as small as 64), 2408.07764
(constant-overhead MSD with q = 2^10 AG codes), 2608.09727 (Gong–Pattison–Rall–Wills: GF(2^s) codes →
qubit CS/CCZ distillation; 4 CS → 1 CS at distance 2), 2502.01864 (addressable CCZ via GF(2^t) diagonal
gates).  All treat the qudit as GF(2^k) ≅ (Z_2)^k — qubits in disguise (BRG's coprime/arithmetic regime).
**No published MSD or triorthogonal code over the cyclic ring Z_{2^k}; no qudit magic-state cultivation.**

## 6. Consequences for this program (see UNIVERSALITY_CHECK.md)

1. The universality certificate for a ring protocol is sector escape of the *induced logical* gate, not
   level-3-ness or odd third differences.  Neither the Phase-5 ring-CCZ nor the `[[8,2,2]]` carrier logical
   passes it.
2. Entangling targets (ring-CS, ring-CCZ, across-block CS-like escapers) are legitimate distillation
   targets but are exactly synthesizable from Clifford + single-qudit level-3 gates (W_m = level-3 module),
   so they are an efficiency choice — identical to the qubit situation (CS = 3 T, CCZ = 7 T).
3. The constructive route is the Haah–Hastings template over Z_4 with the S1–S3 screens and the parity
   transmission condition (`rgt.axis_parity_report`): codes whose logical-X axes have an odd number of odd
   entries.
