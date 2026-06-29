# Attack Vectors toward Low-Overhead Single-Qudit Prime-Power MSD

*Status 2026-06-29.  An 8-vector workflow (17 agents) with adversarial verification, + independent
re-verification of the load-bearing result.  This is a MENU of directions, not a resolution.*

## ⚠️ Headline correction (and a follow-up correction, both verified)

> **The "Z₄ zero-divisor ⇒ distance d=2" barrier (`M4_FINDINGS.md`) is a CSS / diagonal-transversal
> ARTIFACT, not ring-intrinsic.**  Dropping CSS and using the full **non-CSS symplectic** stabilizer
> formalism over Z₄ gives **distance-3 codes** — verified by full centralizer enumeration: a concrete
> non-CSS Z₄ code with **zero weight-1 logicals** (the 2·eₚ mechanism only fires on all-even columns).

> **🔧 FOLLOW-UP CORRECTION (this is a k-count fix): the verified code is `[[5,1,3]]₄` (k=1), NOT
> `[[5,2,3]]₄`.**  For a stabilizer code |C|/|S| = d^{2k}, so |C|/|S| = 16 = 4² ⇒ **k=1** (code dim =
> dⁿ/|S| = 4⁵/4⁴ = 4 = d¹).  Hence **γ = log(5/1)/log(3) ≈ 1.46, not 0.83.**  Since γ<1 needs n/k<d
> and here k=1 forces γ = log(n)/log(d) ≥ 1, **breaking the d=2 cap does NOT by itself give γ<1** —
> γ<1 still requires **k≥2 with d≥3** (unsettled for non-CSS Z₄).

> **V7 first result (validated even-d machinery, this turn): the verified [[5,1,3]]₄ admits NO
> transversal single-qudit level-3 gate** (0 of 780 tested: the 12 diagonal gates + semi-Clifford
> C₁DC₂).  A clean null for *this one* code.  (A broader multi-code search TIMED OUT on slow n=6
> dense codespace builds — inconclusive; so the V7 evidence is this single-code null + the structural
> argument that random non-CSS d≥3 codes lack the triorthogonality-analog structure, NOT a completed
> broad search.)  Machinery now correct: the even-d −I obstruction (X^aZ^b has order 8 when a·b is odd)
> is fixed by the **signed stabilizer** S′ = e^{−iπ(a·b)/4}X^aZ^b (validated; `tests/test_noncss_v7.py`).

> **V1 (GF(4)) is the field/2-qubit COLLAPSE — confirmed dead for the single-qudit goal.**  GF(4)'s
> additive group is (Z₂)², so its Pauli group is the 2-qubit group (shift order 2, fails anti-collapse),
> and diag(1,1,1,i) = controlled-S is literally a 2-qubit entangling gate.  It is the known field route,
> not genuine single-(cyclic-Z₄)-qudit magic.

**Net:** the d=2 cap is genuinely a CSS artifact (real result), but γ<1 single-qudit MSD is **still open
and still hard** — it now requires a non-CSS Z₄ code with **k≥2, d≥3, AND a transversal level-3 gate**
(the triorthogonality-analog), none of which the random codes provide.  V1 (field) is closed.

Concrete verified witness (full brute-force, min symplectic weight in centralizer\stabilizer = 3):
```
[[5,2,3]]_4 non-CSS Z_4 stabilizer code, |S|=256, |C|=4096
  stab (a|b):  [1,2,1,1,2] | [1,0,3,3,1]
               [1,1,1,2,1] | [2,1,3,0,2]
               [3,2,0,1,2] | [0,1,1,1,0]
               [2,1,2,0,2] | [2,1,3,1,3]
  non-CSS (mixed X/Z), distance 3, no weight-1 logical.
```

## Portfolio table (verdicts AFTER adversarial verification)

| Vector | Escapes barrier how | Verdict | Effort | Decisive next experiment |
|---|---|---|---|---|
| **V7 — non-CSS / non-diagonal ring stabilizer codes** | leaves the diagonal-CSS reduction; 2·eₚ becomes a *detectable error*, d≥3 codes exist (verified) | **LIVE (promising)** | medium | do the found d≥3 non-CSS Z₄ codes admit a transversal strict-level-3 (semi-Clifford C₁DC₂) gate inducing a non-Clifford logical? |
| **V1 — GF(4) field + separable (non-CCZ) gate (the hinge)** | field has no zero divisor (γ<1, d≥3); CS=diag(1,1,1,i) is a genuine non-Clifford level-3 GF(4) magic gate; strict triorthogonality forces *separable* T^⊗k | **LIVE (uncertain)** — proposal's "dead" REFUTED | medium | apply transversal CS-type gate to a small γ<1 GF(4) punctured-RM code; verify induced logical is non-Clifford, separable, AND non-factoring; compute distance |
| **V8 — magic-state cultivation (decouple gate from distance)** | d=2 gate confined to a tiny injector; distance grown by a high-distance Z₄ storage code (HGP: [[13,1,3]]₄,[[25,1,4]]₄,[[41,1,5]]₄, verified cyclic-Z₄ logical) | **LIVE (uncertain)** — wins on a *different metric* | medium | analytic first: does the lattice-surgery merge seam share the injector's all-even columns (→ weight-2 joint logical)? if not, simulate the gauge-fix + accepted-state overhead |
| V2 — larger d=8,16 | doesn't escape: order-2 zero divisor d/2 exists in every Z_{2^k} | **dead** (confirmed; cap still 2 at d=8, 0 codes with d≥3) | low | — |
| V3 — Galois ring GR(4,m) | doesn't escape: char-4 ⇒ nilpotent 2, 2²=0 verbatim | **dead** (1 honestly-flagged untested point: all-even-column *forcing* over GR never computed) | low | (if revisited) GR(4,m) triorthogonal-CSS distance + whether residue field still forces an all-even column |
| V4 — code-switching / gauge-fixing | only defeats Eastin-Knill; T still transversal in one gauge, capped at d≤2 | **dead** (residual non-CSS crack folds into V7) | low | — |
| V5 — addressable / non-uniform transversal | attacks only the rate ceiling; distance is gate-independent | **dead** — proposal's "promising" REFUTED | medium | — |
| V6 — CCZ→T catalysis | source distance from field code, inject phase at end | **dead** — P=diag(1,1,1,−1) provably unreachable from linear CCZ deposits; both gates are ordinary 2-qubit Cliffords in the CCZ frame | medium | — |

## The three live vectors (idea · enabling fact · first experiment)

**V7 — non-CSS ring stabilizer codes (sharpest, ring side; the hard part is already done).**
A semi-Clifford gate C₁DC₂ (C₁≠C₂) need not preserve a CSS code, so the "all-even column ⇒ 2·eₚ logical"
engine never fires.  *Enabling fact (verified):* d=3 non-CSS Z₄ codes exist at n=5,6,7 (full symplectic
formalism: form a·b′−a′·b mod 4, centralizer via the Howell right-kernel, distance = min symplectic weight
in C\S).  *First experiment:* extend the transversality check to general G=C₁DC₂ (explicit stabilizer
conjugation, small n) and run it on the already-found d≥3 codes.  Coexistence of d≥3 with a transversal
level-3 non-Clifford logical = the cap is broken (even at γ>1 it is a structural breakthrough + a scaling
target toward the verified γ≈0.83 parameters).  No coexistence across all small d≥3 codes = the no-go
upgrades from "CSS class" to "any ring stabilizer code."

**V1 — GF(4) hinge (cleanest route to actual γ<1, field side).**
A strictly-triorthogonal punctured-RM code over GF(4) already has γ<1 and d≥3 (no zero divisor).
*Enabling facts (both newly verified, two myths fell):* (1) for any strict-triorthogonal field code the
logical cross-tensor is C_abc = −[a=b=c], so the induced gate is **automatically separable** T^⊗k (no
logical CCZ); (2) **CS = diag(1,1,1,i) is a genuine non-factorizing non-Clifford level-3 GF(4) magic gate**
(the "additive order 2 ⇒ no genuine d=4 magic" argument is a non-sequitur).  *First experiment:* apply the
transversal CS-type gate to one small γ<1 GF(4) triorthogonal code; verify codespace-preserving + induced
logical non-Clifford + separable + (crucially) does NOT silently factor into two internal-qubit gates.

**V8 — magic-state cultivation (low-risk hedge, different metric).**
Confine the d=2 transversal constraint to a tiny injector ([[8,1,2]]₄, where d=2 is fine because injection
error is post-selected), then grow logical distance via a high-distance Z₄ storage code needing no
transversal T.  *Enabling fact (verified):* Z₄ HGP codes reach unbounded distance with a genuinely cyclic
Z₄ logical.  *First experiment (near-free, analytic):* check whether the merge seam shares the injector's
all-even columns.  *Caveat:* wins on qubits·rounds/accepted-state, NOT on γ=log(n/k)/log d.

## Recommended sequencing

1. **V7 first** — its hard sub-result (working symplectic code, concrete d≥3 codes) is already done; the
   remaining test (transversal level-3 on those codes) is small, finite, and binary.
2. **V1 in tandem** — the field-side mirror; both verified sub-results (separability theorem + genuine CS
   gate) just need composing on one real γ<1 code.  Cleanest path to actual γ<1.
3. **V8 as the parallel hedge** — near-free analytic gate first; delivers a usable protocol even if V7+V1 fail.

V7 (ring) and V1 (field) test the SAME reframed thesis from two sides: *distance ≥ 3 and a genuine single-
qudit transversal magic gate can coexist once you drop diagonal-CSS*.  Together they will settle whether
low-overhead single-qudit prime-power MSD is reachable or the no-go is genuinely structural.
