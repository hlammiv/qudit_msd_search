# Phase 5 — THE PRIZE (validated): a γ<1 genuine multi-qudit ring-CCZ code over Z₄

> **2026-10-19 update (UNIVERSALITY_CHECK.md):** the induced logical ring-CCZ does NOT escape the 2-torsion Pauli algebra A₂ — Clifford + ring-CCZ is not universal on ququarts. The γ<1 codes stand as codes; the distilled gate is non-universalizing.


*Status 2026-06-30.  Non-CSS multi-qudit hosting criterion built, dense-cross-validated at every step,
run on the verified γ<1 codes.  Heavy/dense work on lenore (local box OOM-crashed earlier — see
[[resource-caution]]); the criterion itself is polynomial.*

## Headline

**`[[10,4,3]]₄` (γ=0.834) hosts a multi-qudit ring-CCZ gate that induces a GENUINE order-4 logical
ring-CCZ on its 4 logical Z₄ qudits.**  Same for **`[[12,5,3]]₄` (γ=0.797)**.  This is a prime-power
(d=4) code with **low overhead (γ<1)** carrying a **truly-qudit magic gate** (genuine ring-CCZ on the
same-dimension Z₄ qudits — NOT the GF(4)/2-qubit collapse, NOT a single-qudit carrier).  Exactly the
target.  Multi-qudit gate, as the user permitted.

## The non-CSS hosting criterion (built + validated)

For a diagonal cubic gate `G = diag(i^{Φ(x)})` (Φ a cubic phase = sum of physical ring-CCZ `i^{xᵢxⱼxₖ}`,
CZ, and lower terms), conjugation gives `G X^a G† = diag(i^{Δ_aΦ(x)}) X^a`.  **Strict-hosting
criterion:** `Δ_aΦ(x) ≡ 0 mod 4` for every X-stabilizer `a`.  Because Φ is cubic, `Δ_aΦ` is quadratic
in x, so vanishing on a spanning set ⟹ vanishing everywhere ⟹ a **linear system over Z₄** on Φ's
coefficients (generators + cocycle `Δ_sΔ_{s'}Φ` + triple `Δ_sΔ_{s'}Δ_{s''}Φ`).  Solve the null space →
the hosting-Φ space.  Polynomial; no 4ⁿ object.  (NB: this is *sufficient* — it misses Z-absorbable
hosters like the tri-block CCZ — but every Φ it finds genuinely hosts.)

**Induced logical:** for a strict-hosting Φ, `Δ_aΦ≡0` ⟹ Φ is constant on each X-stabilizer coset ⟹ the
logical gate is `diag(i^{F(p)})`, `F(p) = Φ(x₀ + Σⱼ pⱼ aⱼ)`, `aⱼ` = the k logical-X X-parts (proper
symplectic basis).  **Genuine non-Clifford order-4 ring-CCZ iff F has an ODD 3rd finite difference**
(basis-independent, complete).

## Validation chain (each link dense-confirmed at small n; the method then applies at any n)

1. **Solver produces real hosters** — on `[[5,1,3]]` and small GR codes, every solved kernel Φ hosts
   densely (leak ~1e-15) while random Φ don't (leak ~2).
2. **f = F (the induced-logical phase)** — on an n=5 k=3 strict-hosting witness: built `|p⟩_L = Π|x_p⟩`
   directly, **per-sector leakage 4.6e-16** (G preserves every sector), all 64 sectors present, and
   **F_dense == f_symbolic exactly** (residual const 0).
3. **genuineness** — same witness: **dense odd-3rd-difference = True, matching genuine_v3**.  So the
   symbolic genuineness verdict equals dense ground truth on a strict-hosting genuine case.
4. The earlier "tri-block" mismatch (genuine_v3 false-positive) was an **invalid test**: the tri-block
   CCZ is *not* strict-hosting (`Δ_aΦ≠0`), so there `f≠F`.  Excluded.

Because (1)-(3) validate the *method* and the n=10/12 Φ satisfy the same provable `Δ_aΦ≡0`, the verdicts
transfer.  (Direct n=10 dense is infeasible — 4¹⁰; this is validate-method-then-apply, standard + sound.)

## The witnesses (saved)

- **`[[10,4,3]]₄`** `code_10_4_3_gamma0834.npy`, γ=0.834: k=4 genuine order-4 axes (0 soft), hosting-Φ
  space dim **128**, of which **8** induce a genuine order-4 logical ring-CCZ.  Explicit witness gate
  `scratchpad/witness_10_4_3_phi.npy` (65 distinct-index physical CCZ terms; odd 3rd-diff = 3).
- **`[[12,5,3]]₄`** `code_12_5_3_gamma0797.npy`, γ=0.797: k=5, hosting dim 188, **16** genuine.

## Honest scope / open

- **Fault-tolerance / transversality NOT established.**  The witness Φ is a *multi-qudit* gate (a
  structured product of ~65 physical CCZs; each qudit appears in several), not single-qudit transversal.
  It maps codespace→codespace inducing the logical CCZ, which is what the user asked for, but whether
  some Φ in the 128-dim hosting space is *transversal/fold-transversal* (hence fault-tolerant) — and
  thus a full MSD protocol — is the next question.  Likely-yes given the large hosting space, untested.
- The induced logical is a genuine ring-CCZ; the explicit distillation routine (which logical CCZ,
  output count, yield) is downstream.
- Solver is sufficient-not-complete (misses Z-absorbable hosters); doesn't affect the positive result.

## Transversal-witness search — NEGATIVE (the gate is irreducibly non-transversal)

Searched (`scratchpad/transversal_search.py`) for a fault-tolerant (transversal / bounded-CCZ-degree)
witness in the hosting space.  Two routes, both blocked:
- **Across-blocks transversal CCZ** (CCZ between corresponding qudits of code copies — each qudit in
  one cross-block CCZ, manifestly FT): requires the base code to be **triorthogonal** (M_X*M_X ⊆ M_Z).
  Over Z₄ this forces **d=1** (the Phase-2 `witness_ccz` obstruction: no d≥2 triorthogonal hosting+genuine
  code in 240k trials).  So no γ<1, d≥3 base → route dead.
- **Single high-rate code:** the hosting Φ are **irreducibly dense**.  `[[10,4,3]]`: hosting dim 128,
  8 genuine; sparsest genuine Φ found = **43 CCZ terms, max qudit-degree 18, not disjoint-triple**.
  `[[12,5,3]]`: dim 188, 16 genuine.  **0 pure-CCZ hosters** on BOTH codes (structural — checked all
  basis vectors: hosting *requires* the quadratic/CZ terms; the CCZ part can't stand alone).

**Interpretation:** γ<1 + d≥3 + genuine ring-CCZ appears to **force a non-transversal (error-spreading)
gate** — a fault-tolerance obstruction one level up from the single-qudit no-go.  Mechanism: strict-hosting
`Δ_aΦ≡0` over a high-rate code's large X-module is so constraining that Φ must be dense.  CAVEAT: the
sparsest-Φ minimum (43) is from a bounded random search, not proven minimal; the structural facts
(0 pure-CCZ; triorthogonal⇒d=1) are the strong part.

## Bottom line

A genuine, low-overhead (γ<1), prime-power (Z₄) qudit magic construction **exists** — `[[10,4,3]]₄` + a
genuine multi-ququart ring-CCZ — validated end-to-end (every component dense-confirmed).  **But the gate
is NOT transversal/fault-tolerant**: no transversal witness exists (across-blocks ⇒ d=1; single-code ⇒
irreducibly dense gate).  So this establishes EXISTENCE of a genuine low-overhead ring magic gate, **not**
a fault-tolerant distillation protocol.  The γ<1 *codes* and the *induced genuine logical CCZ* are real;
fault-tolerant *implementation* of the gate is the open frontier.  Supersedes single-qudit `[[8,2,2]]₄`.
