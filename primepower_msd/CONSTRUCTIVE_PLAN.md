<!--
Generated 2026-06-29 by a 15-agent constructive-planning workflow (design -> adversarial critique -> plan).
PI-VERIFIED addenda (independent checks this turn):
 * Centralizer/engine bug flag (Critique 1: "7/8 odd-X centralizer rows fail commutation") = REFUTED.
   Audited [[5,1,3]]_4: 6 centralizer generators, all 4096 centralizer elements commute with every
   stabilizer (0 failures), stabilizers isotropic, distance=3 confirmed. The non-CSS d=3 result HOLDS.
 * omega=-1 "soft collapse": our verified [[8,1,2]]_4 gate diag(1,-i,1,i) IS genuine level-3 +
   passes basic anti-collapse (order-4 clock), but its magic content is Z_2-sublevel (omega=-1), i.e.
   weaker than a full ring resource. The hardened certificate (logical commutator phase = i, omega=i)
   is the stronger test adopted in Phase 0a. This refines, does not refute, the existence result.
 * Single-qudit diagonal cubic T at d>=2 is OBSTRUCTED (Critique 1 proof) -> reinforces the V7 null
   with a mechanism; do not chase single-qudit diagonal T into the distance frontier.
-->

# Constructive Research Plan: A Genuine Ring-Qudit Magic Resource over Z_4

## 1. Executive Summary

Seven designs were proposed; the adversarial critiques killed two outright, gutted the central mechanism of three more, and converged — independently, from three different starting points — on a single structurally correct lever. The honest state after reconciliation:

**The two most promising constructive approaches (post-critique):**

1. **Galois-ring GR(4,r) Reed–Muller, non-CSS, targeting a genuine ω=i multi-ring-qudit CCZ (target B).** This is the *merge point* of the salvaged Design 5 and Design 6. Both critiques independently demanded the same pivot: abandon full-ring Z_4^m evaluation (degenerate Gram, Critique 5) and abandon welding-a-fixed-d=2-seed (no distance lever, Critique 6), and replace both with polynomial evaluation over the Galois ring GR(4,r), whose trace form Tr: GR(4,r)→Z_4 is **non-degenerate**. This is the *only* surviving proposal with a genuine distance-scaling lever (polynomial degree → d→∞, the Saha–Prakash/Bravyi–Haah mechanism) on a form that is not 2-torsion-degenerate, while keeping single coordinates as genuine order-4 Z_4 qudits.

2. **Fixed-gate constraint-solve on a fully unit-supported module (salvaged Design 4).** Keep the genuine parts — DXD† = X·S with S a genuinely-quadratic ζ_16 phase, and the reduction of transversality to a finite phase/module condition — but **delete the 2-torsion "free lunch"** (proven to produce order-2 = qubit logicals). This is the tractable near-term *witness generator*: it can produce a small [[n,1,3]] or [[n,2,3]] unit-supported code to falsify or confirm the vacuity question before any RM lifting.

**Single recommended starting point:** Do **Phase 0** first — it is cheap, finite, and gates everything: (0a) harden the anti-collapse certificate to test the **logical commutator phase**, not substrate clock order; (0b) run the finite **gate-existence enumeration** to settle whether an ω=i carrier even exists; (0c) run the **vacuity existence probe** for any small non-CSS d≥3 code hosting a level-3 anti-collapse logical. Then commit Phase 1 to **GR(4,r) Reed–Muller (approach 1)** as the primary scaling bet, with approach 2 run in parallel as the witness generator.

The deprioritized/killed designs are §5.

---

## 2. The Multi-Qudit-CCZ Resolution (decided)

**Decision: the multi-qudit CCZ is a BLESSING — but only the ω=i version. Target B (a genuine multi-ring-qudit CCZ) is the primary near-term deliverable; single-qudit T (target A) is the optional downstream via catalysis.**

Justification, reconciling the designs against their critiques:

- **Single-qudit T via a diagonal level-3 gate is structurally obstructed at the distance frontier.** Critique 1 proves it: forcing the cubic diagonal gate diag(i^{x³}) transversal *requires* even X-parts (the x²-coefficient 3a is non-Pauli unless a is even), and even-X stabilizers at d≥2 force every Z_i² into S, which forces every logical X-part even, which forces the logical symplectic form **even** ⇒ logical commutator i^{even} ∈ {±1} ⇒ the forbidden (Z_2)² collapse. Distance and anti-collapse are **mutually exclusive** under that gate. So we do not chase a single-qudit diagonal T into the d≥2 region.

- **The CCZ is the natural output of any triorthogonal-type construction** (the user's framing is correct), and it can be made genuinely ring. The carrier's cubic finite-difference is what produces it.

- **The hard line within "embrace CCZ": ω=i, not ω=−1.** Critique 3 is decisive: the *verified* gate t = diag(1,−i,1,i) has a natural triple 3rd-difference = 2 mod 4 ⇒ i² = −1, yielding diag((−1)^{x₁x₂x₃}) — **the qubit CCZ embedded in Z_4 qudits**. The code qudits stay order-4, but the *magic content lives only in the Z_2 sublevel* — a real but strictly weaker resource and a soft-collapse. A genuine ring CCZ requires the triple 3rd-difference to be an **odd unit of Z_4 (value i)**, which needs an **order-8 (ζ_8) phase carrier with 3rd-difference = 2 mod 8**. Whether such a carrier that is *simultaneously* single-ququart level-3 + T-condition + pairwise-Clifford exists is a **finite question** (Phase 0b) and must be settled before anything else.

- **Conjugation does not produce a CCZ.** Critique 2 kills Design 2's framing: Clifford-conjugating a code+gate that induces T⊗T leaves a logical action that is Clifford-*equivalent* to T⊗T; CCZ is not Clifford-equivalent to T⊗T, so the "ring-CCZ" was a basis artifact. A genuine CCZ must be built on **one** isotropic module that the entangling level-3 gate actually normalizes — never by relabeling two T-blocks.

- **CCZ→T conversion survives on the ring.** The prime-route CCZ→T no-go was a *field artifact* (GF(4) trace-CCZ has phase order 4, T has order 4 — and Critique 3/7 corrected the arithmetic: the *ring* CCZ phase has order 8, T has order 4, and 8 is a multiple of 4, so catalysis 8→4 is order-compatible). So target A remains reachable downstream, but we do not make the program depend on it.

**Target hierarchy:** (B) genuine ω=i multi-ring-qudit CCZ, certified by the hardened logical-commutator certificate → primary. (A) single-ring-qudit T via Z_4 catalysis/code-switching of the certified CCZ → optional follow-on. We **force single-qudit T directly only** in the purely-quadratic-phase, no-cross-term regime — which is distance-capped, so it is not a scaling route.

---

## 3. Promising Approaches: Recipe, Key Math, Milestones, Decisive Experiment

### Approach 1 (PRIMARY): GR(4,r) Reed–Muller, non-CSS, ω=i ring CCZ

**Why this and not its two parents:** Design 5 (full-ring Z_4^m evaluation) is dead as posed — Critique 5 verified the Gram B = E·Eᵀ ≡ 0 mod 4 for all m≥2, making self-adjointness *vacuous*, isotropy *automatic* (selects nothing, hence the universal distance-1), and the level-3 bookkeeping (which lives mod 8) *invisible* to a mod-4 form. Design 6 (welding) is dead — Critique 6 showed a fixed d=2 seed has no distance lever (+1 per weld while n grows comparably ⇒ d bounded ⇒ n/k can't fall below 3) and the weld-distance lemma is unsound (merging coordinates generically caps/reduces distance). **Both critiques pointed to the same fix:** evaluate over GR(4,r), whose trace form is non-degenerate, restoring meaningful Gram, RM distance scaling with degree, and a non-vacuous isotropy condition — while single coordinates remain genuine order-4 Z_4 qudits (GR(4,r) additive structure is (Z_4)^r, order 4, **not** the order-2 GF(4) collapse).

**Construction recipe:**
1. Fix GR(4,r); evaluate Z_4-polynomials/monomials at GR(4,r) points (Hammons–Kerdock–Preparata lineage). Single code symbols are Z_4.
2. Build the non-CSS symplectic generator rows (a|b) = (ev(f) | ev(M f)) with M **B-self-adjoint w.r.t. the non-degenerate trace Gram** B(f,g)=Tr Σ f g. Non-scalar M ⇒ genuinely non-CSS (escapes the diagonal/CSS d=2 artifact).
3. STABILIZER = high-degree band (so logical cosets are low-degree, few, high-weight — distance from the RM degree gap, now *real* because the Gram is non-degenerate). LOGICALS = centralizer/stabilizer.
4. CARRIER = the ω=i order-8 phase gate selected in Phase 0b, applied transversally. Codespace preservation = the ring-symplectic triorthogonality condition (below).

**Key math to derive (the load-bearing deliverable — currently an IOU; Critique 7 is right to flag it):**

The **ring-symplectic triorthogonality condition** for a diagonal carrier U = diag(ζ_{2^j}^{f(v)}), via the layered finite differences of f over (Z_4)^k — and it must hold over **all of (Z_4)^k, not the {0,1}^k corner** (Critique 3's key flaw: Z_4 carries at v=2,3 generate uncontrolled higher terms the binary corner misses):

- **Linear / T-layer** Δ_a f: must vanish (mod 2^j) / land in S → the "T-condition" (Gauss-sum/weight identity on unit-support coordinates).
- **Bilinear / pairwise layer** Δ_a Δ_{a'} f: must reduce to a logical Clifford (CZ-admissible) — this layer couples **both a and b** of generator pairs (the non-CSS symplectic twist; absent in CSS triorthogonality).
- **Trilinear / cubic layer** Δ_a Δ_{a'} Δ_{a''} f: must be **0 mod 2^j on every triple that involves a stabilizer** and **an odd unit (= value i, requiring 3rd-difference = 2 mod 8) on the pure-logical triple**. This is the Z_4/Z_8 generalization of the Bravyi–Haah GF(2) triple-product condition.

Two facts to honor: (i) **work mod 8 minimum** — the level-3 induced-Clifford bookkeeping is invisible mod 4 (Critiques 3, 5). (ii) The **CSS Bravyi–Haah-Z_4 carry route is the one demonstrably-working seed** for the cubic layer ("odd number of positions where three X-supports are all units mod 4," Critique 3) — but CSS caps d=2 over Z_4, so the central theoretical risk is bridging this guarantee onto the **non-CSS** GR code that supplies distance. **Deriving this bridge is Milestone 1.**

**Milestones:**
- M1: Derive and *prove* the non-CSS ring-symplectic triorthogonality condition at mod-8 precision over GR(4,r), reproducing the verified k=1 codes as a sanity check.
- M2: Prove a GR(4,r) union-of-2-adic-supports distance bound that *beats d=2 at k≥2* (the explicit open frontier; if it cannot, the family is dead — pivot carrier).
- M3: Construct a concrete [[n,2,3]] or [[n,3,3]] code hosting the ω=i carrier.

**Decisive first experiment (after Phase 0):** at the smallest GR(4,r) with a non-degenerate trace Gram (start r small, n≈16–64), build the high-degree-band non-CSS module, apply the Phase-0b ω=i carrier via `transversal_test` (leakage 0), extract the induced k-qudit logical, and require: (i) `is_level3_n` True, `is_clifford_n` False; (ii) the diagonal-phase **3rd finite-difference = odd unit of Z_4 over all (Z_4)^k** (genuine ω=i CCZ, not ω=−1); (iii) **hardened anti-collapse**: logical X̄–Z̄ commutator symplectic form ODD for every logical pair; (iv) d≥3 by bounded symplectic-weight enumeration on centralizer\S. PASS on (n, k≥2, d≥3) ⇒ first γ<1 genuine ring-qudit magic code. KILL signal: M ranges over all self-adjoint twists but the distance band never breaks d=2 at k≥2 (M2 fails) ⇒ GR(4,r) carrier insufficient.

### Approach 2 (PARALLEL — witness generator): Fixed-gate constraint-solve on a fully unit-supported module

**Salvaged Design 4.** Keep: D(1,0) is the correct magic generator (DXD† = X·S, S genuinely quadratic ζ_16); transversality reduces to a finite phase/module condition (probe-verified diagonal defect). **Delete the 2-torsion free-lunch** — Critique 4 verified that a 2-torsion logical X-support (a_i∈{0,2}) has clock order **exactly 2** ⇒ that logical *is a qubit* ⇒ fails anti-collapse, and carries no ζ_16 magic. So magic + anti-collapse ⇒ **every logical must be unit-supported**, which costs symplectic weight identically to the prime route — there is no free overhead win.

**Recipe (honest version):** require all logical X and Z unit-supported; solve the linear self + bilinear layers by GF(2)-Gaussian **with bounded backtracking** (the "Hensel-is-Gaussian" claim oversells — the cubic/decoupling layer gives quadratic obstructions and the variety can be singular/empty; Critique 4). Before solving, **derive the honest ring weight/Singleton bound for unit-supported logicals** — the real skipped question of whether (k≥2, d≥3) is feasible at all. Near-term deliverable: host D on a small unit-supported [[n,1,3]] (or push k=2, d=3) and **measure γ honestly** — even γ≥1 with a real transversal level-3 gate on a hardened-anti-collapse-passing code is a genuine structural advance past dead random search, and directly attacks the vacuity risk.

**Role:** This is the cheapest path to a *witness* that the constraint system {d≥3, transversal level-3, k≥2, anti-collapse} is non-empty, which Phase 0c also targets. Use it to feed/validate the M1 condition for Approach 1.

### Approach 3 (CONDITIONAL — single entangling module): salvaged Design 2 building block

Only the *building block* survives: CSUM(D⊗D)CSUM† is a genuine 2-qudit level-3 diagonal gate. Deploy it on **one** isotropic Z_4 module it actually entangles, solving the bilinear/trilinear cross-orthogonality (the T2 layer) **from scratch** over Z_4, verifying the induced logical against the CCZ phase polynomial (not against any conjugate of T⊗T). Keep this in reserve as an alternative carrier family for Approach 1 if the diagonal ω=i carrier from 0b does not exist — an *entangling* carrier is another way to source the cubic layer.

---

## 4. Recommended Sequenced Program

**Anti-collapse check at every step** = the *hardened* certificate (0a): logical X̄–Z̄ commutator symplectic form ODD (= phase i), **plus** for CCZ targets the triple 3rd-difference an odd unit over all (Z_4)^k. Substrate clock order 4 alone is necessary but **not sufficient** and must never be accepted as the certificate.

### Phase 0 — Cheap finite gates (do first; days)

- **0a. Harden the certificate + audit the engine.** Implement the logical-commutator-phase witness. Re-audit every prior "anti-collapse PASS." Brute-force cross-check `right_kernel_howell` centralizer outputs (Critique 1 flagged 7/8 audited odd-X "centralizer" rows failing explicit commutation recheck — a load-bearing bug that would invalidate any positive claim).
  - PASS: certificate distinguishes genuine order-4 logical commutator from {±1} collapse; centralizer machinery verified. KILL the program until fixed if the engine is unsound.

- **0b. Gate-existence enumeration (finite).** Enumerate q: Z_4→Z_8 (and →Z_16): find any diagonal single-ququart gate with level_of = 3, non-Clifford, **3rd-difference = 2 mod 8 (ω=i)**, vanishing linear (T-condition), pairwise Clifford.
  - PASS: ≥1 ω=i carrier exists → carry into Phase 1.
  - KILL: none exists → the genuine ω=i single-diagonal-carrier route is dead; the only honest single-carrier deliverable is the ω=−1 (Z_2-content) CCZ, which we **reject as a soft-collapse**. Pivot to Approach 3 (entangling carrier) for the cubic layer.

- **0c. Vacuity existence probe (brute force, small n).** Find ANY non-CSS Z_4 code, smallest n, even k=1, with d≥3 hosting a transversal level-3 gate (from 0b) inducing a level-3 logical passing hardened anti-collapse.
  - PASS: constraint system non-empty → invest in scaling (Phase 1).
  - KILL: empty at all small n → deep obstruction; reconsider whether non-CSS d≥3 can host *any* level-3 gate (would be a clean ring no-go worth publishing).

### Phase 1 — First construction (target B, k=2–3, d=3)

- Commit to **Approach 1 (GR(4,r) RM)** primary; run **Approach 2 (unit-supported constraint-solve)** in parallel as witness generator.
- M1: prove the non-CSS ring-symplectic triorthogonality condition (mod 8). Gate: must reproduce the verified k=1 codes.
- M2: prove the GR(4,r) 2-adic union-of-supports distance bound breaks d=2 at k≥2. **This is the make-or-break theorem.** If false, the GR carrier is dead.
- Run the decisive experiment (§3, Approach 1). PASS = first γ<1 genuine ring CCZ.
  - Pass/kill: leakage 0 AND ω=i over all (Z_4)^k AND hardened anti-collapse AND d≥3 ⇒ PASS. Any failure ⇒ diagnose which triorthogonality layer broke and feed back to M1.

### Phase 2 — Scale to k≥2, d≥3 ⇒ γ<1, then optional A

- Scale GR(4,r) (increase r and the degree band) to grow k and d simultaneously; track γ = n/k vs d honestly at each step with the hardened certificate.
- Optional target A: catalyze the certified ω=i ring CCZ → single-ring-qudit T via the Z_4 CCZ→T gadget (order 8→4 compatible). Construct, do not assert.

---

## 5. Deprioritize / Kill (with reasons)

- **Design 1 (even-X diagonal cubic gate diag(i^{x³})) — KILL. Collapse risk realized.** Critique 1's chain is essentially a proof: transversality *forces* even X-parts; even-X at d≥2 forces all Z_i²∈S, forcing even logical X-parts, forcing an **even logical symplectic form** ⇒ commutator ∈{±1} ⇒ the exact (Z_2)² collapse the project forbids. Genuine ring (odd form) survives only at d=1 — matching the design's own d=1 probe. Distance and anti-collapse are mutually exclusive here. *This was the most confident design ("medium risk, distance is the open lever") and the critique downgraded it to broken/deprioritize — the sharpest reversal in the set.* Salvageable lesson: the certificate must check the logical commutator, not clock order (→ 0a).

- **Design 2's conjugation framing — KILL (gives nothing).** Clifford-conjugation is a code equivalence; induced logical is Clifford-equivalent to T⊗T, never a genuine CCZ; and even granting d=3, n=16/k=2 ⇒ γ>1. Only the building block survives (→ Approach 3).

- **Design 4's 2-torsion "free-lunch" rate — KILL (collapse).** 2-torsion logicals have clock order exactly 2 = qubits = the (Z_2)² collapse. No free overhead. The constraint-solve survives only on **fully unit-supported** modules (→ Approach 2).

- **Design 5's full-ring Z_4^m evaluation — KILL (degenerate).** Gram ≡ 0 mod 4 ⇒ self-adjointness vacuous, isotropy automatic, distance universally 1, level-3 bookkeeping invisible. Pivot to GR(4,r) (→ Approach 1).

- **Design 6's welding — KILL (no distance lever; unsound lemma).** A fixed d=2 seed cannot generate a distance-growing family; the "d₁+d₂−deficit" weld lemma is false (merging coords caps/reduces distance); and target B is not easier than A (both need the same self/no-leakage condition). The one real asset — the genuine 4-torsion associator C∈{0,4,8,12} — is preserved by moving to GR(4,r) RM (→ Approach 1).

- **Designs 3 and 7 — not killed but reframed.** Design 3's ω=i claim was unverified (natural coupling gives ω=−1) and its non-CSS M2 mechanism is unsupported (k=1 self-phase is no evidence for 3-qudit cross terms); its *correct* content (the ω=i/ω=−1 distinction, the CSS-carry cubic route, the full-(Z_4)³ certificate) is absorbed into Phase 0b and the §3 math. Design 7's carrier-conflation (polynomial degree of the target phase ≠ arithmetic locality of the physical gate — in Bravyi–Haah the carrier is single-qubit T, cubicness comes from the *code*) is corrected; its strongest contribution — **front-load a vacuity existence probe before any RM lifting** — is adopted as Phase 0c, and its "lift-and-Howell-repair" step is replaced by the intrinsic mod-8 derivation (M1), because Howell canonicalization preserves neither triorthogonality nor distance.

---

### Bottom line

Start with **Phase 0** (hardened certificate → gate enumeration → vacuity probe): three cheap, finite, decisive tests that either unblock or kill the program before any heavy construction. The single recommended constructive bet is **GR(4,r) Reed–Muller, non-CSS, targeting a genuine ω=i multi-ring-qudit CCZ** — the convergent salvage of two independently-critiqued designs and the only route with a real distance lever on a non-degenerate form and genuine order-4 qudits. The central honest risk is **Milestone 2** (proving the GR distance bound breaks d=2 at k≥2) and the still-open **bridge theorem (M1)** carrying the CSS triorthogonality gate-guarantee onto the non-CSS GR module. Until 0a hardens the certificate and 0c exhibits a non-empty witness, treat every prior "anti-collapse PASS" — including ones in the verified record — as provisional.