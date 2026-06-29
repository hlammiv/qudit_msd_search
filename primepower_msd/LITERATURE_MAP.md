<!--
Generated 2026-06-28 by a 40-agent literature+planning workflow (~1.1M tokens).
All load-bearing claims were adversarially verified; see confidence flags / Verified references.
-->

# Literature Map: Single-Qudit Magic State Distillation for Prime-Power Dimensions d = p^e

*Scope: the open problem of a genuinely single-qudit (monolithic, cyclic Z_{p^e}) non-Clifford resource and a distillation code that produces it. Field GF(p^e) constructions are treated as the settled negative baseline. All load-bearing claims have been cross-checked against the adversarial verification verdicts; corrected or down-weighted claims are flagged inline.*

---

## 1. Executive summary

Single-qudit prime-power MSD is **plausible but unproven**. Nothing in the surveyed literature constructs it, and nothing rules it out. The situation as of mid-2026 is sharper than "open": the **target gate is now explicitly named** in the literature, while the **distillation code that would produce it does not exist**.

Two facts reframe the project:

1. **The target resource is identified.** Borda–Rincón–Galindo (arXiv:2512.20787, *verified*) prove that for d = p^m with m ≥ 2 the universal non-Clifford resource is a **single-qudit diagonal phase gate** T_s = Σ_k e^{2πik/s}|k⟩⟨k| whose phases live over the **ring Z/dZ**, not a multi-qudit CCZ. They prove universality (with Clifford) but give **no distillation protocol and no code** — that is precisely the gap.

2. **A correction to the precision target.** The verdicts confirm that the even-d Clifford phase gate is S = diag(ω_{2d}^{x²}), so on Z_4 the **8th-root quadratic phase diag(ζ_8^{x²}) is Clifford, not magic**. The prior project assumption ("8th roots give level-3 for d=4") is wrong. A level-3 single-ququart gate needs precision *above* the 8th-root quadratic: either a **16th-root quadratic** (diag ζ_16^{x²}) or a **non-additive cubic** (diag ζ_8^{x³}).

The single most promising route is **route R1 below**: a ring-triorthogonality analog over Z/p^m, evaluated over the *full ring* (not the Teichmüller set), engineered so that full-ring power-sum / Gauss-sum cancellations make a transversal **single-qudit T_s** appear. The pivotal enabling fact — verified by direct computation — is that **over Z_{2^k} (k ≥ 2) the squaring and cubing maps are genuinely non-additive**, unlike Frobenius over GF(2^e). This is the exact structural reason the field collapse need not recur. It is best pursued in tandem with **route R2**, a synthesis-first computational search (reusing the existing qmsd toolkit) that could settle existence with a single small example.

Honest caveat: two hard constraints bound any payoff. Free stabilizer codes over chain/Frobenius rings cannot beat residue-field codes on distance (so the win must be the *gate*, not the code), and the standard tensor/field Pauli group for p^r factorizes into r p-qudits (so the cyclic clock Pauli group is mandatory and is unclassified).

---

## 2. The precise open problem

**Field vs ring.** GF(p^e) is a field; Z_{p^e} is a chain ring with zero divisors and a cyclic additive group. The published prime-power MSD route is field-based and produces a structured **multi-qudit** trace-CCZ |x,y,z⟩ → ω^{Tr(xyz)}, ω a p-th root, which decomposes via a normal/trace basis into e smaller qudits (d=4 → two qubits). This is "Clifford in the dimension it lives in" and is *not* the target. The open problem requires the **ring Z_{p^e}** with a cyclic clock-shift Heisenberg–Weyl group.

**Single- vs multi-qudit.** "Single-qudit" means the distilled output is one monolithic d-level system and its non-Clifford resource is a single-system diagonal level-3 gate (Howard–Vala / BRG T_s), not an entangling gate across subsystems. Every recent transversal-non-Clifford construction (Section 3.4) yields a *multi-block* CCZ / C^{m-1}Z instead.

**What "magic" means for even / prime-power d.** The clean odd-prime-power resource theory (discrete Wigner function, Hudson theorem, mana, contextuality = negativity) **does not transfer to even d**: Zhu (arXiv:1504.03773, *verified*) proves no discrete Wigner function is Clifford-covariant in any even prime-power dimension (the permutation-symmetric / 2-design representation exists only for odd prime powers plus d=2 and **d=8**, the Hoggar-lines exception). So mana/Wigner-negativity are unavailable as the magic witness for d=4,8,16. Magic is nonetheless well-defined and distillable in even d (qubit MSD is the existence proof), but it must be quantified by **Heisenberg-Weyl-based monotones** (stabilizer Rényi entropy, robustness of magic), whose monotonicity is proven only for prime d and is **open for composite/ring d=4,8**.

**The required precision (corrected).** For cyclic Z_{2^k}, the level grading is set jointly by root-of-unity precision and polynomial degree (Cui–Gottesman–Krishna, *verified*, for the prime case). The even-d Clifford S already consumes 2d-th-root (= 8th-root for d=4) precision at quadratic order. Therefore a genuine level-3 single-ququart gate needs **16th-root quadratic precision or a non-additive cubic at 8th-root precision**.

---

## 3. State of the art, by theme

Notation per entry: **[field/ring/qubit/general]** · **[single/multi-qudit]**. ✅ = independently verified in the verification pass; ⚠ = scout-asserted only (plausible, not independently verified here); ❌ = attribution/ID problem, see flag.

### 3.1 Ring Clifford hierarchy and single-qudit non-Clifford structure

- **Borda–Rincón–Galindo, arXiv:2512.20787** ✅ — [ring] · [single-qudit]. For d=p^m (m≥2), universality with Clifford is achieved by the single-qudit diagonal T_s = Σ_k e^{2πik/s}|k⟩⟨k| over Z/dZ. *Correction:* T_s is non-Clifford iff s ∤ K_d with K_d = d (odd d) or 2d (even d) — so "Clifford iff s|2d" holds **only for even d** (the project's d=4,8,16); for odd prime powers it is s|d. Universality additionally needs the density bound s > ~1.57(d−1). **No distillation code.**
- **Cui–Gottesman–Krishna, arXiv:1608.06596** ✅ — [general] · [single-qudit]. Diagonal hierarchy gates for **prime** d are p^m-th roots of unity raised to a polynomial of the basis label; level set by m and degree. *Key obstruction:* for p^r with the **standard tensor/field Pauli group**, the hierarchy is isomorphic to r separate p-qudits — so single-qudit prime-power magic **cannot** live on the standard Pauli group; the cyclic Z_{p^r} clock group (unclassified) is mandatory.
- **de Silva–Lautsch, arXiv:2501.07939** ✅ — [general] · [single-qudit]. Completely solves the single-qudit Clifford hierarchy **only for prime d** (full hierarchy, not just diagonal). No prime-power / cyclic-ring claim — the gap is genuine.
- **Howard–Vala, arXiv:1206.1598** ✅ — [field] · [single-qudit]. Explicit single-qudit T for prime d via a cubic phase using 12^{-1} and 2^{-1} mod d; requires Z_d a field. The arithmetic obstruction at d=p^e is exact: **2 is not a unit in char 2** (12 = 2²·3), so the closed form is undefined for powers of two and three.
- **Rengaswamy–Calderbank–Pfister, arXiv:1902.04022** ✅ — [ring phase / qubit carrier] · [multi-qubit]. Level-k diagonal gates = quadratic forms (symmetric matrices) over Z_{2^k}. *Confirmed false-friend:* Z_{2^k} is the **phase ring of an n-qubit register**, not a single 2^k-level qudit; v2 erratum restricts to **2-local and certain higher-locality** gates. The correct ring-precision dictionary, but not a single-qudit theory.

### 3.2 Codes over rings (substrate, no magic gate)

- **Hammons–Kumar–Calderbank–Sloane–Solé, math/0207208** ⚠ — [ring] · [no gate]. Foundational Z_4-linearity (Kerdock/Preparata over GR(4,m)), Teichmüller/Gray-map machinery. The Gray map is exactly the collapse mechanism that turns a ring code into a binary multi-qubit code.
- **Andriatahiny et al., arXiv:1801.05114** ✅ — [ring] · [no gate]. Generalized Reed–Muller over GR(p^s,r): free module, dual law, distance all identical to the field GRM. *Confirmed:* evaluation is at the **Teichmüller set**, and the cancellation (Prop 7.1) is the **multiplicative geometric series** Σ_i ξ^{ij}=0, *not* a genuine Z_{p^s} additive power-sum. This is why a naive lift regenerates the field trace-CCZ.
- **Gluesing-Luerssen–Pllaha, arXiv:1710.09884** ✅ (❌ misattributed in scouts to "Dutta/Dastbasteh-Klappenecker") — [ring] · [no gate]. *Corrected substance:* Thm 5.4 **proves** a free ring stabilizer code's relative distance ≤ the residue-field code's (cannot outperform). Conj 5.5 is the **reverse/equality** direction (no underperform), **free codes only**. Non-free codes behave erratically (Ex 6.3 over Z_8) and are an **open loophole**. Takeaway: any ring payoff must come from the *gate*, not code parameters.
- **Gunderman, arXiv:2501.04888** ✅ — [ring] · [no gate]. Stabilizer codes over genuine composite rings Z_d **with zero divisors** (the correct substrate, Z_4/Z_8). Contains no magic-state / non-Clifford-gate result — a substrate to build on, not a solution.
- Supporting ring/character-sum tooling (⚠, several without arXiv IDs — see §7 flags): GRM over Z_{p^s} (Bhaintwal–Wasan, DOI only), Gauss/Jacobi sums over Galois rings (1603.02018, 2001.04028), Galois-ring bases (1410.0289), Witt-ring / power-sum divisibility (2210.12433, 2304.07605).

### 3.3 Magic theory for non-prime d

- **Zhu, arXiv:1504.03773** ✅ — [general] · [no gate]. No Clifford-covariant discrete Wigner function in any even prime-power dimension; 2-design representation only for odd prime powers + d=2, **d=8** (Hoggar). Kills mana/negativity for d=4,8,16; flags d=8 as structurally special.
- **Veitch–Mousavian–Gottesman–Emerson, arXiv:1307.7171** ⚠ — [field] · mana defined only for odd dimension.
- **Howard–Wallman–Veitch–Emerson, arXiv:1401.4174** and **Delfosse et al., arXiv:1610.07093** ⚠ — contextuality = Wigner-negativity only for odd local dimension; fails in even d via state-independent (Mermin-square) contextuality.
- **Delfosse et al. (rebits), arXiv:1409.5170** ⚠ — even-d quasiprobability *can* be salvaged by restricting the subtheory; template for a Z_{2^k} frame despite Zhu's no-go.
- **HW-based monotones:** stabilizer Rényi entropy (Leone et al., 2106.12587 ⚠; qudit/CV extension 2408.15161 ⚠), robustness of magic (Heinrich–Gross, 1807.10296 ⚠). Dimension-agnostic but **monotonicity proven only for prime d**.
- **Bound magic:** Zurel–Davis (2602.22336 ⚠), Prakash et al. (1905.00392 ⚠) — magic-presence ≠ distillability, a caveat that also threatens ring d=4,8.

### 3.4 Recent transversal non-Clifford constructions (the field/multi-qudit baseline)

All produce **multi-block** resources; none reaches a single-qudit T.

- **Cervia–Lamm–Liu–Murairi–Zhu, arXiv:2512.21874** ✅ — [field] · [multi-qudit]. Good triorthogonal codes over F_{2^{2m}}; transversal gate is the **trace-CCZ (−1)^{Tr(xyz)}** over GF(q), smallest q=64, e.g. [[42,14,6]]_64. Does not disclaim single-qudit magic — the multi-qudit outcome is *structural, not stated*.
- **Golowich–Guruswami, arXiv:2408.09254** ⚠ — [field] · [multi-qudit]. Asymptotically good codes, transversal inter-block CCZ over prime-power q. (Author/coauthor attribution varies across scouts — Guruswami vs "Ting-Chun Lin"; treat as Golowich–Guruswami.)
- **Nguyen, arXiv:2408.10140** ⚠ (❌ misattributed as "Wills–Hsieh–Tan" in one scout) — [qubit] · [multi-qudit]. Good binary codes with transversal CCZ.
- **He–Vaikuntanathan–Wills–Zhang, arXiv:2502.01864 / 2507.05392** ⚠ — [qubit] · [multi-qudit]. Addressable transversal CCZ on chosen logical triples; "addressable orthogonality."
- **Guémard, arXiv:2510.19809; Gasnier–Guémard, arXiv:2606.27211** ⚠ — [field/group-algebra] · [multi-qudit]. Qudit codes with transversal C^{m-1}Z / C^mZ across m blocks.
- **Wills–Hsieh–Yamasaki, arXiv:2408.07764** ⚠ — [qubit] · γ=0 constant-overhead MSD, qubit only.
- **Lin (sheaves), arXiv:2410.14631** ⚠; **Kobayashi–Zhu–Hsin, arXiv:2511.02900** ⚠ — cohomological/topological origin of transversal CCZ; reinforces that the resource is degree-m cup-product → multi-block.
- **Haruna, arXiv:2602.14499** ✅ (partial) — [qubit] · Bockstein-type obstruction to transversal diagonal-gate implementability. *Corrected:* the framework is **qubit-only, powers-of-two only** (coefficient tower Z_{2^m}→Z_{2^{m+1}}); it does **not** admit Z_{p^e} coefficients, and the literal Z_2→Z_4 Bockstein lands in H^{n+1}(·,Z_2), not Z_4. Usable as a lift test for char-2 only, not the general ring.

### 3.5 Prime-dimension single-qudit baseline (what to replicate)

- **Campbell–Anwar–Browne, arXiv:1205.3104** ⚠; **Krishna–Tillich, arXiv:1811.08461** ⚠; **Saha–Prakash, arXiv:2510.10852** ⚠ — [field] · [single-qudit]. Punctured Reed–Muller over F_p → transversal single-qudit T, γ<1 for all primes. The structure the ring route must reproduce one ring-level up.
- **Prakash et al., arXiv:2003.07164, 2403.06228** ⚠; **de Silva, arXiv:2011.00127** ⚠; **Zurel–Jana–de Silva, arXiv:2603.18560** ⚠ — prime-only single-qudit magic states / teleportation / QR-code MSD.

### 3.6 Hardware (genuinely ring Z_d)

- **Brock et al., arXiv:2409.15065** ⚠ — [ring] · [single-qudit]. GKP qutrit and ququart in one oscillator, beyond break-even; logical ops are clock-shift X_d, Z_d with primitive d-th root — a genuine cyclic Z_d ring qudit. No magic state.
- **Transmon ququart** (2304.11159, 2303.04796, 2212.04496 ⚠): arbitrary single-qudit diagonal phase via virtual-Z is hardware-trivial; the missing piece is a *distillable* encoded version. Caution: device decomposition into two qubits is the GF(4) collapse to avoid.
- **GKP obstruction:** Hastrup–Andersen (2009.05309 ⚠) — the cubic phase gate fails as a GKP non-Clifford T even at infinite squeezing → magic-state injection is the only route. **Molecular spin** single-qudit FT memory (2307.10761 ⚠) — no logical non-Clifford gate.

---

## 4. Ranked candidate attack routes

### R1 — Full-ring triorthogonality over Z/p^m targeting T_s (promise: HIGH)
- **Idea.** Define a triorthogonality analog on codes evaluated over the **full ring Z_{2^k}^m** (every point, additive character ω^x), with a divisibility/power-sum condition forcing a single physical diag(ζ^{f(x)}) (16th-root quadratic or 8th-root non-additive cubic) to be transversal and to realize a logical single-qudit T_s.
- **What must be true.** A nontrivial cubic/quadratic ring identity (the analog of Σ_x x^a = 0 mod p) must hold with simultaneously good rate (γ<1) and self-orthogonality; the surviving phase must be provably level-3 non-Clifford on the cyclic clock.
- **Why it might work.** Over Z_{2^k} squaring/cubing are **non-additive** (✅ verified), so the field collapse need not recur; full-ring quadratic Gauss sums equal (1+i)√(2^k) and **can vanish** mod m≡2 (✅ verified) — a cancellation lever the field lacks; BRG fixes the exact target gate.
- **Why it might fail.** All existing GRM dual/distance theorems are Teichmüller-based (✅) and inapplicable to full-ring evaluation — duality and distance must be rederived from scratch; the good-rate-vs-triple-cancellation tradeoff is unproven; a small-prime degeneracy could still collapse the phase.
- **Decisive test.** Construct any small full-ring code over Z_4 and check (i) self-orthogonality under the ring symplectic form, (ii) that transversal diag(ζ_16^{x²}) or diag(ζ_8^{x³}) fixes the codespace, (iii) that the induced logical gate is a non-Clifford T_s (commutator into level 2, not level 1). Even one [[n,1]]_4 example confirms existence.

### R2 — Synthesis-first computational search for transversal T_s over Z_4 / GR(4,m) (promise: HIGH)
- **Idea.** Fix T_s as the target; search small stabilizer codes over Z_4 (and GR(4,m)) for one admitting it transversally or addressably, reusing the existing qmsd discovery toolkit retargeted from field-CCZ to ring-T_s.
- **What must be true.** A computable ring-symplectic stabilizer formalism (available — Hostens–Dehaene–De Moor, Gunderman) plus a transversality checker for diag(ζ^{f(x)}).
- **Why it might work.** One example settles the central existence question without needing a closed-form ring-RM theory first; mirrors how punctured-RM was found in the prime case.
- **Why it might fail.** The search space is larger over rings; the field-derived Singleton/triorthogonality feasibility filter does not transfer, so the search may be broad or return only Clifford-degenerate gates.
- **Decisive test.** A brute/heuristic search over Z_4 codes up to moderate n returning a code with a verified non-Clifford transversal T_s — or a systematic null result that sharpens the obstruction.

### R3 — Determine the cyclic-Z_{2^k} single-qudit hierarchy (enabling theory) (promise: MEDIUM-HIGH)
- **Idea.** Classify which diagonal gates diag(ζ_{2^{k+j}}^{f(x)}) on the cyclic clock are level-3 and non-Clifford — the object de Silva–Lautsch and CGK leave open for prime powers.
- **What must be true.** A worked Clifford/Heisenberg–Weyl group over Z_{2^k} (exists: de Beaudrap 1102.3354 ⚠, Appleby quant-ph/0412001 ⚠, with the even-d 2d-doubling) and a commutator computation.
- **Why it might work.** This is a finite, concrete computation; it directly answers whether diag(ζ_8^{x³}) (which is {1,ζ_8,1,ζ_8³} on {0,1,2,3}) is genuinely level-3 or secretly semi-Clifford.
- **Why it might fail.** The cyclic constraint may force additivity/Clifford for all low-degree forms, which would itself be an important no-go.
- **Decisive test.** Compute U P U† for the candidate U and the clock-shift Paulis; check level-2 membership of the result and level-1 non-membership of the original.

### R4 — Port RCP symmetric-matrix forms to a single cyclic qudit (promise: MEDIUM)
- Reuse the only correct ring-precision dictionary, but the carrier is intrinsically multi-qubit / 2-local (✅ v2 erratum); folding to one cyclic clock via Gray/normal-basis **reintroduces the tensor factorization** the project forbids. Decisive test: whether a 1×1 symmetric form over Z_8 on one clock register stays non-Clifford and transversal — likely not without new structure.

### R5 — Char-2 magic monotone for the success metric (mandatory tooling) (promise: MEDIUM)
- mana/Wigner negativity are provably unavailable (✅ Zhu). Build/prove a HW-based monotone (SRE, robustness) over the cyclic Z_{2^k} clock group. Decisive test: prove monotonicity under ring-Clifford operations for d=4, or exhibit a counterexample. Without this, no distillation threshold can be certified. The **d=8 Hoggar-lines** exception (✅) may give a privileged witness worth separate examination.

### R6 — Gunderman beyond-integral-domain codes + bolt-on ring triorthogonality (promise: MEDIUM)
- Genuine Z_4/Z_8 stabilizer codes with zero divisors exist (✅) but carry no gate; supply the substrate for R1/R2. Risk: Frobenius/character-sum degeneracy collapses any added cubic phase.

### R7 — Trichotomy-guided universality filter (promise: MEDIUM)
- Enumerate single-ququart level-3 gates, keep those passing the BRG universality criterion (s ∤ K_d **and** s > ~1.57(d−1)), then seek MSD codes for them. Risk: a universal single-qudit gate may still demand an entangling partner for *distillation*.

### Dead ends (with reason)
- **Field GF(p^e) RM/RS/AG codes for single-qudit magic** — structurally forced to the multi-qudit trace-CCZ / C^{m-1}Z carrying only p-th roots (✅ 2512.21874). Settled negative baseline.
- **Teichmüller-evaluated GRM over rings** — cancellation is the multiplicative geometric series (✅ 1801.05114), regenerating the field collapse.
- **Gray-map / normal-basis Z_4 codes; tensoring p^e → e prime qudits** — reintroduce the multi-qubit factorization (✅ CGK standard Pauli group ≅ r p-qudits).
- **CCZ-to-T catalysis / code-switching / synthillation to localize a field-CCZ** — phase-order obstruction: Clifford + field-CCZ generate only order-p phases and catalysis is reversible (Gidney–Fowler 1812.01238 ⚠, ADCP 1403.2734 ⚠), so p^e-th-root precision cannot be manufactured.
- **mana / discrete-Wigner negativity at even d** — no Clifford-covariant Wigner function exists (✅ Zhu).
- **8th-root quadratic phase as the d=4 magic gate** — diag(ζ_8^{x²}) **is the Clifford S** (✅ verified). Memory corrected.

---

## 5. Known obstructions and no-gos (with reason)

1. **Standard-Pauli factorization** (✅ CGK 1608.06596). For p^r with the tensor/field Pauli group, the Clifford hierarchy ≅ r separate p-qudits. *Reason:* the Pauli/clock group is X(v)=⊗X^{v_i} over Z_p. *Consequence:* single-qudit prime-power magic must use the cyclic Z_{p^r} clock group, which is unclassified.
2. **Field trace-CCZ collapse** (✅ 2512.21874). GF(p^e) triorthogonality → ω^{Tr(xyz)}, ω a p-th root, multi-qudit, Clifford within each factor. *Reason:* field-linear trace carries only p-th roots; Frobenius x→x^p is additive in char p.
3. **Even-d Clifford precision** (✅ verified). S = diag(ω_{2d}^{x²}); the 8th-root quadratic on Z_4 is Clifford. *Reason:* even-d Heisenberg–Weyl requires the 2d-doubling, absorbing 2d-th-root quadratic phases into the Clifford group.
4. **No even-d Clifford-covariant Wigner function** (✅ Zhu 1504.03773). *Reason:* permutation symmetry uniquely fixes the Wigner function; the fixing fails in even prime power (except the d=8 Hoggar case). *Consequence:* mana unavailable; rebuild the witness.
5. **Free-ring-code distance bound** (✅, corrected — Gluesing-Luerssen–Pllaha 1710.09884). Free ring stabilizer codes cannot *outperform* residue-field codes (proven, Thm 5.4); equality conjectured (Conj 5.5). *Consequence:* the payoff must be the gate, not code parameters. *Loophole:* non-free codes are uncovered and behave erratically.
6. **Howard–Vala arithmetic obstruction** (✅ 1206.1598). The prime closed form needs 12^{-1}, 2^{-1}; 2 is not a unit in char 2. *Consequence:* the prime formula cannot be reused for powers of two; a new ring construction is required.
7. **Generalized Eastin–Knill** (Zeng–Cross–Chuang, PRA 78, 012353, ⚠ no arXiv): no subsystem stabilizer code has a universal transversal set on even one logical qudit, any d. *Consequence:* distillation (not free transversality) is mandatory, justifying the MSD framing.

---

## 6. Open questions / gaps in the literature

1. **Existence (central).** Does ANY stabilizer/CSS-type code over Z/p^m (or a Galois ring) admit a transversal or addressable **single-qudit T_s** with phase a genuine 2^{k}- or 2^{k+1}-th root for d=2^k — and if so, can γ<1? Existence, not just overhead, is unknown.
2. **Cyclic-clock hierarchy.** Extend the de Silva–Lautsch / CGK single-qudit diagonal hierarchy from prime d to the **cyclic ring Z_{p^e}**. Is diag(ζ_8^{x³}) or diag(ζ_16^{x²}) on Z_4 provably level-3 non-Clifford?
3. **Full-ring RM theory.** Duality and minimum-distance theorems for RM-type codes evaluated over the **full ring** Z_{2^k}^m (additive character) — none exist; all are Teichmüller-based.
4. **Ring triorthogonality identity.** Is there a full-ring power-sum / Gauss-sum cancellation Σ x^a = 0 mod 2^k (and weight-2/3 analogs) giving good rate AND the cubic cancellation? 2-adic valuations are known but not packaged as a triorthogonality condition.
5. **Char-2 monotone.** Is SRE / robustness a magic monotone under ring-Clifford operations for composite d=4,8? Does single-ring-qudit **bound magic** exist?
6. **Non-free loophole.** Does the free-code distance bound extend to non-free ring codes, or is there a genuine non-free advantage?
7. **d=8 special structure.** Does the Hoggar-lines exception (Zhu) give a privileged single-qudit magic witness at d=8 unavailable for d=4,16?
8. **BRG universality → distillation.** Does the universal single-qudit T_s admit a code-based (distillable) implementation, or does the absence reflect a deeper transversality no-go over Z/p^m?

---

## 7. Verified references

**Independently verified in the verification pass (arXiv real, substance confirmed).** These are safe to cite with the stated correction.

- arXiv:2512.20787 — Borda, Rincón, Galindo, *Quantum Universality in Composite Systems: A Trichotomy of Clifford Resources* (2025). ✅ Single-qudit ring T_s, universal, **no distillation code**. **Correction:** Clifford iff s|K_d (K_d=d odd, 2d even); universality needs s > ~1.57(d−1).
- arXiv:1608.06596 — Cui, Gottesman, Krishna, *Diagonal gates in the Clifford hierarchy*, PRA 95, 012329 (2017). ✅ p^r standard Pauli ≅ r p-qudits; does not cover cyclic Z_{p^r} clock.
- arXiv:2501.07939 — de Silva, Lautsch, *The Clifford hierarchy for one qubit or qudit* (2025). ✅ Prime-d only.
- arXiv:1206.1598 — Howard, Vala, *Qudit versions of the qubit π/8 gate*, PRA 86, 022316 (2012). ✅ Needs 12^{-1}, 2^{-1}, Z_d a field.
- arXiv:1902.04022 — Rengaswamy, Calderbank, Pfister, *Unifying the Clifford Hierarchy via Symmetric Matrices over Rings*, PRA 100, 022304 (2019). ✅ Z_{2^k} = phase ring of an n-qubit register; v2 erratum: 2-local + certain higher-locality; **not** a single 2^k-qudit theory.
- arXiv:1710.09884 — **Gluesing-Luerssen, Pllaha**, *On Quantum Stabilizer Codes derived from Local Frobenius Rings*. ✅ **Flag: scouts misattributed this to "Dutta / Dastbasteh-Klappenecker."** Thm 5.4 proves free ring distance ≤ field; Conj 5.5 (equality) is free-codes-only; non-free open.
- arXiv:1801.05114 — Andriatahiny, Ratahirinjatovo, Andrianalisefa, *Generalized Reed-Muller codes over Galois rings* (2018). ✅ Teichmüller evaluation; multiplicative-geometric-series cancellation, not a Z_{p^s} additive power-sum.
- arXiv:2602.14499 — Haruna, *Homological origin of transversal implementability of logical diagonal gates in quantum CSS codes* (2026). ✅ (partial) Bockstein obstruction real but **qubit / powers-of-two only — does NOT admit Z_{p^e} coefficients**.
- arXiv:1504.03773 — Zhu, *Permutation Symmetry Determines the Discrete Wigner Function*, PRL 116, 040501 (2016). ✅ No even-prime-power Clifford-covariant Wigner function; odd prime powers + d=2, d=8.
- arXiv:2512.21874 — Cervia, Lamm, Liu, Murairi, Zhu, *Magic State Distillation using Asymptotically Good Codes on Qudits* (2025). ✅ Trace-CCZ over F_{2^{2m}}, q=64, [[42,14,6]]_64; multi-qudit, no single-qudit disclaimer.
- arXiv:2501.04888 — Gunderman, *Beyond Integral-Domain Stabilizer Codes* (2025). ✅ Composite-ring (zero-divisor) stabilizer codes; no magic gate.

**Verified algebraic / number-theoretic facts (no arXiv; checked by direct computation).**
- Over Z_{2^k} (k≥2), x→x² and x→x³ are **non-additive**; over GF(2^e) Frobenius is additive. ✅
- Even-d Clifford phase gate S = diag(ω_{2d}^{x²}); diag(ζ_8^{x²}) on Z_4 is Clifford. ✅
- Quadratic Gauss sum mod 2^k = (1+i)√(2^k) (an 8th-root-phase object); Gauss sums over Z/mZ vanish when m≡2 mod 4; field Gauss sums have modulus √q and never vanish. ✅

**Cited by scouts but NOT independently verified here (plausible, use with caution).** Prime-baseline and tooling: 1205.3104, 1811.08461, 2510.10852, 2011.00127, 2003.07164, 2403.06228, 2603.18560, 1209.2426, 1910.09333, 2001.04887, 1709.08658, 2107.09684. Recent multi-qudit transversal: 2408.09254, 2408.10140, 2502.01864, 2507.05392, 2510.19809, 2606.27211, 2408.07764, 2410.14631, 2511.02900. Ring codes: math/0207208, 1704.06375, 1603.02018, 2001.04028, 1410.0289, 2210.12433, 2304.07605. Magic theory: 1307.7171, 1401.4174, 1409.5170, 1610.07093, 2106.12587, 2408.15161, 1807.10296, 2602.22336, 1905.00392. Stabilizer-over-ring: quant-ph/0408190, 1102.3354, quant-ph/0412001, 2209.01449. Catalysis/switching: 1812.01238, 1403.2734, 1311.0879, 1703.03860, 1606.01906, 1606.01904, 1603.03948, 2305.07720, 2407.05683, 1803.03228. Hardware: 2409.15065, 2304.11159, 2303.04796, 2212.04496, 2009.05309, 2409.05455, 2307.10761. Other: 1605.07639, 2605.04758.

**Flagged citation problems (do not cite as-is).**
- **arXiv:1710.09884** — authors are **Gluesing-Luerssen & Pllaha**, NOT Dutta/Dastbasteh-Klappenecker (scout error); and the "does not outperform" result is *proven* while the *equality* is the conjecture, free-codes-only.
- **arXiv:2408.10140** — *Good binary quantum codes with transversal CCZ gate* is by **Quynh T. Nguyen**, not "Wills–Hsieh–Tan."
- **arXiv:2408.09254** — Golowich–Guruswami (one scout wrote "Golowich, Ting-Chun Lin"); confirm coauthor before citing.
- **No arXiv ID located** (cite by journal/DOI only, or treat as unverified): Nadella–Klappenecker *Stabilizer Codes over Frobenius Rings* (ISIT 2012); Kumar–Helleseth–Calderbank *Weil sums over Galois rings* (IEEE 1995); Bhaintwal–Wasan *GRM over Z_q* (DCC 2010, DOI 10.1007/s10623-009-9315-x); Zeng–Cross–Chuang (PRA 78, 012353); Conrad expository Gauss-sum notes; Brock et al. Yale qudit-QEC (Nature 2025, s41586-025-08899-y).
- **"de Silva-Lautsch cover prime-power via tensor decomposition"** — *down-weighted*: the paper is prime-d only and merely remarks a p^t qudit can be identified with t p-qudits (the multi-qudit picture); it makes no cyclic-Z_{p^t} claim.
- **Project-memory claim "8th roots give level-3 for d=4"** — **refuted**; corrected to require 16th-root quadratic or 8th-root non-additive cubic.

---

*Bottom line: the field route is a confirmed dead end for single-qudit magic; the ring route is genuinely open with the target gate (BRG T_s) now specified and the key algebraic lever (full-ring non-additivity + vanishing ring Gauss sums) identified. The decisive next step is to establish existence — via the cyclic-clock hierarchy computation (R3) and a full-ring triorthogonality construction or computational search (R1/R2) — while rebuilding a char-2 magic monotone (R5) to make any result certifiable.*