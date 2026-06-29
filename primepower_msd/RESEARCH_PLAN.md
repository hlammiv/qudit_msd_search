<!--
Generated 2026-06-28 by a 40-agent literature+planning workflow (~1.1M tokens).
All load-bearing claims were adversarially verified; see confidence flags / Verified references.
-->

# Authoritative Plan — Single-Qudit Prime-Power Magic State Distillation over the Ring Z_{2^k}

> ## ✅ RESOLVED 2026-06-29 — SINGLE-QUDIT d=4 MSD EXISTS (see `M3_FINDINGS.md`)
> A verified **`[[9,1,≥2]]₄`** CSS code over cyclic Z₄ admits a **transversal genuinely-quadratic level-3
> gate** (antidiff(S³Z¹)) inducing the strict level-3 logical gate **diag(1,−i,1,i)** — confirmed by dense
> codespace simulation in C^(4⁹) (leakage 1.3e-15; 42 tests pass). The n≤7 obstruction below (M2/M2b/R1)
> was a **finite-size artifact**; the M3 scalable Howell kernel reached n≥9 where witnesses live.
> Remaining work is *efficiency, not existence*: low overhead (γ<1), exact ring distance/A_d, the
> distillation routine, and d=8/16/32.
>
> ## ⛔ LOW OVERHEAD — γ<1 is BARRIERED for this construction (see `M4_FINDINGS.md`)
> Minimal code is **[[8,1,2]]₄** (the [[9,1,2]] minus its free |0⟩ ancilla; dense-verified, γ=3.0); k≥2
> codes are direct sums (γ stays 3.0). **γ<1 is a hard structural barrier** over cyclic Z₄, caused by the
> **zero divisor 2**: distance is capped at 2 (weight-1 Z-logical 2·eₚ; d≥3 never seen in ~470k+
> constructions) and the separable (single-qudit) rate at n/k≥8. The open route to γ<1 is a **GF(4) field
> punctured-RM** code with a forced separable (non-CCZ) transversal gate — reaches γ<1 or yields a clean
> no-go (γ<1 single-qudit ⇒ zero-divisor ring ⇒ d=2). Publishable either way.
>
> **STATUS 2026-06-29 — Phase 0 (M0+M1+K1), M1.5, M2 COMPLETE. See `PHASE0_FINDINGS.md`, `M2_FINDINGS.md`. 42 tests pass.**
> - **K1 CLEARED** (d=4,8,16,32): single-qudit C₃\C₂ magic gates exist, exact count d(d−1). Target gate is
>   the wraparound-corrected quadratic `antidiff(S^a Z^b)` (d=4: `diag(1,ζ₁₆,i,ζ₁₆¹³)`), NOT the naive
>   `diag(ζ_{4d}^{x²})` (level ≥4) and NOT BRG's `T_s` (universal but not level-3). [corrects §1/§4 below]
> - **M1.5 DONE:** `ringlinalg.py` Howell-normal-form over Z_{2^k}, validated vs brute force.
> - **M2 RESULT (a characterized NULL):** transversal level-3 single-ququart gates DO act on cyclic-Z_4
>   CSS codes (the **non-free** loophole is essential — free codes give nothing), but in the weakly-self-dual
>   k=1 family (n≤7, 443 codes, complete 12-gate level-3 family, uniform + {0,1}-addressable) they occur
>   ONLY on **distance-1 (unencoded)** codes. **No distance-≥2 distillation code found** — the protected
>   codes reject every transversal level-3 gate. Not a no-go proof; bounded to M_X=M_Z, n≤7, d=4.
> - **NEXT (M2b/M3):** general CSS (M_X≠M_Z) [most promising]; M3 Howell-scaled distance to reach larger/
>   punctured n; the analytic ring-triorthogonality construction (R1); or a no-go proof for this family.

*Principal-investigator synthesis of the math-first (A), engineering-first (B), and skeptic-first (C) drafts, reconciled against the verified literature map. Where the drafts disagree, the resolution and its rationale are stated inline. Target: `/home/hlamm/Desktop/QC/prime_msd/primepower_msd/` (empty); reuse base: `/home/hlamm/Desktop/QC/prime_msd/qmsd/` (prime-only, 18 modules, `galois`-backed).*

---

## 1. Problem statement & definition of success

**Problem.** Construct (or prove impossible) a stabilizer/CSS code on a **monolithic cyclic qudit** of dimension d = p^e (priority d = 4, then 8, 16) whose distilled output carries a **single-system, non-Clifford, Clifford-hierarchy level-3 diagonal gate** (a Borda–Rincón–Galindo T_s over the ring Z/dZ), and a distillation routine that purifies the corresponding magic state. The carrier Heisenberg–Weyl group is the **cyclic clock group** X|j⟩=|j+1 mod d⟩, Z|j⟩=ζ_d^j|j⟩ with the even-d 2d-phase doubling — **not** the tensor/field Pauli group (CGK 1608.06596: the latter splits p^r into r p-qudits, which is the forbidden collapse).

**What counts as success — three nested certificates, all mandatory for the headline claim:**

1. **Existence (transversality).** An explicit `[[n,k]]_4` ring stabilizer code (k≥1) and a physical diagonal gate U = ⊗_i diag(ζ^{f(x_i)}) that fixes the codespace and induces a logical diagonal gate U_L.
2. **Magic (level-3).** Iterated-commutator certificate: U_L P U_L† ∈ level-2 (Clifford) for every clock-shift Pauli P, but U_L ∉ level-2. The phase precision must exceed the even-d Clifford floor: **16th-root quadratic diag(ζ_16^{x²}) or 8th-root non-additive cubic diag(ζ_8^{x³})** (the 8th-root quadratic diag(ζ_8^{x²}) is the Clifford S — the old memory is corrected here).
3. **Non-reducibility (anti-collapse).** No single-qudit Clifford basis change maps U_L's carrier (cyclic Z_4 Pauli group) onto a Z_2×Z_2 tensor Pauli group under which U_L factorizes into a 2-qubit gate. This is the certificate that separates a genuine result from "qubit MSD in GF(4) disguise."

**Full success** additionally adds: (4) a distillation step that strictly improves output fidelity (rules out *bound* magic), certified by (5) a char-2 Heisenberg–Weyl magic monotone.

**A clean no-go is also success.** A proof that every transversal low-degree diagonal level-≤3 phase on a cyclic Z_4 (and Z_8) stabilizer code is Clifford or factorizes is publishable and terminates the program early. The skeptic framing (Draft C) is adopted: existence/non-existence with a certificate is the deliverable; overhead γ<1 is a secondary prize.

**Honest payoff bound (all drafts agree).** Free ring stabilizer codes cannot beat residue-field codes on relative distance (Gluesing-Luerssen–Pllaha 1710.09884, Thm 5.4). Therefore **the novelty is the gate, not the code parameters.** Every theorem is designed to route novelty through the phase; γ<1 is pursued only via the non-free-code loophole and only after existence is settled.

---

## 2. Ranked attack routes — primary and fallback

| Rank | Route | Role |
|---|---|---|
| **PRIMARY** | **R1 — full-ring triorthogonality over Z_{2^k} targeting T_s**, executed on a computational scaffold | The theory that, if it closes, gives a code family + closed form. **Sequenced behind, and shadowed by, R3+R2** (see §3): the census de-risks the target gate, the search can short-circuit to an existence witness before the closed-form theory is finished. |
| Enabling | **R3 — cyclic-Z_{2^k} single-qudit hierarchy census** | The cheapest decisive test. Runs first. Confirms or kills the target gate with days of finite computation. |
| Co-primary | **R2 — synthesis-first computational search** for transversal T_s over Z_4 / GR(4,m) | Settles existence with one example; reuses the most qmsd code. Runs immediately after R3 clears a target. |
| Mandatory tooling (deferred) | **R5 — char-2 HW magic monotone** | Required to *certify distillation*, not to establish existence. Built only after a state to measure exists. |
| **FALLBACK** | **Pivot to d=8** (Hoggar-lines special structure, Zhu 1504.03773) or **publish the d=4 no-go** | Fires if the d=4 gate target collapses (K1) or the search returns a clean null (K2). |
| Deprioritized | R4 (RCP fold-to-one-qudit), R6 (Gunderman substrate bolt-on), R7 (universality filter) | R6 supplies substrate to R1/R2; R4/R7 low EV. |
| **Dead — do not revisit** | Field GF(p^e) RM/RS/AG; Teichmüller-evaluated GRM; Gray-map/normal-basis folding; CCZ→T catalysis; mana/Wigner at even d; 8th-root quadratic as the magic gate | Confirmed negatives in the literature map. |

**Resolution of the central draft disagreement.** Draft A wants the math program (T2–T4) funded up front; Draft C refuses to fund heavy theory before existence. **C's sequencing wins; A's theorem statements survive as the *specification* for the search and the eventual paper.** Concretely: the full-ring Reed–Muller duality/distance theory (the months-long "reinvent from scratch" item) is **not started until the existence gate M2 is positive.** Until then, R1 contributes only its *cancellation-identity spec* (T2 below), which the R2 search consumes as a filter.

---

## 3. Phase 0 — cheap feasibility gates (run FIRST, with pass/kill criteria)

These are the consensus across all three drafts. Each is days, not weeks, and each can kill the program before any heavy build.

**G0 — Substrate calibration (≈1 day).** Build the cyclic Weyl/Clifford oracle and reproduce the three verified algebraic facts:
- x↦x² is non-additive on Z_4 (the anti-collapse lever vs Frobenius).
- diag(ζ_8^{x²}) is certified **Clifford** (level ≤2) by the oracle.
- Z_4 quadratic Gauss sum Σ_{x∈Z_4} ζ_16^{x²} = (1+i)√4 = (1+i)·2.
- **PASS:** all three reproduce. **KILL/FIX:** if the oracle reports diag(ζ_8^{x²}) as level-3, the oracle is wrong — stop and fix; nothing downstream is trustworthy.

**G1 — Single-ququart hierarchy census (R3, ≈days).** Enumerate the (small, finite) single-ququart Clifford group, close it under the commutator map, and classify diag(ζ_16^{x²}) and diag(ζ_8^{x³}) by level.
- **PASS:** at least one candidate is certified level-3 non-Clifford on the cyclic clock. Record the complete list of level-3 single-ququart diagonal phases.
- **KILL (criterion K1):** if *no* low-degree diagonal phase on one cyclic ququart is level-3 non-Clifford → single-*ququart* diagonal magic does not exist. Pivot to d=8 or publish the no-go. **Do not proceed to code search.**

**G2 — Anti-collapse pre-check on the bare gate (≈1 day, runs with G1).** For each level-3 candidate from G1, run the anti-collapse certificate (§1.3) on the *bare single-qudit gate*: does a single-ququart Clifford basis change factorize it into a 2-qubit gate?
- **PASS:** at least one level-3 candidate is irreducible. **KILL/PIVOT:** if every level-3 candidate factorizes, the d=4 resource is qubit MSD in disguise even before coding — pivot to d=8.

Only after G0–G2 pass do we build ring linear algebra and run the existence experiment.

---

## 4. The mathematical program (dependency-ordered), centered on Z_{2^k}

Theorem statements follow Draft A; **funding gates** follow Draft C (pre-/post-M2). Each theorem names the obstruction it defeats.

**L0 — Cyclic clock Heisenberg–Weyl/Clifford group over Z_{2^k}** *(pre-M2).* Construct X, Z with the even-d 2d-doubling (Appleby/de Beaudrap convention). Lemma: the diagonal Clifford generator is S = diag(ζ_{2^{k+1}}^{x²}); it normalizes the clock Pauli group. *Defeats §5.1 (standard-Pauli factorization) by never leaving the cyclic group, and §5.3 by nailing the precision floor.*

**T1 — Canonical single-qudit magic gate and its level** *(pre-M2; this is G1).* Claim: T_{2^k} := diag(ζ_{2^{k+2}}^{x²}) is level-3 non-Clifford (UXU† lands in level 2, not 1); equivalently the non-additive cubic diag(ζ_{2^{k+1}}^{x³}) is level 3. For d=4: T_4 = diag(ζ_16^{x²}) = diag(1, ζ_16, i, ζ_16^9), an s=16 BRG gate (16∤8 ⇒ non-Clifford; 16 > 1.57·3 ⇒ universal with Clifford). *Method:* finite-difference calculus on the clock; finite and exact. *Defeats §5.6 (Howard–Vala needs 2⁻¹, undefined in char 2) and refutes the old "8th-root = level-3" memory.*

**T1′ — Anti-collapse / irreducibility theorem** *(pre-M2; this is G2, elevated from Draft C).* For each level-3 gate, no single-qudit Clifford conjugation maps the cyclic Z_4 Pauli group to a Z_2×Z_2 tensor Pauli group under which the gate factorizes. *This is the operationalization of CGK obstruction §5.1 and the single most important guard against self-deception (kill-criterion K3). The other two drafts underweight it; it is promoted here to a first-class, mandatory theorem.*

**T2 — Full-ring power-sum / Gauss-sum cancellation identity** *(spec pre-M2; full proof post-M2).* Characterize when Σ_{x∈Z_{2^k}} ζ_{2^{k+2}}^{a x² + b x} (and weight-2/3 cross sums) vanishes or equals a pure Clifford phase. *Lever:* full-ring squaring/cubing are non-additive (✅), and quadratic ring Gauss sums vanish for moduli ≡2 mod 4 (✅) — a cancellation the field lacks (field Gauss sums never vanish). *Defeats §5.2 (field trace-CCZ) and the Teichmüller dead-end:* evaluation is over the **full ring** at every point with additive character, so the cancellation is the genuine Z_{2^k} additive quadratic Gauss sum, not the multiplicative geometric series of Andriatahiny et al.

**T3 — Ring-triorthogonality ⇒ transversal logical T** *(post-M2 closed form; pre-M2 only as the search predicate).* For generator G over Z_{2^k}, the transversal U = T_{2^k}^{⊗n} accumulates exponent Σ_i x_i² = Σ_{a,b} u_a u_b S_{ab}, S_{ab}=Σ_i G_{a,i}G_{b,i} in Z_{2^{k+2}} with explicit 2-adic carries. **Ring-triorthogonality:** (self) each X-stabilizer row b has Σ_i b_i² ≡ 0 (mod 2^{k+2}) and pairwise Σ_i b_i g_i ≡ 0 (mod 2^{k+1}); (logical) the surviving quadratic form equals the logical T form. Because the magic phase is quadratic, the binding condition is a **pair/quadratic** ring identity with a half-precision (mod 2^{k+1}) cross term — the carry terms 2Σ_{i<j} are exactly where ring ≠ field. Carry the cubic variant (genuine triple condition) in parallel. *Defeats §5.5 by routing all novelty through the phase.*

**T4 — Existence: full-ring Reed–Muller construction + worked example** *(post-M2).* A small full-ring RM-type generator over Z_4 (low-degree polynomials evaluated at all of Z_4^m, additive character) that is ring-triorthogonal per T3, yielding a certified `[[n,1]]_4`. If the closed form resists, the R2 search supplies the witness. *Defeats §6.3 (no full-ring RM theory exists) by rederiving only the minimal duality for one example.*

**T5 — Distance / overhead** *(post-M2).* Ring symplectic min-distance; honest γ. Expect γ≥1 generically (free-code ceiling); investigate the non-free loophole (Gluesing-Luerssen–Pllaha Ex 6.3 over Z_8) as the only γ<1 route.

**T6 — Char-2 magic monotone + distillation certification** *(post-M2, R5).* HW monotone (stabilizer Rényi entropy / robustness) over the clock group; prove monotonicity under ring-Clifford ops for d=4 or exhibit a counterexample (itself a no-go). Examine the d=8 Hoggar exception as a privileged witness. *Defeats §5.4 (Zhu: no even-d Wigner) by replacing the witness.*

---

## 5. Implementation plan — modules, lift map, libraries, milestone ladder

**The one hard library gap (Draft B's load-bearing finding).** `galois==0.4.6` provides GF(p^e) but **cannot represent the ring Z_{2^k}**. Every `.null_space()`/`.left_null_space()`/`matrix_rank` call in `triorthogonal.py`, `mindist.py`, `distance.py` assumes field division and breaks over a ring with zero divisors. **All ring-linear algebra must be rebuilt on the Howell normal form** (the canonical row form over the principal-ideal ring Z_{p^e}; gives row-span, kernel, dual, membership). No clean PyPI package exists — write it.

**Resolution of the Howell-form sequencing dispute.** Draft B makes Howell form *the* foundational first deliverable; Draft C defers it. **Resolution: Howell form is NOT on the critical path for Phase 0** — G0/G1/G2 are dense complex/cyclotomic 4×4 linear algebra needing only `numpy` and exact root-of-unity arithmetic. Howell form becomes foundational the moment we touch ring *codes* (M2 onward). Build it as **M1.5**, between the census and the existence experiment, so it is never built speculatively but is ready before it blocks.

### Module map (new `primepower_msd/` ← lifted from `qmsd/`)

| New module | Lifts | Role | Phase |
|---|---|---|---|
| `ring.py` | `field.py` (`field_power_sum`, `GFp`) | Z_{2^k} scalars (plain Python int mod 2^k — not `galois`), units vs zero-divisors, 2-adic valuation, Z_{2^k}→Z_{2^{k+2}} lift with explicit carries | M0 |
| `weyl.py` | *(new)* | cyclic clock X_d, Z_d with even-d 2d-doubling; Pauli group; ring symplectic form | M0 |
| `clifford_ring.py` | *(new)* | Clifford generators (S=diag(ζ_{2d}^{x²}), Fourier H, multiplier M_a); `level_of(U)` via iterated commutator membership | M0 |
| `single_qudit_gate.py` | *(new)* | candidate diag(ζ_N^{f(x)}); `is_clifford`, `is_level3`, `certify_magic` (T1/G1), **`anticollapse_certificate`** (T1′/G2) | M0–M1 |
| `ring_power_sum.py` | `field.py`, `pnomial.py` | full-ring Gauss/power sums (T2 engine) | M1 |
| `ringlinalg.py` | replaces `galois` calls in `triorthogonal.py` | **Howell/Smith normal form**: `row_span`, `kernel`, `dual_generator`, `in_rowspan`, `rank_profile` | M1.5 |
| `ring_stabilizer.py` | *(new; Hostens–Dehaene–De Moor / Gunderman 2501.04888)* | ring symplectic stabilizer formalism; codespace projector | M2 |
| `transversal_oracle.py` | `oracle.py` | given code + U, certify U fixes codespace; extract induced logical U_L | M2 |
| `ring_codes.py` | `reedmuller.py::rm_generator`, `triorthogonal.py::build_*` | **full-ring** (NOT Teichmüller) evaluation generator over Z_d^m; CSS builder | M2 |
| `ring_triorthogonal.py` | `triorthogonal.py` (its `% p` already works for composite modulus; conditions are new) | `is_ring_triorthogonal` (mod-2^{k+2} pair + mod-2^{k+1} cross); search predicate | M2 |
| `ring_search.py`, `ring_sampling.py` | `search.py` (joblib `manhattan_sweep`/`random_search`), `sampling.py` | retarget orchestration to ring codes + transversal-T_s objective (R2). **Singleton/r_max filter does NOT transfer** | M2 |
| `ring_distance.py`, `ring_mindist.py` | `mindist.py` (MITM int64 syndrome bijection ports; "nonzero"→"unit/zero-divisor"), `weightdist.py` (MacWilliams holds — Z_{2^k} is Frobenius). `distance.py::delta_p` does **not** transfer | ring symplectic min-distance, ring A_d | M3 |
| `ring_monotone.py` | *(new, R5)* | SRE/robustness over clock group; distillation threshold | M4 |
| `ring_codes_meta.py`, `ring_distill.py` | `codes.py` (`Code` dataclass + `gamma` reusable **verbatim**, add `ring=` field), `distillation.py` (`nbar_T`/`cost`/`delta_out` formulas reuse; the (p−1)/p depolarizing constants must be re-derived for the Z_d Pauli channel) | M4 |

`asymptotics.py` reused unchanged for γ_0 once T5 fixes the rate–distance tradeoff. `structured_ad.py`, `puncture.py` adapted opportunistically.

### Milestone ladder (test-first, with acceptance gates)

- **M0 — substrate calibrated (week 1).** ACCEPT iff the G0 calibration triple passes (x² non-additive; diag(ζ_8^{x²}) certified Clifford; Gauss sum =(1+i)·2). Smallest case: 4×4 matrices.
- **M1 — single-ququart census + anti-collapse on the bare gate (week 1–2, GATE G1+G2).** ACCEPT iff the complete list of level-3 single-ququart phases is produced AND ≥1 is irreducible (T1′). KILL on K1/K3.
- **M1.5 — Howell-form correctness (week 2).** ACCEPT iff `kernel`/`dual_generator` satisfy `G @ dual.T == 0 (mod 2^k)` with matching Howell-rank dimensions on randomized Z_4/Z_8 matrices, plus a golden test vs a known Z_4 self-dual code. *Property-based randomized tests gate everything downstream — Howell bugs poison the whole search.*
- **M2 — EXISTENCE (week 2–4, GATE G2/THE decisive experiment).** ACCEPT (positive) iff one explicit Z_4 code passes all three §1 certificates (transversality + level-3 + anti-collapse), search-found acceptable. ACCEPT (negative) iff exhaustive search to n≈8 plus a structural argument yields a clean no-go. *Bound n hard to avoid a muddy timeout result; use the lenore 32-core box for n=7,8.*
- **M3 — distance + A_d + bound-magic check (week 4–6, only if M2+).** ACCEPT iff `ring_distance` returns certified d and A_d (MITM vs MacWilliams cross-checked) AND a one-shot distillation step shows δ_out<δ_in (rules out bound magic, K4).
- **M4 — monotone + overhead (months, gated).** ACCEPT iff `ring_monotone` certifies monotonicity for d=4 (or documents a counterexample) AND γ is reported honestly against the residue-field ceiling.

---

## 6. The minimal d=4 (Z_4) proof-of-concept, end to end

The smallest object that demonstrates single-ququart MSD:

1. **Carrier.** One cyclic ququart: X_4, Z_4 with the 2d=8 phase doubling (`weyl.py`).
2. **Resource gate.** T_4 = diag(ζ_16^{x²}) = diag(1, ζ_16, i, ζ_16^9). Certified by `single_qudit_gate.certify_magic`: T_4 X_4 T_4† = X_4·diag(ζ_16^{−2x+1}) ∈ level-2, while T_4 ∉ level-2 ⇒ **level 3** (T1). Cross-check the cubic diag(ζ_8^{x³}) as backup.
3. **Irreducibility.** `anticollapse_certificate(T_4)` returns no Clifford basis change factorizing T_4 into a 2-qubit gate (T1′) ⇒ genuinely single-ququart, not GF(4).
4. **Code.** A small `[[n,1]]_4` ring stabilizer/CSS code (full-ring evaluation over Z_4^m, `ring_codes.py`), satisfying ring-triorthogonality (`ring_triorthogonal.is_ring_triorthogonal`): self-orthogonal mod-2^{4} pair + mod-2^{3} cross conditions.
5. **Transversality.** `transversal_oracle` certifies U = T_4^{⊗n} fixes the codespace and extracts U_L.
6. **Logical magic.** `level_of(U_L)` re-certifies U_L as a level-3 single-ququart T_s; `anticollapse_certificate(U_L)` re-certifies irreducibility on the *logical* clock.
7. **Distillation.** One round on noisy input ququart magic states through the code; `ring_distill` + `ring_monotone` show output magic > input below a threshold (rules out bound magic).

Steps 1–3 are Phase 0 (no codes). Steps 4–6 are M2 (the existence headline). Step 7 is M3–M4. **Steps 1–6 alone constitute the publishable result even if γ≥1.**

---

## 7. Risks, kill-criteria, and the new-vs-reduction test

**Kill-criteria (any fires → stop that route):**
- **K1 (cyclic-clock collapse, highest, tested first at M1).** No single-ququart diagonal phase is level-3 non-Clifford → single-ququart diagonal magic impossible. Pivot to d=8 or publish no-go.
- **K2 (search null).** Exhaustive Z_4 search to n≈8 returns only Clifford or factorizing gates → sharpened no-go; pivot to d=8/d=9. *Escalate n only if M1 guarantees a viable target gate exists.*
- **K3 (anti-collapse failure = NOT new).** Every candidate logical gate is Clifford-equivalent to a 2-qubit CCZ/CZ → reduction to known qubit MSD. **This is the most likely way to fool ourselves; T1′ is mandatory before any claim.**
- **K4 (bound magic).** A non-Clifford state exists but no Clifford circuit distills it (Zurel–Davis, Prakash) → keep the gate result, drop the MSD claim.
- **K5 (free-code ceiling, not fatal).** γ<1 unreachable on free codes (Gluesing-Luerssen–Pllaha) → demote to "existence + suppression," chase the non-free loophole as stretch. Do not promise γ<1 from free codes.

**Soft risks:** Howell-form bugs (mitigate: M1.5 golden + property tests before any code theory); 2d-doubling convention errors (caught by M0 calibration); MacWilliams Frobenius-character subtlety (validate by brute-force weight enumeration on a tiny code before trusting A_d); carry-arithmetic in the Z_{2^k}→Z_{2^{k+2}} lift (fix the canonical {0,…,2^k−1} lift and verify every phase claim by brute-force complex evaluation before trusting symbolic mod arithmetic); search blow-up (hard n cap + lenore 32-core box).

**How to tell a NEW result from a reduction (the decision rule).** A result is new **iff** it passes all three: (a) the non-Clifford resource lives on the **cyclic Z_4 clock** Pauli group; (b) its precision is **strictly above** the 8th-root quadratic Clifford floor (16th-root quadratic or 8th-root non-additive cubic); (c) it passes the **anti-collapse certificate** — no Clifford basis change factorizes it to a 2-qubit gate. Failing (c) is the GF(4) collapse (qubit MSD in disguise). Carrying only p-th-root phases is the field trace-CCZ baseline. Both are NOT new.

---

## 8. Recommended first concrete step (this week)

**Build `primepower_msd/{ring.py, weyl.py, clifford_ring.py, single_qudit_gate.py}` and run M0 then M1** — the cyclic-clock substrate, its calibration triple, and the R3 single-ququart hierarchy census including the anti-collapse certificate.

All three drafts converge on this, and it is correct: it is finite exact 4×4 cyclotomic linear algebra, needs **no** ring code theory and **no** Howell form, reuses none of the fragile `galois` paths, and within days returns a decisive verdict. Concretely:

1. Implement X_4, Z_4 with the 2d=8 doubling; Clifford generators S=diag(ζ_8^{x²}), Fourier H, multiplier M_a; and `level_of(U)` testing whether iterated commutators land in the Pauli group (level 1) / Clifford group (level 2) / neither, via exact 8th/16th-root arithmetic.
2. **Gate with the calibration test:** `level_of(diag(ζ_8^{x²}))` must return Clifford (validates the oracle and refutes the old memory).
3. Run `certify_magic(diag(ζ_16^{x²}))` and `certify_magic(diag(ζ_8^{x³}))`, then `anticollapse_certificate` on each level-3 hit.

**Branch decision:** if a certified, irreducible level-3 single-ququart gate exists → green-light M1.5 (Howell form) and the M2 existence experiment on the now-confirmed target phase. If both candidates collapse to Clifford or factorize → fire K1/K3, pivot to d=8 (Hoggar special structure), and begin drafting the no-go.

**Files.** New work: `/home/hlamm/Desktop/QC/prime_msd/primepower_msd/` (empty). Reuse for style/skeleton (not ring math): `/home/hlamm/Desktop/QC/prime_msd/qmsd/triorthogonal.py` (the `% p` condition pattern already works for composite modulus; the `galois` puncture/shorten/dual do not), `qmsd/oracle.py` (transversality-checker skeleton), `qmsd/mindist.py` (MITM syndrome bijection to lift), `qmsd/weightdist.py` (MacWilliams to re-derive over the Frobenius ring), `qmsd/codes.py` (`Code` dataclass + `gamma`, reusable verbatim). Ring arithmetic is hand-rolled on `numpy` int mod 2^k — `galois==0.4.6` covers GF(p^e), not Z_{2^k}.