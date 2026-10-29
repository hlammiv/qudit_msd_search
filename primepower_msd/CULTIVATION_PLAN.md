# Ququart magic-state CULTIVATION — briefing, plan, and TODOs

*Handoff document for a fresh session. Self-contained: read this + the linked findings docs and you can
start work immediately. Created 2026-06-30 after the rate-vs-transversality wall was confirmed (#2, #4).*

---

## 0. How to use this document

This is the chosen path (idea **#1**) toward the north star: **a low-overhead, genuine ququart (d=4,
ring Z₄) magic-state-distillation protocol.** Sections 1–3 are context (why cultivation, what it is, how
it maps to ququarts). Section 4 lists reusable ASSETS (verified codes/gates/tools with file paths).
Section 5 is the PLAN (phases). Section 6 is the ordered TODO list — start there once you've read 1–5.
Section 7 = open questions/risks. Section 8 = hard constraints (resource limits, genuineness cert).

Memory index: `lowoverhead-ququart-msd-ideas.md` (the 6-idea menu), `primepower-msd-program.md` (full
program history). Findings: `PHASE1`–`PHASE5_FINDINGS.md` in this directory.

---

## 1. Why cultivation (the context in one screen)

The program's goal is single-/few-ququart MSD with low overhead, on a **genuine** monolithic Z₄ ring
qudit (clock-shift Pauli, X order 4) — NOT the GF(4)/2-qubit "collapse" (which is what most naive
constructions and Henry's own arXiv:2512.21874 produce). "Genuine" is certified by the **anti-collapse
cert**: a logical Pauli pair with ODD symplectic commutator (= order-4, true ququart), not order-2.

**What is settled (don't re-derive):**
- **γ<1 genuine d=3 ring CODES exist & are verified:** `[[10,4,3]]₄` (γ=0.834), `[[12,5,3]]₄` (γ=0.797).
- **A genuine ring-CCZ magic gate exists** and these codes *host* one (Phase 5, dense-validated end-to-end).
- **BUT the hosting gate is irreducibly NON-TRANSVERSAL** (`[[10,4,3]]`: 0 pure-CCZ hosters, sparsest
  ~43 CCZ / qudit-degree 18). Confirmed with a corrected hosting solver (#2) — not a solver artifact.
- **The transversal route is blocked by Z₄ arithmetic:** transversal CCZ needs a triorthogonal base, and
  over Z₄ triorthogonality forces **d=1** — true for random *and* structured (Reed–Muller) codes (#4).
- **NET: a *transversal* fault-tolerant γ<1 ring-CCZ ququart code very likely does NOT exist** (the
  rate-vs-transversality wall, fundamental to the Z₄ zero-divisor 2).

**Why cultivation:** it does not need a transversal gate. Fault tolerance comes from the *cultivation
procedure* (prepare magic at small distance → check/post-select → grow distance), so the wall is
irrelevant. This is the path that reaches the north star regardless. It is also where qubit MSD got its
current lowest overhead (Gidney–Shutty–Jones, "Magic state cultivation," arXiv:2409.17595, 2024).

---

## 2. What magic-state cultivation IS (qubit version)

Classic MSD: encode many noisy magic states in a distillation code, apply a *transversal* non-Clifford
gate, measure stabilizers, output fewer better ones. Overhead exponent γ = log(n/k)/log(d).

Cultivation flips the order:
1. **Prepare** the magic state on a SMALL code (small distance) by a cheap, NON-fault-tolerant but
   POST-SELECTED procedure (inject the gate where it's cheap).
2. **Check / verify**: measure stabilizers (and/or a redundant check); on a bad syndrome, REJECT
   (post-selection). Acceptance is high when physical error is low.
3. **Grow the distance** (code switching / lattice surgery / morphing) so the now-good state is protected
   — paying the expensive distance ONLY AFTER the magic is in and verified.

Overhead metric is **spacetime volume** (not single-shot γ). The win: you never apply the expensive gate
at high distance; post-selection cleans the cheap small-distance preparation; distance-growth is Clifford.

---

## 3. How it maps to ququarts (Z₄) — and what's new

| ingredient | qubit | ququart (Z₄) status |
|---|---|---|
| small code hosting the magic gate | Steane/color code + transversal T | **HAVE IT:** `[[8,2,2]]₄` hosts the genuine ω=i carrier at d=2 (Phase 1); genuine ring-CCZ on small codes (Phase 5) |
| the magic gate / state | T (or CCZ) | **HAVE IT (verified):** ω=i carrier `diag(ζ₁₆^{0,1,4,13})` and `CCZ_ring=diag(i^{xyz})` (level-3, order-4) |
| post-selection check | stabilizer meas. | standard; needs Z₄ syndrome extraction (straightforward) |
| **distance growing** | color-code morphing / lattice surgery | **THE NEW PIECE — undeveloped for Z₄.** This is the core research task. |

**Key realization:** the Phase-1–3 obstruction ("diagonal ω=i carrier hosting ⟺ d≤2") does NOT block
cultivation — we deliberately prepare the magic at d=2 (where `[[8,2,2]]₄` hosts it) and grow d afterward.
Cultivation *uses* the small-d hosting we already have.

**Target magic state (to decide in Phase 0):** simplest is the single-ququart genuine ring-T (ω=i) state;
alternative is the 3-ququart ring-CCZ resource state (our Phase-5 result). Single-ququart T is the cleanest
to inject & cultivate; pick it unless catalysis (CCZ→T over Z₄, order-8→4 compatible) argues otherwise.

---

## 4. ASSETS (verified — reuse, don't rebuild)

**Codes (in `primepower_msd/` unless noted):**
- `witness_8_2_2_hosts_carrier.npy` — `[[8,2,2]]₄`, k=2, d=2, GENUINE, **hosts the diagonal ω=i carrier**
  (dense-verified). The natural cultivation seed.
- `code_10_4_3_gamma0834.npy` — `[[10,4,3]]₄`, γ=0.834, genuine d=3 (storage/output target).
- `code_12_5_3_gamma0797.npy` — `[[12,5,3]]₄`, γ=0.797, genuine d=3.
- `M2_code_12_2_3.npy` / `code_12_2_3_d3_k2.npy` — `[[12,2,3]]₄` non-CSS d=3.
- `scratchpad/genuine_*.npy`, `scratchpad/VG_*.npy` — assorted verified witnesses.

**Gates (verified genuine level-3):**
- ω=i single-ququart carrier: `diag(ζ₁₆^q)`, q=[0,1,4,13] = `diag(1,ζ₁₆,i,ζ₁₆¹³)`; 3rd finite-difference
  = 4 mod 16 → genuine (ω=i). (ω=−1 = "soft collapse", forbidden.)
- `CCZ_ring = diag(i^{x_i x_j x_k})` — multi-ququart, level-3, **order-4** (verified `scratchpad/dexp_gate.py`).

**Tools (in `primepower_msd/` and `scratchpad/`):**
- `noncss_transversal.py` — dense codespace projector + transversal-gate test (signed stabilizers; n≤6
  dense). Ground truth for small codes.
- `ringlinalg.py` — Howell normal form, kernels, module ops over Z₄ (scalable, no dⁿ).
- `clifford_ring.py` — single-qudit Clifford-hierarchy level oracle (`level_of`, `is_clifford`).
- `scratchpad/dexp_gate.py` — MULTI-qudit hierarchy oracle (`level_of_m`) — checks level + order of
  multi-ququart gates (distinguishes genuine order-4 from collapsed order-2).
- `scratchpad/cubic_solve2.py` — CORRECTED hosting solver (Φ const on logical-rep cosets; the strict
  `cubic_solve.py` is over-restrictive — use the v2).
- `scratchpad/genuine_v3.py` — induced-logical genuineness (odd 3rd finite-difference; complete,
  basis-independent). `scratchpad/logical.py` — proper logical symplectic basis.
- Anti-collapse cert: odd logical commutator — see `p2_engine.py` / `logical.py`.
- GR(4,2) machinery (repo root): `gr_engine.py`, `gr_evrm.py`, `gr_kd.py`, `gr_graphhunt.py`, `gr_gamma.py`.
- Remote compute: **lenore** (`ssh -p 60022 lenore_remote`, 32 cores / 125 GB, venv at
  `~/cap_validate_work/venv`); project synced to the SAME path `/home/hlamm/Desktop/QC/prime_msd`.

---

## 5. THE PLAN (phases)

**Phase 0 — Literature + target lock (cheap, do first).**
Read cultivation papers; lock the target magic state and the overhead metric (spacetime volume +
acceptance rate). Decide single-ququart ring-T vs 3-ququart ring-CCZ.

**Phase 1 — Cheap post-selected preparation on a small code.**
Take `[[8,2,2]]₄` (d=2, hosts ω=i). Build the prep circuit: encode logical |+⟩-analog, apply the ω=i
carrier (transversal at d=2, which IS allowed there), measure stabilizers. Characterize, by dense sim
(n=8, 4⁸ — borderline; do on lenore), the **acceptance rate** and **post-selected logical fidelity** vs a
single-qudit depolarizing rate p. This quantifies the "cheap magic" starting point.

**Phase 2 — The verification / check that drives post-selection.**
Design the error-detecting measurement (Z₄ syndrome extraction; possibly a redundant "double-check" like
qubit cultivation's morphing check). Compute post-selected logical error vs p and the acceptance cost.

**Phase 3 — DISTANCE GROWING over Z₄ (the make-or-break novelty).**
Develop a Clifford procedure that maps the small-d magic state into a larger-d Z₄ code preserving the
logical state. Candidate mechanisms to survey & prototype:
  (a) **Z₄ code switching** between two CSS-over-Z₄ codes sharing the logical (gauge-fixing).
  (b) **Z₄ lattice surgery / merging** of the small code into a patch of a larger Z₄ topological/HGP code.
  (c) **Concatenation / grafting**: embed `[[8,2,2]]₄` logical into `[[10,4,3]]₄`/`[[12,5,3]]₄` (our γ<1
     codes) as the storage layer, growing distance from 2→3+ by a Clifford encoder.
  (d) **Morphing** (color-code-style) if a Z₄ color/topological code is available.
Start with (c) (we have both codes) or (a); these reuse existing assets. (b)/(d) need new Z₄ topological
constructions (bigger lift, parallels the untested HGP idea).

**Phase 4 — Overhead analysis.**
Spacetime volume, acceptance-weighted cost, output infidelity. Compare to (i) Henry's prime-dimension
MSD, (ii) qubit cultivation, (iii) the γ-based bound. Show the spacetime-volume win.

**Phase 5 — End-to-end + verification.**
Assemble the smallest full cultivation cycle; verify by dense/stabilizer sim where feasible (lenore);
confirm the OUTPUT state passes the anti-collapse (genuine order-4) cert.

---

## 6. TODOs (ordered — start here)

1. **[research] Pull cultivation literature.** Gidney–Shutty–Jones arXiv:2409.17595 ("Magic state
   cultivation") + its lineage (Li 2015 post-selected injection; Chamberland–Noh; color-code morphing).
   Extract: the exact cult. cycle, the check structure, how distance is grown, the overhead accounting.
   (Use the `deep-research` skill or WebSearch.) Output: a 1-page mechanics summary appended here.
2. **[decision] Lock the target magic state** (single-ququart ring-T vs 3-ququart ring-CCZ) and the
   overhead metric. Record rationale.
3. **[build] Phase-1 prep+post-select on `[[8,2,2]]₄`.** Dense sim on lenore: acceptance rate &
   post-selected logical fidelity vs depolarizing p. Reuse `noncss_transversal.py`. (RAM: 4⁸ dense is
   ~heavy — keep on lenore, never local; see §8.)
4. **[build] Phase-2 check.** Z₄ syndrome extraction + post-selection logic; logical error vs p.
5. **[research+build] Phase-3 distance-growing prototype.** Start with mechanism (c): a Clifford encoder
   mapping the `[[8,2,2]]₄` logical magic into `[[10,4,3]]₄`/`[[12,5,3]]₄`. Prove (symbolically, via the
   stabilizer/symplectic tools) the encoder is Clifford and preserves the logical state. This is the
   research crux — budget the most time here.
6. **[verify] Anti-collapse the OUTPUT.** Confirm the cultivated, grown state is still genuine order-4
   (odd logical commutator) — the magic must not collapse during growth.
7. **[analysis] Phase-4 overhead.** Spacetime volume + acceptance; compare to qubit cultivation & the
   γ-bound. Decide if the win is real.
8. **[write] Update `PHASE6_FINDINGS.md` + memory** at each milestone (per the program's documentation
   discipline).

---

## 7. Open questions / risks

- **Z₄ distance-growing is undeveloped — the main research risk.** No Z₄ code switching / lattice surgery
  / morphing in the literature (as far as the program knows). Phase 3 is genuinely novel; it may need a
  new Z₄ topological or CSS-pair construction. Mechanism (c) (graft into our γ<1 codes) is the lowest-risk
  starting point because both codes are in hand.
- **Does the genuine (anti-collapse) property survive growth?** Clifford distance-growing should preserve
  the logical magic, but the order-4/odd-commutator structure must be checked end-to-end (TODO 6).
- **Is the overhead actually low?** Acceptance rate must stay high; the spacetime-volume win depends on it.
- **Single-ququart vs CCZ:** if single-ququart ring-T can't be cleanly cultivated, the 3-ququart ring-CCZ
  (Phase-5 result) is the fallback target — but injection/catalysis is more involved.
- **Connection to the untested HGP idea:** mechanisms (b)/(d) overlap with the homological/hypergraph-
  product-over-Z₄ idea (the one structured family never tried). If Phase 3 pushes toward topological Z₄
  codes, revisit HGP.

---

## 8. Hard constraints (read before running anything)

- **RESOURCE CAUTION (a local OOM crash happened in this program).** Local box = 15 GB, shared with the
  user's apps. **NEVER run dense 4ⁿ simulations locally for n≥6** (the RAM spiker). Do all dense work on
  **lenore** (125 GB) via `ssh -p 60022 lenore_remote ... ~/cap_validate_work/venv/bin/python`. Always:
  put `timeout N` INSIDE the ssh command (not just ssh-client side, else the remote job orphans);
  `ulimit -v`; `OMP_NUM_THREADS` small; check `free -h` first. See memory `resource-caution.md`.
- **GENUINENESS IS MANDATORY.** Every magic claim must pass the anti-collapse cert (odd logical
  commutator / order-4). The ω=−1 ("soft") and GF(4) ("field collapse") versions are FAILURES, not magic.
- **Polynomial > dense.** Prefer the symplectic/Howell tools (`ringlinalg`, `logical`, `cubic_solve2`,
  `genuine_v3`) over dense simulation; only validate small cases densely.
- Don't trust workflow/agent claims unverified — independently re-derive load-bearing results (this
  program caught several false positives that way).

---

---

## 9. PHASE 0 — COMPLETED 2026-06-30 (literature pulled; target & metric locked)

### A. How qubit cultivation actually works (Gidney–Shutty–Jones, arXiv:2409.17595)
Three stages on a **2D color/surface-code substrate** with matching decoders:

1. **INJECTION** (d=3 color code, *fault* distance 1). A unitary injection circuit prepares the color-code
   |T⟩ with **the T gate applied in the MIDDLE of the circuit so a Z stabilizer protects it** (catches
   X/Y errors on the gate itself). The state lives in a distance-3 code but has fault distance 1.

2. **CULTIVATION** (grow the *fault* distance by checks). Repeated **"double-checks"** = a forward check +
   its time-reverse. Each: apply T† to data (turns the T-check into an X-parity measurement), entangle
   each data qubit with a |+⟩ ancilla via CNOTs, collapse ancillas along a spanning tree, measure X
   parity, then **time-reverse the whole thing** (the reverse restores helpers to |+⟩ → many flag checks).
   **FULL post-selection** (discard on ANY detection event) ⇒ retry cost grows **exponentially** with gate
   count. Code grows between checks (Bell-pair prep); ~3 stabilizer cycles between checks.

3. **ESCAPE** (the state is now "too good for its code"; one un-postselected round would damage it).
   **GRAFTING**: transform the small color code into a large *matchable* surface code preserving fidelity.
   Initialize a large surface-code region alongside; measure the combined "grafted" stabilizers
   (color + surface) for O(d) rounds; then drop the color-code stabilizers → fully matchable surface code.
   Post-selection is NOT exhaustive (circuit too big): "**desaturation**" (post-select only specific
   color-region detection events so residual errors have ≤2 detections → a matching decoder works) +
   a per-shot "**complementary-gap**" confidence threshold (discard low-confidence shots; tunable
   cost/fidelity trade-off).

**Overhead metric = spacetime volume** = Σ_shots (active qubits × rounds), survival-weighted (NOT γ).
**Headline:** d=5 cultivation at 10⁻³ depolarizing noise → **2×10⁻⁹** logical error, ~99% discard
(~100 attempts/kept shot); an order of magnitude fewer qubit·rounds than prior MSD. (At 5×10⁻⁴ → 4×10⁻¹¹.)

### B. Ququart (Z₄) mapping — what we have vs what's new

| stage | qubit mechanism | ququart (Z₄) adaptation | status |
|---|---|---|---|
| **substrate** | 2D color/surface code + matching | **Z₄ qudit surface/color code** (qudit toric codes over Z_d EXIST) | NEW but a known family — build/adopt it |
| injection | unitary circuit, T in middle, Z-stab protected | inject genuine ω=i ring-T `diag(ζ₁₆^{0,1,4,13})` in the middle, protected by a Z₄ stabilizer | **HAVE the gate**; need the Z₄ injection circuit |
| cultivation | double-checks (T†→X parity + flags), full postselect | Z₄ double-checks: ring-T†→Z₄-parity meas., controlled-shift (Z₄ "CNOT") entangling, flag checks, time-reverse | NEW: design the Z₄ double-check |
| **escape** | grafting color→surface; desaturation + gap | graft small Z₄ code → larger Z₄ surface code; Z₄ matching/decoding | **NEW — THE CRUX** (= Phase 3) |
| output | |T⟩ at high d | genuine ring-T at high d; **anti-collapse cert MANDATORY on output** | verify order-4 |

**KEY:** the Phase-1–3 obstruction (diagonal ω=i can't be *hosted* at d≥3) is **IRRELEVANT** here —
injection applies the gate at fault-distance 1 and cultivation grows the *fault* distance via checks; the
gate is never hosted transversally on the big code.

**KEY UPDATE to Phase 3 (§5):** the most faithful path uses a **Z₄ surface/color code substrate** (qudit
topological codes exist) + a **Z₄ grafting / lattice-surgery** primitive — MORE concrete than the original
"graft into `[[10,4,3]]`" suggestion. The γ<1 codes remain a possible escape *target*, but a 2D Z₄ surface
code is the natural one (matchable, scalable distance, standard decoders).

### C. TARGET & METRIC — LOCKED
- **Target magic state: the single-ququart GENUINE ring-T (ω=i) state.** Matches the |T⟩ cultivation
  paradigm and is cleanest to inject. Fallback: the 3-ququart ring-CCZ (Phase-5 result); ring-CCZ→ring-T
  catalysis (order 8→4 compatible) is a bridge if single-ququart T proves awkward.
- **Overhead metric: spacetime volume (active-qudits × rounds, survival-weighted) + acceptance/discard
  rate + output logical error vs physical depolarizing rate p.** NOT γ. Report the cost-vs-fidelity curve.

### D. References for the new session
- **Gidney, Shutty, Jones, "Magic state cultivation: growing T states as cheap as CNOT gates,"
  arXiv:2409.17595 (2024)** — THE paper. Code: github.com/Strilanc/magic-state-cultivation.
- "Efficient magic state cultivation with lattice surgery," arXiv:2510.24615 — escape via lattice surgery
  (alternative to grafting; PDF didn't text-extract — fetch the HTML/abstract).
- "Cultivating T states on the surface code with only two-qubit gates," arXiv:2509.05232.
- "Simulating magic state cultivation with few Clifford terms," arXiv:2509.08658.
- "Efficient Magic State Cultivation," PRX Quantum 10.1103/9kys-3whh.
- **Z₄ SUBSTRATE (not yet pulled — first Phase-1 sub-task):** qudit surface/toric codes over Z_d, qudit
  color codes (search "qudit topological stabilizer code Z_d surface code", Watson/Campbell, Bombín).

### E. Revised immediate next steps (supersede §6 TODOs 1–2, which are now DONE)
- **3a. [research]** Pull/choose the **Z₄ surface or color code** substrate (qudit toric code over Z₄).
  This is the prerequisite for injection; do it before 3b.
- **3b. [build]** Phase-1 injection circuit for the ω=i ring-T on a small (d=3) Z₄ code, post-selected;
  dense sim on **lenore** (never local — §8). Characterize acceptance + post-selected fidelity vs p.
- **3c. [research]** Phase-3 escape = the Z₄ **grafting / lattice-surgery** primitive — the make-or-break.
  Start from the grafting picture above (small code + large Z₄ surface code, measure combined stabilizers
  O(d) rounds, drop small-code stabilizers); or adapt arXiv:2510.24615's lattice-surgery escape.
Then continue with §6 TODOs 4 (check) → 6 (anti-collapse output) → 7 (overhead).

---

*End of briefing. A fresh session: read §1–5 and §9, then execute §9.E (substrate → injection → escape).
The crux is the ESCAPE stage = Z₄ distance-growing (grafting/lattice surgery); injection & checks reuse
the verified ω=i gate and the small-code assets.*

---

## 10. 2026-10-20 — target corrected; the T-CHECK is transversal on Z_{2^k} colour codes (verified)

**Target correction (supersedes §3/§9.C).** Since 2026-10-19 we know (UNIVERSALITY_CHECK.md) that ring-CCZ and
every mod-4 phase are non-universalizing, and the `[[8,2,2]]₄` carrier's *logical* gate is non-escaping. So the
cultivation target is the single-qudit **escaping** T-like state |T⟩ = T|+⟩ with T a strict level-3 diagonal
that escapes all Pauli sectors (d=4: q = [0,1,4,13] mod 16; general 2^k: `cultivation.escaping_carrier(d)`).
`[[8,2,2]]₄` is no longer the seed.

**Check-operator theorem (verified exhaustively, d = 4, 8, 16).** For every strict level-3 diagonal T on Z_d,
T X T† = X·diag(ζ^g) with g a Clifford table c + 2(a x² + 2bx), and
**T escapes ⇔ a is odd.** (48/448/3840 strict level-3 tables at d = 4/8/16: every escaper has a odd and every
non-escaper has a even.) So the ququart analogue of H_XY ∝ SX is the Clifford **X·S_{2d}^{odd}**, and its
+1 eigenstate is |T⟩.

**Substrate (verified).** Signed triangular 6.6.6 colour code over Z_d: X rows = 0/1 plaquettes, Z rows =
plaquettes × τ (τ = ±1 honeycomb bipartition). This gives [[ (3D²+1)/4, 1, D ]]_d (ring distance certified D = 3, 5
for d = 4, 8), with X^{⊗n} = X_L. The transversal product **∏ᵢ T^{τᵢ} X T^{−τᵢ}** (T on one sublattice, T† on
the other, exactly the qubit S/S† pattern) is hosted exactly (symbolic, CSS support mode: D ≤ 9 at d = 4, 8, 16)
and acts as **T_L X_L T_L†** with the *same* carrier. A dense state-vector check on n = 7, d = 4 gives zero
leakage, logical g = [1,3,9,3] (= physical), and T_L|+_L⟩ as a +1 eigenvector. The argument is purely
algebraic: balanced plaquettes and ± overlaps cancel, and the logical exponent is Στ = ±1. It should therefore
hold for any qudit dimension, including odd prime powers, but only 2^k has been verified numerically.
Code: `cultivation.py`; tests: `tests/test_cultivation.py`.

**Consequence.** GSJ's cultivation check transfers verbatim to every 2^k. The check is T† on data, a
controlled-SUM from a Z_d cat state, then T. That needs physical T (a native single-qudit level-3 gate), CSUM, the
Fourier gate, and Z_d cat states. The check has d outcomes (eigenvalues ω^j), so it detects every logical Z^j error,
not just a sign flip. No transversal non-Clifford gate on a large code is ever needed: the
rate-vs-transversality wall and the γ ≈ 2 MSD saturation are both irrelevant here.

**Next steps (replace §9.E):**
1. Z_d stabilizer / Pauli-frame circuit simulator (qudit Gottesman–Knill: CSUM, F, S_{2d}, Pauli noise),
   plus a Clifford proxy for T as in GSJ (replace T by a Clifford whose check is X·S^{odd}; the proxy state is a
   stabilizer state). Validate the proxy against exact state vectors at n = 7.
2. Injection (T in the middle, protected by a Z stabilizer) and the double-check on D = 3, then growth to D = 5.
   Measure the discard rate and post-selected logical error vs p, and compare with qubit cultivation at equal p.
3. Escape: graft the Z_d colour code into a Z_d surface code (both CSS, so the operation is Clifford). Decoding
   uses Z_d matching (qudit surface-code charges).
4. Overhead: spacetime volume vs output error; compare with the γ-based MSD, where the best certified option is
   now γ_T = 1.995 (Sidon-82 [[8110,82,10]]₄, d_Z certified 2026-10-20; HAILMARY_FINDINGS §3b).

## 11. Stage-1 simulations (2026-10-20)

**Tools.** `zdsim.py`: vectorised Z_d Pauli-frame simulator (R, F, F†, S_{2d}^a, CSUM, M; DEP1/DEP2/X/Z noise;
detectors; fault injection; determinism check), dense-validated for d = 2, 3, 4, 8. `cult_circuit.ColorCult`:
unitary encoder → [stabilizer round, 2 T-checks (n-qudit cat, flags)] × rounds → final stabilizer round, with all
detectors post-selected and ideal min-weight decoding at the end. T is replaced by the Clifford proxy S^σ (GSJ);
detectors are deterministic iff inj = σ. Sweep: `cult_sweep.py`; data in `runs/cultivation/`.

**Single faults.** No undetected harmful single fault for D = 3, 5 and d = 2, 4 (up to 194k faults each).

**Monte Carlo, D = 3** (20 qudits, 100 noise locations, rounds = 2, 10⁷ shots per point, uniform depolarizing p):

| d | p = 5e-4: accept / p_L | 1e-3 | 2e-3 |
|---|---|---|---|
| 2 | 0.865 / 5.8e-6 (50) | 0.748 / 1.9e-5 (144) | 0.560 / 7.0e-5 (392) |
| 4 | 0.853 / 5.7e-6 (49) | 0.727 / 2.3e-5 (170) | 0.528 / 8.7e-5 (460) |
| 8 | 0.848 / 4.0e-6 (34) | 0.719 / 2.1e-5 (154) | 0.516 / 7.7e-5 (395) |

p_L ∝ p² (fault distance 2 for this naive stage-1 circuit). Ququarts and quoctits match qubits within statistics;
acceptance is 1–4 points lower (more non-trivial Pauli outcomes per location).

**Monte Carlo, D = 5, direct encoding** (56 qudits, 292 noise locations, 10⁶ shots): accept 0.65 / 0.42 / 0.18
(d = 2), 0.62 / 0.39 / 0.15 (d = 4), 0.61 / 0.37 / 0.14 (d = 8) at p = 5e-4 / 1e-3 / 2e-3. There were 0 logical
failures at p ≤ 1e-3 (p_L ≲ 8e-6 at 95%) and 1–3 failures at 2e-3. Direct D = 5 encoding discards too much, so
GSJ-style growth from D = 3 is needed.

**Proxy validation (exact T vs Clifford proxy, dense, D = 3, n = 7).** Data-level model: an ideal measurement of
C = A†X^{⊗n}A, A = ∏ G^{−τᵢ}, with Pauli faults in 4 slots per data qudit (before A, after A, after the controlled
X, after A†), then post-selection on the check and on the trivial syndrome. Scripts: `scripts/cult_proxy/`; data:
`runs/cultivation/proxy_validate_d2.json`, `proxy_slots_d4.json`.

The results:
- **Low fault weights.** Single faults (84 / 420) and all fault pairs (d = 2) are undetected with weight 0 in both
  models.
- **Logical errors before the check.** The exact check accepts a logical X error with probability 1/2 and *projects
  it back onto |T_L⟩*. The proxy rejects it. Both give zero undetected error; the proxy is pessimistic in
  acceptance.
- **Weight-3 faults.** Only faults on logical supports matter: by cleaning, any other weight-3 fault is trivial on
  the code. Summing u over all Paulis and slot patterns gives:

| d | Σu exact T | Σu proxy | Σaccept exact | Σaccept proxy |
|---|---|---|---|---|
| 2 | 448 | 448 | 728 | 896 |
| 4 | 1704 | 1984 | 2280 | 3200 |

The distribution differs. The exact check is worse when two of the three faults lie after the controlled X
(e.g. d = 4 pattern (0,2,2): 36 vs 0), and better when exactly one does (102 vs 174). Under uniform depolarizing
data noise the totals agree exactly at d = 2 and the proxy over-estimates by 16% at d = 4. So the Clifford proxy is
a sound (conservative) estimate for qudits too. Caveats: data-level only (cat/ancilla faults are Clifford and
already exact in the frame sim), a single check, uniform slot weights.

**Next:** T-in-the-middle injection; D = 3 → 5 growth (code deformation, which is not MSD-style code switching,
but flagged against the user's "avoid code switching" preference); grafting into a Z_d surface code; spacetime
overhead.

## 12. D = 3 → 5 growth by code deformation (2026-10-20, `cult_growth.py`)
The D=3 colour code sits in a corner of the D=5 lattice.  New qudits are prepared in a |0⟩/|+⟩ split chosen so that
the D=5 logical representatives restrict correctly (2280 of 4096 splits feasible; default zset [17,18]).  The D=5
plaquettes are then measured.  Deterministic outcomes become detectors.  The 10 random ones (3 X, 7 Z) are
gauge-fixed by feed-forward destabilisers that commute with x̄ and z̄ (new `zdsim` op `FF`).

Tests (`tests/test_cult_growth.py`): determinism at d = 2, 4, 8; controls (no FF, injection mismatch); an infeasible
split; all single faults handled (d = 2: 7815 faults; d = 4: 117 375).

Sweeps (lenore; `runs/cultivation/sweep_grow35.json`, `sweep_grow35_d4_long.json` (10⁸ shots), `sweep_D5_r1.json`), d = 4:

| p | growth 3→5: accept, p_L | direct D5 (1 round): accept, p_L |
|---|---|---|
| 5e-4 | 0.672, 1.5e-7 (10 / 67M) | 0.752, 0 / 2.3M |
| 1e-3 | 0.452, 9.1e-7 (41 / 45M) | 0.565, 2.4e-6 (4 events) |
| 2e-3 | 0.204, 7.3e-6 | 0.320, 1.4e-5 |

Growth has ≈ 2–3× lower p_L at ≈ 20 % lower acceptance.  The 10⁸-shot growth run confirms the 1e-3 point
(9.1e-7 ± 16 %) and resolves 5e-4 at 1.5e-7 ± 32 % (p_L ∝ p^2.6 between the two points).  The grown states feed the distillation stages of
`MEAS_DISTILL.md` §5d.

**d = 8** (`sweep_grow35_d8_long.json`, 2×10⁷ shots per point, lenore, 2026-10-29): growth 3→5 acceptance 0.662 / 0.439 at p = 5e-4 / 1e-3, failures 0 / 13.2M accepted (< 1.7e-7) and 3 / 8.8M (3.4e-7, 90% CL 0.9–8.8e-7) — consistent with the d = 4 values above at the available statistics; the growth circuit, checks and signs are the same for every 2^k.
