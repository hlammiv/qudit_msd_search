# Prime-Power-Qudit (Galois-Qudit) Generalization of the Punctured-RM MSD Construction

**Question.** arXiv:2510.10852 (Saha–Prakash) builds qudit magic-state-distillation (MSD)
codes by puncturing Reed–Muller codes over the prime **field** F_p, yielding triorthogonal
CSS codes with a transversal non-Clifford (T-type) gate, restricted to prime dimension p.
Can the *whole* construction be lifted to prime-power dimension q = p^e (d = 4, 8, 9, 16, 25, 27)
using RM / Reed–Solomon codes over the Galois **field** GF(q)?

**One-sentence answer.** The *classical / coding* layer generalizes cleanly over any GF(q),
and a *transversal non-Clifford gate does exist* on the Galois-qudit (GF(q)-linear = e prime-p
qudits) — but it is the **3-qudit trace-CCZ |x,y,z⟩→ω^{Tr(xyz)}**, not a single-qudit Howard–Vala
T-gate (which **degenerates to a Clifford** over GF(q)). This object is **already published**,
the small hardware-native cases **collapse** (d=4 = two-qubit MSD, d=9 = two-qutrit MSD), and
native single-system d=p^e devices are **ring Z_{p^e} qudits the field construction cannot serve**.
Net verdict: **the construction is well-defined but largely reduces to known results.**

---

## Per-layer summary

### Layer 1 — Clifford hierarchy / transversal-T (`generalizes: with-caveats`; verified: **solid**)

- **Field layer generalizes.** The prime-p transversal-T proof uses primality *only* through
  the power-sum identity ∑_{x∈F} x^a = 0 unless (q−1)|a. This holds verbatim over any GF(q)
  (verified, q=4,8,9,16,25,27). Consequently triorthogonality of RM_q(r,m) for **3r < m(q−1)**
  holds field-by-field (verified directly: pair- and triple-sums vanish; the boundary case
  GF(4),m=2,r=2 with 3r=6=m(q−1) **fails** the triple condition, matching the strict inequality).
- **A transversal non-Clifford gate exists — on the Galois qudit.** Cui–Gottesman–Krishna
  (arXiv:1608.06596) characterize diagonal Clifford-hierarchy gates and state explicitly that
  prime-power p^r works because its Clifford group/hierarchy is **isomorphic to that of r prime-p
  qudits** "with a standard choice of Pauli group." Verified that the additive group of GF(p^e)
  is (Z_p)^e (every nonzero element has additive order p), so a GF(q) clock/shift factorizes into
  e prime-p clock/shifts: **a GF(q) qudit literally = e prime-p qudits.**
- **Caveat (load-bearing).** This is the *Galois* (GF(q)-linear) qudit, **not** a single monolithic
  **ring Z_{p^e}** qudit. The RM/RS machinery lands on the field, hence on the Galois qudit; it does
  **not** hand you a transversal T on a true d=p^e cyclic-clock device (that needs the ring picture,
  Rengaswamy–Calderbank–Pfister arXiv:1902.04022, binary-ring case).

### Layer 2 — Concrete small dims d=4, 8, 9 (`generalizes: partially`; verified: **solid**)

- **Single-qudit T BREAKS** for a sharp reason: the level-3 magic phase is a **ring Z_{p^m}** object
  (8th roots in char 2, 9th roots in char 3), not a field GF(q) object. RM-over-GF(q) triorthogonality
  only ever produces field-linear phases ω^{Tr(...)} with ω a p-th root — **all Clifford** (verified:
  every ±1 diagonal gate on the GF(4) two-qubit qudit is Clifford; genuine magic needs 8th roots).
- **Frobenius degeneracy** (verified): x→x³ is additive in char 3 and x→x² additive in char 2, so the
  obvious single-qudit cubic phase ω^{Tr(c·x³)} collapses to a Clifford/separable gate. d=9 is **not**
  "two independent qutrits running qutrit-MSD," but the naive field-cubic phase gives nothing.
  *(Prime dimension is special: a field element lifts canonically to 0..p−1, so ∑xyz over F_p carries
  into the Z_{p²} phase that creates magic; GF(p^e) has no ring-compatible integer lift.)*
- **What survives:** the **multi-qudit trace-CCZ** route. ω^{Tr(xyz)} is genuinely trilinear in three
  independent GF(q) elements (no monomial collapse), and triorthogonality ∑_i(xyz)_i=0 over GF(q) is
  converted to a Z_p phase cancellation by the trace: Tr(0)=0. Golowich–Guruswami (arXiv:2408.09254)
  build transversal CCZ over **arbitrary prime power q including q=2** via RS/AG codes.
- **Classical layer verified:** character sum (= q−1 = −1 when (q−1)|a) and generalized RM duality
  RM_q(r,m)^⊥ = RM_q(m(q−1)−r−1, m) confirmed exhaustively for q=4 (m=2) and q=9 (m=2).
- *Nuance:* "all ±1 gates Clifford" is **formalism-dependent** — true on the two-qubit/Galois qudit
  (the right reading), false on the cyclic Z_4 qudit.

### Layer 3 — Novelty / usefulness / hardware (`generalizes: partially`, `promising: false`; verified: **solid**)

- **Mostly NOT novel.** (a) Galois-qudit RM/stabilizer codes are an established catalogued family
  (Error Correction Zoo). (b) "triorthogonality over GF(q) → transversal CCZ" is already published:
  Golowich–Guruswami (arXiv:2408.09254); Cervia–Lamm–Liu–Murairi–Zhu (arXiv:2512.21874) *explicitly*
  generalize triorthogonal matrices to Galois qudits q=2^m and decompose each into m qubits via a
  normal/trace basis; Krishna–Tillich (PRL 123, 070507) use RS-over-GF(q) for transversal CCZ with
  γ→0. (c) The only unwritten piece is the *specific* 2510.10852 punctured-GRM sublog recipe — incremental.
- **d=4 collapse: YES.** The normal-basis map (explicit in 2512.21874) sends each GF(4) qudit to 2 qubits
  and the trace-CCZ to structured qubit CCZ/CZ. So d=4,8,16 → multi-qubit MSD and d=9,27 → multi-qutrit MSD.
- **Usefulness: no proven advantage** over e separate prime-p qudits at the hardware-native q. The GF(q)
  code's logical gate is the *entangled* trace-CCZ — useful only if you actually run GF(q)-arithmetic
  algorithms. γ→0 requires q→∞ (RS/AG), where the *alphabet*, not the punctured-RM trick, does the work.
- **Hardware: decisive.** Native single-system d=p^e devices are **ring Z_{p^e}** qudits: transmon ququart
  phase gate Z_4=diag(1,i,−1,−i) (arXiv:2304.11159, 2303.04796), GKP ququart in one oscillator
  (Brock et al., Nature s41586-025-08899-y / arXiv:2409.15065), spin-3/2 molecular qudits. The GF(q)
  construction needs GF(q)-linear qudits = e prime-p qudits with imposed field structure. So it serves
  **qubit/qutrit registers, not native ring ququarts.** The paper's omission of prime powers is *principled*.

---

## Verdict

### (1) Is the GF(q) triorthogonal-MSD construction WELL-DEFINED — and where does it break?

**Well-defined as a classical/Galois-qudit object; the single-qudit T-gate is where it breaks.**

- ✅ **Coding layer (clean):** power-sum ∑_{x∈GF(q)} x^a = 0 unless (q−1)|a; triorthogonality of
  RM_q(r,m) iff **3r < m(q−1)**; duality RM_q(r,m)^⊥ = RM_q(m(q−1)−r−1, m). All verified in
  python+galois for q = 4, 8, 9, 16, 25, 27. Because GF(q) is a field, every step of the prime-p
  proof goes through monomial-by-monomial.
- ✅ **A transversal non-Clifford gate exists** on the Galois qudit (= e prime-p qudits): the
  **trace-CCZ ω^{Tr(xyz)}**, level-3, via the CGK isomorphism (prime-power = e prime-p qudits).
- ❌ **Where it breaks — the ring-vs-field / transversal-T question.** The single-qudit Howard–Vala
  π/8 gate does **not** lift. Magic requires **Z_{p^m} ring precision** (8th roots char 2, 9th roots
  char 3), but RM-over-GF(q) only produces **field-linear** phases ω^{Tr(·)} with ω a p-th root,
  all of which are **Clifford**. The natural cubic monomial ω^{Tr(c·x³)} additionally **degenerates**
  because cube/square are Frobenius-**linear** in char 3 / char 2. GF(p^e) has no ring-compatible
  integer lift, so the field-linear triorthogonality and the ring precision needed for level-3 never
  align *on a single qudit*. The construction is therefore well-defined only as a **multi-qudit
  trace-CCZ** code, not a single-qudit-T code; and it does **not** serve a monolithic ring Z_{p^e}
  qudit at all.

### (2) Is it NOVEL, or does it reduce to known results?

**Largely reduces to known results.**

- The Galois-qudit / GF(q)-triorthogonal / transversal-CCZ framework is **already published**
  (Golowich–Guruswami 2408.09254; Cervia–Lamm et al. 2512.21874 with explicit Galois→qubit reduction;
  Krishna–Tillich PRL 123 070507; ECZ Galois-qudit RM entry).
- **d=4 = GF(4) = two qubits**, so a GF(4) triorthogonal MSD code **is** a two-qubit MSD code via the
  normal-basis decomposition; d=8,16 → multi-qubit MSD, d=9,27 → multi-qutrit MSD. The Galois view
  adds a *structured/constrained* logical gate, not a new resource.
- The only genuinely-unwritten item is the *specific* punctured-GRM-over-GF(p^e) sublog recipe of
  2510.10852 — an incremental instantiation whose asymptotic payoff (γ→0) is already known to come
  from large alphabet, not from this construction.

### (3) Is it USEFUL (overhead / hardware)?

**No demonstrated advantage; hardware mismatch is decisive.**

- No proven yield (γ) advantage over the honest competitor (e decomposed prime-p qudits running
  prime-qudit MSD) at the hardware-native set {4,8,9,16,25,27}. γ→0 needs q→∞.
- The produced magic resource is the **trace-CCZ over GF(q)** — an *entangled* combination of
  prime-qudit gates, useful **only** for GF(q)-arithmetic quantum algorithms, and it outputs a
  **Galois (e prime-qudit register)** resource, not a native ring-ququart resource.
- Native single-system prime-power hardware (transmon/GKP ququart, molecular spin-3/2) are **ring
  Z_{p^e}** qudits; the **field-based RM/RS machinery cannot address them.** Serving them would
  require a separate **ring (Galois-ring / chain-ring) MSD theory** — a different, harder paper.

**Overall: `likely-reduces-to-known`.** The clean math generalizes, but the result is mostly
already in the literature, the small dims collapse to multi-qubit/qutrit MSD, and the
hardware-relevant native devices are ring qudits the construction cannot serve.

---

## If pursued anyway — concrete construction sketch + qmsd additions

The only defensible deliverable is a **clarifying / negative empirical benchmark** documenting the
reduction, not a new capability. Because `galois` already supports GF(p^e), the code changes are small.

**Construction (what is actually well-defined):**
1. Work over GF(q), q=p^e. Build generalized-RM generators by evaluating monomials
   x_1^{a_1}···x_m^{a_m} (0 ≤ a_i ≤ q−1, ∑a_i ≤ r) at all q^m points of GF(q)^m.
2. Triorthogonality holds iff **3r < m(q−1)**; dual is RM_q(m(q−1)−r−1, m). Puncture/shorten exactly
   as in the prime case to get the CSS code (G0 ⊆ G' ⊆ G0^⊥).
3. The transversal logical gate is the **3-qudit trace-CCZ** |x,y,z⟩→ω^{Tr(xyz)} (ω = e^{2πi/p}),
   NOT a single-qudit T. Triorthogonality ⇒ Tr(∑_i(xyz)_i)=Tr(0)=0 ⇒ logical phase preserved.
4. For reporting, decompose each Galois qudit into e prime-p qudits via a trace/normal basis to expose
   the equivalent multi-qubit/qutrit code (this is the "reduces to known" demonstration).

**qmsd additions (mirroring the existing `field.py` / `reedmuller.py` / `triorthogonal.py` split):**
- `field.py`: generalize `field_power_sum(p,a)` → `field_power_sum(q,a)` returning q−1 iff (q−1)|a,
  and `GFp(p)` → `GFq(q)` thin wrapper over `galois.GF(q)` (works for prime powers).
- `reedmuller.py`: replace the `p` threshold m(p−1) with **m(q−1)** in `r_max`, `r_tilde`; let the
  per-variable degree cap be q−1 (galois reduces mod x^q−x). Duality test already structurally identical.
- `triorthogonal.py`: `is_triorthogonal(basis, q)` over GF(q) — pair/triple entrywise sums in-field
  (already multilinear; just swap the field).
- **New** `galois_gate.py`: implement the trace-CCZ logical-gate check (and the Frobenius-degeneracy
  guard that flags single-qudit cubic phases as Clifford) so the tooling does not mislabel it a T-gate.
- **New** benchmark in `asymptotics.py`/`distillation.py`: tabulate γ = log_d(n/k) and single-round cost
  for q∈{4,8,9,16,25,27}, compared against (a) e-copies prime-p MSD and (b) RS-over-GF(q). Expected
  outcome: small-q field structure does **not** beat the prime-qudit baseline — settles "worth doing."

---

## Verified references (cited only where confirmed real)

- S. X. Cui, D. Gottesman, A. Krishna, *Diagonal gates in the Clifford hierarchy*, PRA **95**, 012329 (2017), arXiv:1608.06596. — diagonal hierarchy gates are p^m-th-roots^poly (ring precision); prime-power p^r ≅ r prime-p qudits (verbatim in full text).
- M. Howard, J. Vala, *Qudit versions of the qubit π/8 gate*, PRA **86**, 022316 (2012), arXiv:1206.1598. — prime-dimension single-qudit T; no prime-power version.
- L. Golowich, V. Guruswami, *Asymptotically Good Quantum Codes with Transversal Non-Clifford Gates*, STOC 2025, arXiv:2408.09254. — transversal CCZ over arbitrary prime power q incl. q=2, via RS/AG codes. *(Authors are Golowich–Guruswami; not "P. Nguyen" — a separate real paper arXiv:2408.10140, The Quynh Nguyen, "Good binary quantum codes with transversal CCZ," is a distinct work.)*
- E. T. Campbell, H. Anwar, D. E. Browne, *MSD in all prime dimensions using quantum RM codes*, PRX **2**, 041021 (2012), arXiv:1205.3104. — prime-qudit RM MSD baseline.
- A. R. Calderbank, E. M. Rains, P. W. Shor, N. J. A. Sloane, *Quantum Error Correction via Codes over GF(4)*, IEEE Trans. Inf. Theory **44**, 1369 (1998), arXiv:quant-ph/9608006. — GF(4)↔binary is Clifford/stabilizer-level only; no magic transport.
- N. Rengaswamy, R. Calderbank, H. D. Pfister, *Unifying the Clifford Hierarchy via Symmetric Matrices over Rings*, PRA **100**, 022304 (2019), arXiv:1902.04022. — ring (Z_{2^k}) picture; the framework needed for a true monolithic ring qudit.
- J. C. M. de la Fuente / Cervia, Lamm, Liu, Murairi, Zhu, *Magic State Distillation using Asymptotically Good Codes on Qudits*, arXiv:2512.21874 (2025). — explicit Galois-qudit (q=2^m) triorthogonal generalization + transversal CCZ + normal-basis Galois→qubit decomposition.
- A. Krishna, J.-P. Tillich, *Towards Low Overhead Magic State Distillation*, PRL **123**, 070507 (2019), arXiv:1811.08461. — RS-over-GF(q) transversal CCZ/U with γ→0 as q grows.
- *Quantum Codes with Addressable and Transversal Non-Clifford Gates*, arXiv:2502.01864 (2025). — diagonal hierarchy gates from polynomials over F_q.
- S. Prakash et al., *Qutrit and ququint magic states*, arXiv:2003.07164. — state of the art is prime d=3,5; negative result for d>3.
- N. de Silva, O. Lautsch, *The Clifford hierarchy for one qubit or qudit*, arXiv:2501.07939. — single-qudit hierarchy structure (prime-dimension scope).
- S. Prakash, S. Saha, *Low Overhead Qutrit Magic State Distillation*, Quantum **9**, 1768 (2025), arXiv:2403.06228. — prime-qutrit MSD baseline (d=9,27 competitor).
- N. Earnest / Seifert et al., transmon ququart, arXiv:2304.11159; *Emulating two qubits with a four-level transmon*, arXiv:2303.04796. — native d=4 transmon = ring Z_4 qudit.
- A. Z. Brock et al., *Quantum error correction of qudits beyond break-even*, Nature **642** (2025), s41586-025-08899-y / arXiv:2409.15065. — logical qutrit & ququart in a single GKP oscillator (ring qudit).
- S. Saha, S. Prakash, *Sublogarithmic Distillation in all Prime Dimensions using Punctured Reed-Muller Codes*, arXiv:2510.10852 (2025). — the target paper; prime p only, prime-power omitted.
- Error Correction Zoo, *Galois-qudit / prime-qudit Reed–Muller code* entries — confirm the Galois-qudit RM family is established and punctured GRM gives transversal Clifford-hierarchy gates + MSD.

*All math claims (power-sum, triorthogonality 3r<m(q−1), RM duality, additive group (Z_p)^e, Frobenius
degeneracy, GF(4)/qutrit Clifford-vs-magic) independently reproduced in python+galois v0.4.6 across all
three layers, including the boundary-failure case.*
