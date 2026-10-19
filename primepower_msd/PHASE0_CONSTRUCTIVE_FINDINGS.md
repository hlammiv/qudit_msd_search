# Phase 0 execution — gates for the GR(4,r) constructive program

*Status 2026-06-29.  The three cheap, finite, decisive gates from `CONSTRUCTIVE_PLAN.md` §4.*

## 0a — Hardened anti-collapse certificate + engine audit  → SOUND

- **Engine audit PASS.** The Critique-1 "centralizer bug" flag is **refuted**: on [[5,1,3]]₄ all 6
  centralizer generators and all 4096 centralizer elements commute (symplectic form 0) with every
  stabilizer; stabilizers are isotropic; distance 3 confirmed.  `right_kernel_howell` centralizers are
  correct (consistent with the earlier 1200-trial validation).
- **Hardened certificate built.** The genuine-ring test is the **logical symplectic commutator**: a code
  is a genuine Z₄ logical qudit iff there is a logical pair X̄,Z̄ with ⟨X̄,Z̄⟩ **odd** (commutator i,
  order 4), not merely substrate clock order 4.  Both [[5,1,3]]₄ and [[8,1,2]]₄ have **odd logical
  symplectic commutators** ⇒ genuine Z₄ logical structure (PASS).  ((Z₄)² legitimately contains 2-torsion
  elements, so "some order-2 logical elements" is not a collapse; the odd commutator is the real test.)

## 0b — ω=i carrier existence  → **PASS**

Enumerating diagonal single-ququart phase gates U=diag(ζ_N^{q(x)}), level-3 + non-Clifford, by the
induced cubic phase ω = ζ_N^{(3rd finite difference of q)}:

| precision N | ω=i carriers (3rd-diff = N/4) | ω=−1 (soft collapse, 3rd-diff = N/2) |
|---|---|---|
| N=8 (2d) | **0** | 16 |
| N=16 (4d) | **16** | 16 |

- **Genuine ω=i carriers EXIST at 16th-root (4d) precision** — 16 of them, e.g.
  **q=[0,1,4,13] = diag(1, ζ₁₆, i, ζ₁₆¹³)** (the genuinely-quadratic level-3 gate; 3rd-diff = 4 mod 16 →
  ζ₁₆⁴ = i).  None at 8th-root precision (the genuine cubic needs 4d precision).  **The carrier route is
  alive.**
- Confirms the ω=−1 flag: the verified [[8,1,2]]₄ gate diag(1,−i,1,i) has 3rd-diff → ζ₁₆⁸ = −1 (a
  Z₂-sublevel "soft collapse").  Phase 1 should use a ω=i carrier (e.g. [0,1,4,13]) for genuine ring magic.

## 0c — vacuity probe (do RANDOM non-CSS d≥3 codes host an ω=i carrier transversally?)  → NULL (search)

**10 valid non-CSS d≥3 codes × 416 ω=i carriers (diagonal + C₁DC₂) = 4160 transversal tests → 0 hits.**
Consistent with V7 (0/780 on [[5,1,3]]₄, which already included these carriers).  Random non-CSS d≥3 codes
do **not** host an ω=i carrier transversally.  This is **"construct, don't search"** — the host code must be
*built structurally* (GR(4,r), Phase 1), not found by random sampling.  NOT a kill (not exhaustive over all
small codes), but confirms random search is exhausted.

## Phase 0 verdict — UNBLOCKED; commit to the constructive (GR(4,r)) program

| gate | result | meaning |
|---|---|---|
| **0a** | SOUND | engine audited (centralizer correct, "bug" refuted); hardened logical-commutator certificate built; both codes are genuine Z₄ logical qudits |
| **0b** | **PASS** | ω=i level-3 carriers exist (16 at 4d precision, e.g. diag(1,ζ₁₆,i,ζ₁₆¹³)) — the carrier route is alive |
| **0c** | NULL (search) | random non-CSS d≥3 codes don't host the carrier (0/4160) → construct, don't search |

Neither carrier-level gate killed the program; the host code just has to be **constructed**.  The
make-or-break now sits entirely in Phase 1: **Milestone 2** (prove a GR(4,r) distance bound that beats d=2
at k≥2) and the **M1 bridge theorem** (carry the CSS Bravyi–Haah carry-route cubic guarantee onto the
non-CSS GR module, at **mod-8** precision — the level-3 bookkeeping is invisible mod 4).  Recommended next:
derive M1 on the smallest GR(4,r), reproducing the verified k=1 codes as a sanity check, then test the
decisive experiment (build the high-degree-band non-CSS module, apply the ω=i carrier, require leakage 0 +
ω=i over all (Z₄)^k + hardened anti-collapse + d≥3).

