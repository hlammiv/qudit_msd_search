# Session findings — qudit MSD code discovery

Consolidated record of results from the qmsd toolkit work (extends arXiv:2510.10852, Saha–Prakash).
Author context: Henry Lamm (Fermilab); related prior work arXiv:2512.21874.

---

## 1. Tooling built and validated

- **`qmsd` package** — given a prime `p`, redoes the paper's punctured-Reed–Muller search; all 10
  published codes reproduced as unit tests; meet-in-the-middle (MITM) exact minimum distance certifies
  every oracle code (incl. the headline `[[519,106,5]]₅`).
- **MacWilliams-from-small-dual engine** (`qmsd/weightdist.py`) — exact weight distribution of
  `G0^perp` by enumerating the small dual `G0` (`q^{dim G0}` codewords) and Krawtchouk-transforming.
  Returns **certified minimum distance AND `A_d`**. Validated: reproduces all 6 small-dual oracle
  `A_d` exactly (760, 648, 512, 1700, 880, 3972); lifts the MITM `d≤6` cap (RS`[10,6]`/F₁₁ → d=7,
  A₇=1200); closes the "exact `A_d` for large codes" gap (`[[667,62,4]]₃` → A_d=3972). 243 tests pass.
  *Known caveat:* `A_d_logical` is hardcoded `= B_d` (no in-engine `Gp` stabilizer filter) — exact
  across the validated regime, but would over-report `A_d` on a novel code with weight-`d` stabilizers.

## 2. The Singleton feasibility filter (general tool)

γ<1 requires `d > n/k` while Singleton gives `d ≤ D−k+1` (`D = dim RM_p(r_max,m)`). Combining:
> **γ<1 is impossible unless `D > 2√(pᵐ) − 2`.**

First-feasible `(p,m)`: p=3→m=5; **p=5→m=4 (exactly where the paper's `[[519,106,5]]₅` lives)**;
p=7→m=3 (loose); **p=11,13→m=2 (tiny)**; p=17,19→m=2. Necessary-but-loose (punctured RM is far from
MDS), so passing it doesn't guarantee findability — but it cleanly rules out small-m and tells you
which `(p,m)` to search.

## 3. NEW RESULT — first γ<1 code from search: `[[234,55,5]]₁₇` (γ = 0.8997)

- Punctured `RM₁₇(10,2)`; full rank, k=55 logical qudits, **d=5 certified** (MITM + `d_max=7` recheck).
- ~330× smaller than the paper's analytic p=17 code `[[77540,5981,15]]`; below the Reed–Solomon
  threshold (p=17 < 23) so unreachable by m=1 RS. Recorded with its puncture set in `NEW_CODES.md`.
- Found via a **predicted** trend: best-γ at m=2 falls monotonically with p — p=11→1.126,
  p=13→1.059 (one distance short), **p=17→0.900 (crosses)**. The reasoning chain
  Singleton-filter → trend-extrapolation → search → certification predicted exactly where to look.
- p=19 m=2: γ<1 also reachable (lighter-shortened codes have `d≥7`, provable by MITM returning
  "no ≤6 dependency"), but the exact-d codes need higher tooling — not locked to a concrete code.

## 4. Methodological insights

- **MITM reach is bounded by `dim(G0)=D−k`**, not by `n` — the int64 syndrome encoder holds `p^{dim G0}`.
  So high p is computable only at **high k** (heavy shortening, small `dim G0`) — which is *also* where
  γ<1 lives (γ<1 needs only d=4–5 there). Search the high-k window.
- **`d ≥ d_cap+1` is itself a γ<1 certificate**: when the MITM finds no dependency ≤ its cap, that
  *proves* `d ≥ cap+1`, giving a proven upper bound `γ ≤ log(n/k)/log(cap+1)`. Don't discard those.

## 5. Prime-power (Galois) qudits — #4 investigated, verdict: reduces to known

(Workflow `wmri52kwt` → `GALOIS_QUDIT_INVESTIGATION.md`.) The classical/coding layer generalizes
cleanly to any `GF(q)` (power-sum, triorthogonality `3r<m(q−1)`, RM duality — all verified in galois).
But the **single-qudit T-gate does not lift**: magic needs ring `Z_{p^e}` precision while field RM
yields only Clifford p-th-root phases (and the cubic monomial Frobenius-degenerates). What survives is
the **multi-qudit trace-CCZ** on the Galois qudit (= e prime-p qudits). Already published (incl.
arXiv:2512.21874 — the user's group); d=4=GF(4)=two-qubit MSD, etc. The paper's prime-only scope is
**principled**. → See also the single-qudit no-go below.

## 6. Galois single-qudit magic — settled boundary

A `GF(p^e)` **field** code cannot give a useful single-qudit magic state on a native **ring `Z_{p^e}`**
qudit: field transversal phases are p-th roots (Clifford on one Galois qudit), and `GF(p^e)` has no
ring-compatible integer lift (which is what supplies single-qudit magic in the prime case). "Single
Galois qudit" = e prime-p qudits → reduces to multi-qubit/qutrit MSD. Useful single-(ring)-qudit magic
for prime-power d is an open **ring-code** (`Z_{p^e}` / Galois-ring / chain-ring) question, not a
field-code one. (User is pursuing neither ring codes nor field-distill-then-convert at present.)

## 7. m=7 qutrit — the engine works; the wall moved to puncture-set design

- The MacWilliams engine **certifies distances at m=7** (the original "compute wall" is gone).
- But **random puncturing collapses the distance to d=1–2** (`RM₃(9,7)` has ~5.8×10⁹ weight-18
  codewords; heavy random puncturing guts one). The analytic Manhattan family at m=7 gives only γ=1.86.
- **Key structural result:** a **cap puncture gives provable `d≥10`** — a cap meets any 2-flat in ≤4
  points, and the weight-18 min codewords are 2 parallel 2-flats, so ≥10 survive. Hence
  `[[2187−k, k, ≥10]]₃` with **γ<1 for k≥200** (k=210→0.974; k=248 (max cap)→0.893).
- **`d≥10` rigor:** proven for weights `[18,45]` (Kasami–Tokura flat structure + the `5/9·w≥10` bound)
  and `w ≥ k+10` (trivial). Open band `(45, k+10)` reduces to "the support contains ≥10 disjoint lines"
  (true for all flat/algebraic supports; a cap-like low-degree support is not known to exist). The
  **max-cap (k=248) code is fully certifiable on a cluster** (`3²⁶` MacWilliams). Cap construction in
  progress (workflow `wbwcnji1k`).

## 8. Open threads / next steps

- **Construct the AG(7,3) cap** (≥200 for γ<1; ideally 248) → the m=7 qutrit γ<1 code (in progress).
- Close the `(45, k+10)` band ("low-degree RM supports are line-rich") or cluster-certify the max-cap code.
- **Harden the engine's `A_d` filter** (add the `Gp` stabilizer check) before trusting `A_d` on novel codes.
- **(γ, C, A_d) Pareto search** (Tier 4): multi-objective archive instead of a single scalar; the engine
  supplies `A_d` for free at p=3,5.
- Lock the **p=19** family code; revisit p≥23 Reed–Solomon regime (Tier 2).
