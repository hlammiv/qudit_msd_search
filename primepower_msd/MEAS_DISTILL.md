# Measurement-based (HHPW) distillation over Z_{2^k}: γ ≈ 1.42 at order 7 with ≈ 2.1k qudits

*2026-10-20/21.  Code: `meas_distill.py`, `outer_code.py`, `outer_hyper.py` (+`_run`), `inner_search2.py`, `meas_sim.py`; tests: `tests/test_meas_distill.py`, `tests/test_outer_code.py`;
data: `runs/meas_distill/`.*

## 1. Mechanism
Haah–Hastings–Poulin–Wecker (arXiv:1703.07847) distil by *measuring* the check H = T X T† on k noisy magic
states held in the logicals of an inner [[n,k,D]] code (T^{-1} · controlled-X^{⊗n} · T on the data: 2n noisy T
per check), and arrange the checks by an outer parity-check matrix M (inputs × checks, column weight s, row weight
k) with 2|Mv| + |v| ≥ D for all nonzero v.  Lemma 4.4 then gives D-th order suppression with
**T per output = 1 + 2 s n/k** (= D + (D−1)(n/k − 1) for s = (D−1)/2), **γ = log(T/out)/log D → 1** as n/k → 1.
Only *level-2* (Clifford) transversality is needed, so the quadruply-even condition behind the γ_T ≈ 2 limit
of transversal-T codes (HAILMARY_FINDINGS.md) does not apply.

## 2. Z_d transfer (verified for d = 4, 8, 16)
* T X T† = X S_{2d}^{a} · (Pauli), with a odd ⇔ T escapes (CULTIVATION_PLAN §10).
* T^{-1} = N T N† · (Clifford), with N: x ↦ −x; q(−x) + q(x) is a Clifford table, so T and T^{-1} each cost one
  magic state (tested).
* **Signed weakly self-dual codes.** A 0/1 X-stabiliser matrix S with signs τ ∈ {±1}ⁿ, Z-stabilisers S·diag(τ),
  Q(s) = Σ_{i∈s} τᵢ ≡ 0 (mod 2d) and B(s,t) = Σ_{i∈s∩t} τᵢ ≡ 0 (mod d).  Then ∏ T^{τᵢ} X T^{−τᵢ} on all n qudits
  (the all-ones check) is a logical Clifford.  Exact integer balance (Q = B = 0 over ℤ) works for **every** d at once.
* **Normal basis.** If some logical has odd weight (HHPW Thm 3.4), Gram–Schmidt with odd pivots gives logical reps
  N_j with G = diag(odd).  The all-ones check then acts as ∏_j X_j^{c_j} S_j^{Q_j} with c_j = (signed weight of
  N_j)/Q_j odd: an independent escaping check on every logical.  This was verified exactly (hosting + separability of
  the induced logical) on [[7,1,3]], [[15,7,3]] (d = 4, 8), [[31,1,7]] and [[63,21,7]].
* **Sub-checks need care.** A check X^u with non-0/1 u, or on a subset of logicals, needs u-weighted balance
  Σ_{s∩t∩u} τ = 0, which the family below does not have.  Hence every check uses all k logicals (row weight w = k),
  as in HHPW Lemma 4.15.

## 3. Inner codes: signed punctured RM (exact balance, all 2^k)
S = monomials of degree 1..r on F₂^m∖{0}, τ(x) = (−1)^{x₁+…+x_{2r+1}} (m ≥ 2r+1): every product of two generators
has degree ≤ 2r, so its signed sum vanishes over ℤ.  Parameters [[2^m − 1, 2^m − 1 − 2Σ_{i≤r} C(m,i), D]].

**Distance theorem (all Z_{2^k}): D = 2^{r+1} − 1.** For a Z-logical z, take the 2-adic layers of z.  Either
z mod 2 ∉ S_bin, so |z| ≥ d(punctured RM(m−r−1,m)) = 2^{r+1} − 1; or z mod 2 is a nonzero stabiliser, so
|z| ≥ 2^{m−r} > 2^{r+1} − 1; or z ≡ 0 mod 2, and the argument recurses on z/2.  The upper bound comes from (d/2)·(a
binary minimum word).  d_X = d_Z via z ↦ τ∘z.  The Z₄/Z₈ MITM confirms d > 6 for [[63,21,7]] at d = 2, 4, 8.

| inner code | D | n/k | T/out (s=(D−1)/2) | γ |
|---|---|---|---|---|
| [[63,21,7]] | 7 | 3.00 | 19.0 | 1.513 |
| [[127,71,7]] | 7 | 1.79 | 11.7 | 1.265 |
| [[255,183,7]] | 7 | 1.39 | 9.4 | 1.149 |
| [[511,253,15]] | 15 | 2.02 | 29.3 | 1.247 |
| [[1023,673,15]] | 15 | 1.52 | 22.3 | 1.146 |
| [[8191,6007,31]] | 31 | 1.36 | 41.9 | 1.088 |

(D = 3 members [[2^m−1, 2^m−1−2m, 3]] need s = 2, so T/out = 1 + 4n/k, γ ≥ 1.46.)

## 4. Outer codes (the real size constraint)
HHPW's sensitivity proof (Lemma 4.15) uses unique-neighbour leaves of high-girth Tanner graphs, so it holds for
**every Z_d alphabet**: a leaf check sees a single nonzero error value.

* **s = 2 (order 5):** the outer code is a w-regular graph (checks = vertices, inputs = edges).  The exact
  criterion is girth ≥ 5 (double edges, triangles and 4-cycles are the only violating configurations of ≤ 4 inputs;
  tested exhaustively against the girth test).  Inputs = m·w/2 ≥ w(w²+1)/2 (Moore bound).
* **s = 3 (order 7):** a 3-uniform hypergraph with point degree w, free of Pasch-, grid- and other dense
  configurations of ≤ 6 inputs.  Lemma 4.15 guarantees this at Tanner girth ≥ 14; girth ≥ 10 is very probably
  enough.  For w = 21 even girth 10 needs ≳ 12k inputs, so this case is open below 10k.

## 5. Concrete protocol under 10k qudits (any d = 2^k)
* inner [[63,21,7]] (signed RM(2,6)), checks of row weight 21;
* outer: 21-regular girth-5 graph on 722 checks.  The biaffine plane over F₁₉ (19-regular, girth 6) plus a greedy
  girth-preserving 2-factor gives 7581 inputs (`runs/meas_distill/outer_w21_girth5_722.json`, girth checked with
  networkx).  The pure biaffine fallback (q = 23, |X| = |M| = 21) is 21-regular with girth 6 and 10 143 inputs;
* **order-5 suppression, T per output = 1 + 2·2·63/21 = 13, γ = log 13/log 5 = 1.594 (< 1.68)**;
* qudits ≈ 7581 (inputs/outputs) + 63 (inner) + 63 (cat) + stabiliser ancillas ≈ 7.8k; one run outputs 7581
  T states.  Error model: only T gates noisy (standard γ accounting); Clifford noise is not included.

## 5b. Order 7 (s = 3) — current best, any d = 2^k
**Outer hypergraphs** (`outer_hyper.py`): exact sensitivity certifier (DFS freeze/add with pruning; agrees with
brute force on 276 random + 35 perturbed + 8 constructed instances) and `construct2` (parallel root scan + random
edge switching until no root has a violator).  Certified (3, w) codes, lenore, `runs/meas_distill/outer7_w*_n*_s*.json`:

| w | certified inputs (checks) | not reached |
|---|---|---|
| 21 | 1470 (210), 2100 | 1050 (stuck at ~40 violating roots after 180 rounds) |
| 26 | **2080 (240)**, 2340, 2600, 3120, 3640 | — |
| 35 | 4900 | 3500 (running) |
| 40 | 5200, 6400 | — |
| 50 | **9500 (570)**, seeds 1 and 2 (`--no4`: 4-cycle pre-pass first, then 1134 / 1163 switches, 9.3 / 13 h) | 9500 without the pre-pass (2694 violating roots at round 0) |
| 71 | — | 7100, 9940 (round 0 only) |

Minimum size grows ≈ 3–4 w², so w ≲ 50 under 10k inputs.  At w ≥ 40 the dominant violators are Tanner 4-cycles
(two inputs sharing two checks: U = 2, t = 2 → 2U + t = 6 < 7), ≈ 18 n²/(2 m²) of them; `remove_4cycles` switches
them out first (rejecting switches that create new ones), after which the certifier repair converges in ~10 rounds.

**Inner codes with D ≥ 7** (`inner_search2.py`: anneal on signed flat collections, objective = number of
≤ 3 + ≤ 3 syndrome collisions, which is 0 iff d_bin(S^⊥) ≥ 7; aligned forbidden flats cap n):
`runs/meas_distill/inner7b/best_inner7.json` — all exactly balanced, normal, Z₄ d > 6, Z₈ d > 5, D ≥ 7 by the 2-adic theorem.

| inner | T/out = 1 + 6n/k | γ (order 7) | outer w | best certified outer |
|---|---|---|---|---|
| [[63,21,7]] signed RM(2,6) | 19.00 | 1.513 | 21 | 1470 inputs |
| **[[64,26,≥7]]** | **15.77** | **1.417** | 26 | **2080 inputs** |
| [[96,40,≥7]] | 15.40 | 1.405 | 40 | 5200 inputs |
| [[112,50,≥7]] | 14.44 | 1.372 | 50 | **9500 inputs** (fault distance 7 confirmed; ≈ 1927 p⁷, §5c) |

[[63,21,7]] hosting re-verified at d = 4, 8 (with distance), 16, 32, 64 (balance, normal basis, all logical S
exponents odd): `runs/meas_distill/rm63_alld.log`.

**Recommended protocol: [[64,26,≥7]] inner + (3,26) outer on 2080 inputs (240 checks).**  Order 7, 15.77 T per
output, γ = 1.417, ≈ 2080 + 64 + ancillas ≈ 2.2k qudits, valid for every d = 2^k.
[[112,50]] + 9500 reaches γ = 1.372 under 10k qudits, but γ counts T gates per output at acceptance 1; with its larger
leading coefficient and acceptance it is never cheaper per accepted output above p ≈ 10⁻⁶ (§5c), so the 2080 protocol stays
the recommendation.

## 5c. Fault model, exact leading order, Monte Carlo (`meas_sim.py`)
Twirled model, Cliffords perfect: inputs carry Z^e (prob p_in), each T gate is followed by Z^a (prob p_T), values
uniform nonzero.  Z commutes with T, so a layer-1 fault shifts the check outcome and stays on the data; a layer-2
fault only stays on the data.  Bookkeeping is linear over Z_d (outcome, X-syndrome S z, logical shift ⟨N_j, z⟩).
* **Validation:** toy instances ([[63,21,7]] with s = 1; [[15,7,3]] with a 7-regular graph on 10 checks).  Exact
  enumeration equals a numba brute force over all ≤ 3-fault patterns (952 = 952); Monte Carlo matches the
  leading-order prediction (241 vs 271, 167/170 vs 167, 22/36 vs 26/39 at d = 8).
* **Fault distance 7 confirmed exactly** for the real protocols: no input set with |v| + 2·#unique ≤ 6 and no inner
  logical below weight 7.  Weight-7 failure patterns (d = 4):

| protocol | single input + 3 flip pairs | 5-input sets | inner logical in last check | per-output p_out |
|---|---|---|---|---|
| [[63,21]] + 1470 | 1 102 707 270 | 39 942 | 970 064 | ≈ 343 p⁷ |
| [[64,26]] + 2600 | 2 044 723 200 | 66 688 | 825 192 | ≈ 360 p⁷ |
| [[64,26]] + 2080 | 1 635 778 560 | 108 288 | 598 328 | ≈ 360 p⁷ (MC at 1e-6, 1e-5: 0 fails, acceptance 96.8 %, 72 %) |
| [[112,50]] + 9500 | 40 040 448 000 | 942 256 | 98 304 | ≈ 1927 p⁷ (MC at 1e-6, 1e-5: 0 fails, acceptance 87.2 %, 25.4 %) |

  Dominant term = N·(d−1)·n³ (one bad input hidden in each of its 3 checks by a layer-1/layer-2 pair Z^a, Z^{−a} on
  one inner qudit), so p_out ≈ n³ p⁷/(d−1)⁶ (d = 8: ≈ 2 p⁷).
  This is why the lower-γ [[112,50]] + 9500 protocol (`runs/meas_distill/count_112_50_n9500_d4.log`, enumeration
  parallel over root inputs, 24 cores, 1.2 h) is not actually cheaper: its coefficient is 5.4× larger (112³ vs 64³)
  and its 4.2× more fault locations (L = 137 180) cost acceptance — at p = 10⁻⁵ it is 57 T per accepted output vs 22 for
  [[64,26]] + 2080; the two tie only below p ≈ 10⁻⁶ (16.6 vs 16.3), where both are far past any fidelity need.
* **Acceptance = exp(−L p)** with L = N(1 + 2sn/k) fault locations (every fault is detected; MC agrees to 3 digits).
  [[64,26]] + 2600: 96 % at p = 10⁻⁶, 66 % at 10⁻⁵.  The protocol wants inputs at ≲ 10⁻⁵ (for example cultivated T
  states); at p = 10⁻⁴ acceptance collapses (6 % for [[63,21]] + 1470).  At such inputs, output fidelity (≈ 10⁻³³) is far
  beyond need, so the binding constraint is acceptance, not order.

## 5d. Noisy inputs: concatenation, not big outer codes (`msd_planner.py`)
Big order-7 rounds have L ≈ 16·N fault locations, so above p ≈ 10⁻⁵ their acceptance exp(−Lp) collapses.
**Order-3 stages** fix this: inner signed RM(1,m) [[2^m−1, 2^m−1−2m, 3]], outer = complete graph K_{k+1} (inputs = edges,
s = 2; order 3 needs only a simple graph).  Exact d = 4 leading-order coefficients (`runs/meas_distill/order3_d4.json`):

| stage | inputs | T/out | L | p_out per output |
|---|---|---|---|---|
| [[7,1,3]] | 1 | 29.0 | 29 | 2.37 p³ |
| [[15,7,3]] + K₈ | 28 | 9.57 | 268 | 0.84 p³ |
| [[31,21,3]] + K₂₂ | 231 | 6.90 | 1595 | 1.24 p³ |
| [[63,51,3]] + K₅₂ | 1326 | 5.94 | 7878 | 2.46 p³ |
| [[127,113,3]] + K₁₁₄ | 6441 | 5.50 | 35397 | 5.08 p³ |

**Partial discard fails** (discarding only the inputs of a failed check): a single T fault in a failed check spreads
logical shifts onto several of its inputs.  Those inputs then cancel a real error elsewhere in both of that
input's checks.  That makes 2-fault failures (708 on [[15,7,3]] + K₈, brute force), and the MC shows p_out ∝ p².
Repairing it with "distrust checks that contain discarded inputs" cascades through the dense outer graphs.

**Cheapest plans** (T per output, full post-selection per stage, d = 4, Cliffords perfect):

| input p | target 1e-10 | 1e-15 | 1e-20 |
|---|---|---|---|
| 1e-3 | 69 ([[15,7]] → [[127,113]]) | 69 | 69 |
| 1e-4 | 8.1 ([[31,21]]) | 44.5 ([[31,21]] → [[127,113]]) | 44.5 |
| 1e-5 | 6.4 ([[63,51]]) | 9.6 ([[15,7]]) | 21.9 ([[64,26,7]] + 2080) |
| 1e-6 | 5.7 ([[127,113]]) | 5.7 | 16.3 ([[64,26,7]] + 2080) |

Order 7 only wins for very strict targets at very clean inputs.  Asymptotic γ is not the right figure of merit
here: the order-3 stages (γ 1.65–2.0) win on cost per output at realistic targets.

**End to end with cultivation (d = 4, physical p = 10⁻³)**, counted as front-end attempts per final output:

| target | raw injection | D3 cultivation (20 q) | 3→5 growth (56 q) |
|---|---|---|---|
| 1e-10 | 69 | 9.8 | 12.6 |
| 1e-15 | 69 | 47 | **12.6** ([[127,113,3]]) |
| 1e-20 | 69 | 47 | 36 ([[64,26,7]]) |

At physical p = 5×10⁻⁴: 1e-15 costs 7.3 (D3) / 8.2 (growth) attempts per output.  Growth inputs from the 10⁸-shot
run (`runs/cultivation/sweep_grow35_d4_long.json`: 9.1e-7 at 1e-3, 1.5e-7 at 5e-4).

Clifford noise: see §5e.  The earlier proxy (adding g·p_C to the fault rate) missed a first-order term and is superseded.

## 5e. Circuit-level validation (`meas_circuit.py`)
Explicit Clifford circuit for the zdsim Pauli-frame simulator.  T is replaced by identity: twirled T faults are Z^a,
which commute with T, so the Z statistics are exact.  X-type errors are treated conservatively.  Inputs are |+⟩.
Each check has a CSUM-network CSS encoder (an invertible L whose pivot rows are S and the normal logical basis N,
with the check's inputs as the slot pivots), the controlled transversal X¹ from a |+⟩ control read out in the X
basis, the inverse network, and syndrome readout (S-pivots in X, ancillas in Z).  Every detector is deterministic
(checked at d = 4, 8, 16).

* **The fault model is validated exactly.**  Injecting identical fault patterns into the circuit and the linear model
  gives identical verdicts on all ≤ 3-fault patterns of [[7,1,3]] (13 275 at d = 4, 161 315 at d = 8; 85 / 185 failures
  = the exact enumeration).  On [[15,7,3]] + K₈ it gives identical verdicts on all single faults, 120k random 2-
  and 3-fault patterns, and every planted inner-logical and input-triangle structure (632 failures at d = 4 = C₃;
  1184 at d = 8).
* **Clifford noise is first order** (DEP1/DEP2 on every gate, measurement and reset errors; p_in = p_T = 0;
  `runs/meas_distill/circuit_clifford_d4.json`):

| protocol | per-output error / p_C (Z only) | (Z or X, conservative) | Clifford locations per run |
|---|---|---|---|
| [[7,1,3]] | 0.4 | 2.9–3.1 | 58 |
| [[15,7,3]] + K₈ | 1.1 | 5.0–5.4 | 1979 |

  Faults in the final decoder land on outputs, so the last stage's logical Cliffords need **p_C ≲ target / 5**,
  as in qubit distillation.  Earlier stages only need p_C below their own output error.  Clifford locations (dense
  CSUM encoders, ~n² gates per block) outnumber T locations by about 7× and dominate rejection at equal rates.
* **Sparse encoders** (`sparse_encoder.py`).  The dense cost is a k×k mixing of the input qudits forced by the dense
  Gram–Schmidt normal basis.  A sparse normal basis is grown vector by vector (lightest odd-norm solution orthogonal
  to the vectors so far, by random information sets) with pivots chosen so the encoder matrix is nearly triangular,
  then synthesised in Markowitz order.  CSUMs per encoder (`runs/meas_distill/sparse_encoders.json`):

| inner code | dense (d = 4) | sparse d = 4 | sparse d = 8 | Σ basis weight (d = 4) | T per check |
|---|---|---|---|---|---|
| [[7,1,3]] | 19 | 12 | 12 | 3 | 14 |
| [[15,7,3]] | 111 | 68 | 68 | 41 | 30 |
| [[31,21,3]] | 565 | 273 | 340 | 158 | 62 |
| [[63,51,3]] | 2552 | 1142 | 1415 | 606 | 126 |
| [[64,26,7]] | 2010 | 1345 | 1580 | 479 | 128 |
| [[127,113,3]] | 10554 | 3998 | 5371 | 1902 | 254 |
| [[112,50,7]] | 5370 | 4721 | 5407 | 1647 | 224 |

  Validated like the dense circuit (determinism; all ≤ 3-fault patterns of [[7,1,3]] at d = 4, 8, 16 and
  66–89k random + planted patterns of [[15,7,3]] + K₈ at d = 4, 8, 16 agree exactly with the linear model run on
  the same basis).  The logical basis enters the failure count: [[15,7,3]] + K₈ C₃ = 656 vs 632 (d = 4), so the
  sparse basis trades ≈ 4 % in c for 1.6× fewer Cliffords.  Orthogonality is mod d, so the basis is per-d.
  An exact search on [[15,7,3]] shows a fully sparse normal basis needs weights up to 7 (≤ 5 gives at most 5 of 7);
  the sequential weights grow to ≈ n/2 for the last vectors, which is what remains of the k² term.
* Planner (`msd_planner.py`): per-output error c·p^W + 5·p_C and acceptance exp(−L p − L_C p_C).  With
  p_C ≤ target/6 the §5d costs are unchanged; with larger p_C the target cannot be reached.

## 6. Open / next
1. Done (§5d): concatenated order-3 stages; partial discard ruled out.
2. Done: w = 50 outer on 9500 inputs for [[112,50]] (γ 1.372, ≈ 9.6k qudits).  Open: smaller w = 26 outer (< 2080).
3. Done (§5e): circuit-level model; p_C ≲ target/5 for the final stage; sparse encoders 1.5–2.6× fewer CSUMs.  Open: a structured (recursive) sparse normal basis for RM(1,m) to remove the heavy tail.
4. Cultivation acceptance and p_L at d = 4 rest on few events (4 failures at 10⁻³); longer lenore sweeps would tighten §5d.
