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

**Other dimensions (exact counts, same circuits).**  The dominant term scales as (d−1)⁻⁶ at order 7 and the order-3
coefficients as roughly (d−1)⁻³, so every protocol improves with d at fixed p (the T-count, L and L_C do not change):

| protocol | c at d = 4 | d = 8 | d = 16 | source |
|---|---|---|---|---|
| O7 [[63,21,7]] + 1470 | 343.3 | 2.126 | 0.0220 | `count_63_21_n1470_d16.log`, `input_count_d8.log` |
| O7 [[64,26,7]] + 2080 | 360.0 | 2.229 | 0.0230 | `count_64_26_n2080_d{8,16}.log` |
| O3 [[7,1,3]] single check | 3.148 | 0.539 | 0.114 | §5d |
| O3 [[7,1,3]] two checks | 2.370 | 0.397 | 0.083 | `order3_d{8,16}.json` |
| O3 [[15,7,3]] + K₈ | 0.836 | 0.123 | 0.0237 | |
| O3 [[31,21,3]] + K₂₂ | 1.242 | 0.168 | 0.0297 | |
| O3 [[63,51,3]] + K₅₂ | 2.465 | 0.336 | 0.0582 | |
| O3 [[127,113,3]] + K₁₁₄ | 5.076 | 0.637 | 0.1020 | |

`msd_planner.stages_for_d(d)` builds the d = 8, 16 stage tables (`python -m primepower_msd.msd_planner 8`); the
optimal order-3 plans are acceptance-bound, so their T-per-output is the same at d = 8 as at d = 4 (68.8 at
p = 10⁻³ → 10⁻¹⁵) while reaching lower error.  Cultivation data (§5d, 5f) exist for d = 4 only.

## 5d. Noisy inputs: concatenation, not big outer codes (`msd_planner.py`)
Big order-7 rounds have L ≈ 16·N fault locations, so above p ≈ 10⁻⁵ their acceptance exp(−Lp) collapses.
**Order-3 stages** fix this: inner signed RM(1,m) [[2^m−1, 2^m−1−2m, 3]], outer = complete graph K_{k+1} (inputs = edges,
s = 2; order 3 needs only a simple graph).  Exact d = 4 leading-order coefficients (`runs/meas_distill/order3_d4.json`):

| stage | inputs | T/out | L | p_out per output |
|---|---|---|---|---|
| [[7,1,3]] single check (s = 1) | 1 | 15.0 | 15 | 3.15 p³ |
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
* **The first-order term cannot be removed, but it is 14× smaller than the conservative value.**  Single-fault
  census of the circuit (`cliff_census.py`: every Pauli at every DEP1/DEP2/XERR site, weighted by the number of
  bad outputs; reproduces the Monte Carlo coefficients) and an exact statevector simulation of the [[7,1,3]]
  gadget with real T gates on 8 qudits (`dense_gadget.py`: inputs = the +1 eigenstate of the logical check,
  a = Σ_f w_f P(accept | f)(1 − F_out | f)):

| gadget | d | frame model, Z-type a | frame model, any (conservative) | exact statevector |
|---|---|---|---|---|
| [[7,1,3]] single check | 4 | 0.357 | 2.92 | **0.359** |
| [[7,1,3]] single check | 8 | 0.147 | 2.41 | — |
| [[7,1,3]] two checks | 4 | 0.357 | 4.28 | — |
| [[15,7,3]] + K₈ | 4 | 1.02 | 5.22 | — (MC: 1.06–1.10 / 5.0–5.4) |
| [[15,7,3]] + K₈ | 8 | 0.42 | 3.86 | — |

  The exact value equals the Z-type one: every X-type fault counted by the frame model is either detected or
  projected back onto the magic state by the check (the check is a projective measurement onto the correct
  eigenspace, so a single input error is caught or undone), and only Z-type faults survive.  Those all sit in the
  **final decoder of an output's last check**: DEP2 faults on the output wire, or on block qudits whose Z
  propagates onto it through the remaining CSUMs — 91 of the 255 Paulis of one CSUM for [[7,1,3]], i.e. one
  gate's worth.  Any gadget that hands over a *decoded* output has this floor (its last gate can fail), so the term
  is second order only if outputs are consumed while still encoded, which is code switching.  The planner now
  carries a_C per stage (0.36 for [[7,1,3]], 1.02 for K₈, 5 where no census exists).
* Planner (`msd_planner.py`): per-output error c·p^W + a_C·p_C and acceptance exp(−L p − L_C p_C).  With
  p_C ≤ target/6 the §5d costs are unchanged; with larger p_C the target cannot be reached.

## 5f. Physical footprint: cultivation → one small round on a host code (`footprint.py`)
γ and T-per-output ignore where the Cliffords run.  The first-order term forces the last stage's Cliffords onto
logical qudits of a host code with p_L(D) ≲ target / a_C, and at ~2D² physical qudits per logical qudit a
2080-input round is ~10⁵–10⁶ physical qudits.  What fits under 10k is **cultivation (56 physical qudits, p_cult ≈
9×10⁻⁷ at p = 10⁻³) → one single-check [[7,1,3]] round on 8 logical qudits** (1 input + 6 block + 1 control):
its output 3.15 p_cult³ ≈ 2×10⁻¹⁸ is far below any target, so the host distance is set entirely by the Clifford
floor 0.36 p_L(D) < target.  Model: p_L(D) = 0.1 (p/p_th)^{(D+1)/2}, F(D) = 2D², escape of the cultivated state
into the host costs one p_L(D) (the host parameters were assumptions when this section was written; §5g measures
them: p_th ≈ 7.6×10⁻³, A ≈ 0.053, which is the "measured" column; §5h measures the escape itself: 4×10⁻⁴ per
kept attempt at p = 10⁻³, which supersedes the "one p_L(D)" charge and changes the attempts column, see §5h).

| p | target | p_th = 5×10⁻³, A = 0.1: D, physical | p_th = 10⁻², A = 0.1: D, physical | **measured Z₄ host (§5g)**: D, physical | cultivation attempts / output |
|---|---|---|---|---|---|
| 10⁻³ | 10⁻¹⁰ | 25, 10 056 | 17, 4 680 | **19  5 832** | 33 |
| 10⁻³ | 10⁻¹² | 31, 15 432 | 21, 7 112 | **23  8 520** | 33 |
| 10⁻³ | 10⁻¹⁵ | 39, 24 392 | 27, 11 720 | 31  15 432 | 33 |
| 5×10⁻⁴ | 10⁻¹⁰ | 17, 4 680 | 13, 2 760 | 15  3 656 | 22 |
| 5×10⁻⁴ | 10⁻¹² | 21, 7 112 | 17, 4 680 | 17  4 680 | 22 |
| 5×10⁻⁴ | 10⁻¹⁵ | 27, 11 720 | 21, 7 112 | 23  8 520 | 22 |

(With the conservative a_C = 5 the first row is D = 27, 11 720 qudits.)  Every other pipeline is larger: the
order-3 K₈ round needs 37 logical qudits, [[127,113,3]] + K₁₁₄ 6513, the order-7 rounds thousands.  The levers,
in order: the host code's threshold (a factor 2 in p_th is 4 steps of D, i.e. 2× in qudits), the physical error
rate, then cultivation acceptance (attempts only, 0.45 → more would cut the 33).  Cultivation alone cannot reach
10⁻¹⁰ at these p (9×10⁻⁷ after 3→5 growth, ~25× per growth step), so one distillation round is needed.

## 5g. The ring Z_{2^k} surface code: threshold and p_L(D) measured (`surface_zd.py`, `surface_sweep.py`, `surface_fit.py`)
**Code.** Signed Z_d rotated planar surface code, distance D: D² data qudits on a grid, checks on the faces
(checkerboard), weight-2 X checks on the left/right boundary and Z checks on the top/bottom.  Signs X^{(−1)^row} for
X checks and Z^{(−1)^col} for Z checks: any two overlapping checks share one row-pair or one column-pair, so the
signs cancel and the group is abelian over every Z_d (verified d = 4, 8; k = 1; logical Z = a column, logical X = a
row with constant sign, overlapping in one qudit, so the logical is a genuine Z_d qudit).  Syndrome extraction with
one ancilla per check: Z checks by CSUM(data → ancilla, w), X checks by F · CSUM(ancilla → data, w) · F†.  Schedule
X: NW SW NE SE, Z: NW NE SW SE — overlapping checks touch their shared pair in the same order (the interleaved
rounds measure the right operators; checked by gauge-randomised determinism), and each ancilla's hook error (its last
two legs) is perpendicular to the error chains that flip a logical (X chains run left–right, Z chains top–bottom).
The first schedule I tried (the N/Z pattern transposed) had the hooks parallel and gave a D = 3 code with fault
distance 2 (5.6% failure at p = 3×10⁻³ instead of 0.7%) — a useful reminder that the rotated code's schedule is tied
to the boundary assignment.

**Noise** (same convention as §5c/§5e: one "p" per location, uniform over the non-identity Paulis): DEP1 on the data
while the ancillas reset, XERR after every reset and before every measurement, DEP1 after F / F†, DEP2 after every
CSUM.  Memory experiment: D rounds, data prepared and measured in the memory basis, detectors = consecutive
differences (first-round and final-round detectors in the memory basis), failure = any nonzero Z_d shift of the
logical (an X² logical error counts).

**Decoder: 2-adic peeling MWPM.**  The detector-error model is built from the circuit by propagating one basis fault
(X¹ or Z¹ on one qudit) per noise location through the frame simulator; every non-identity Pauli at a location is a
Z_d combination of its basis faults.  Faults with ≤ 2 detectors in a sector are edges carrying the Z_d pattern
(c_a, c_b, o) of a unit flow (charges deposited at the two ends, logical shift); faults with more detectors are
decomposed into their single-generator components (stim's decompose_errors; no probability is dropped: hyperedge
share 0 for every D).  Where two fault classes share an edge with incompatible patterns (boundary degeneracies, e.g.
a data error on the logical column vs. a hook on the far side with the same single detector) the more probable
class fixes the pattern.  Decoding level j = 0..k−1 (d = 2^k): binary MWPM (PyMatching) on the detectors whose charge
has 2-adic valuation j, with weights from the probability of a valuation-j fault on each edge; each matched pair is
walked along a shortest path carrying the Z_d charge (flow w with w·c_a = −q at each step, logical shift += w·o),
leaving at the far end a residue divisible by 2^{j+1} for the next level.  After k levels every charge is zero
(no "leftover" in any run).  Sanity check at code capacity (data noise only, perfect syndrome, d = 4): the D =
7/11/15 curves cross at p ≈ 0.17, i.e. an X-error rate of ≈ 13.6% (12/15 of the depolarising Paulis carry an X
part), above the qubit MWPM value 10.3%, as expected for a qudit code.

**Results** (`runs/surface/zd_d4.jsonl`, `zd_d8.jsonl`; 10⁵ shots at p ≤ 2×10⁻³, 4×10⁴ above; p_fail after D rounds).

| d | basis | crossings (3,5) (5,7) (7,9) | fit A, p_th (sub-threshold, p_fail < 0.1, fails ≥ 5) |
|---|---|---|---|
| 4 | Z (X errors) | 6.6, 8.0, 8.3, (9,11) > 8 ×10⁻³ | 0.046, 7.8×10⁻³ |
| 4 | X (Z errors) | 6.3, 7.6, 8.1, (9,11) 7.7 ×10⁻³ | 0.049, 7.2×10⁻³ |
| 8 | Z | 6.9, 7.5, (7,9) > 8 ×10⁻³ | 0.050, 7.6×10⁻³ |
| 8 | X | 6.7, 7.5, (7,9) > 8 ×10⁻³ | 0.051, 7.0×10⁻³ |
| 16 | Z | 6.1, 8.2, 8.3 ×10⁻³ | 0.053, 7.5×10⁻³ |
| 16 | X | 6.3, 7.6, 7.8 ×10⁻³ | 0.066, 7.4×10⁻³ |
| 32 | Z | 6.3, 7.9, 8.1 ×10⁻³ | 0.053, 7.4×10⁻³ |
| 32 | X | 6.4, 7.6, 7.8 ×10⁻³ | 0.069, 7.4×10⁻³ |

**d = 16** (`zd_d16.jsonl`, 80 points, lenore, 2026-11-01; same shots): p_th unchanged; A drifts up 15–25% from d = 4
to d = 16 (refit on identical point sets: Z 0.046 / 0.050 / 0.053, X 0.053 / 0.058 / 0.066 for d = 4 / 8 / 16).  At the
D that matter the raw failure counts agree point by point: D = 7, 9 below 4×10⁻³ give d16/d4 ratios 0.97–1.33 (Poisson
scatter); the excess sits at D = 3 (up to 1.8× at 10⁻³) and near threshold — the extra valuation levels add a few short
uncorrectable chains.  A 20% larger A moves no footprint D, so the host is d-independent for 4 ≤ d ≤ 16.

**d = 32** (`zd_d32.jsonl`, 80 points, local, 2026-11-02; same shots; made feasible by the exact R1/R2 reduction of
`surface_zd.build_dem`, which enumerates only the d^|R1| Pauli rows of a location that touch sector detectors or the
observable instead of all d^m − 1 — identical DEMs at d = 4, 8, 20–60× faster): the drift has stopped.  A = 0.053 (Z) /
0.069 (X), p_th = 7.4×10⁻³ in both bases, equal to d = 16 within the fit scatter; the D = 7, 9 counts below 4×10⁻³ match
d = 16 point by point (Z, D = 9: 0 / 1 / 7 / 20 / 103 vs 0 / 1 / 8 / 18 / 87).  The host is d-independent for 4 ≤ d ≤ 32
up to a saturated +15–30% in A at the larger d, which moves no footprint D.

Selected points at d = 4 (Z memory): p = 10⁻³: D = 3 9.3×10⁻⁴, D = 5 1.7×10⁻⁴, D = 7 4×10⁻⁵, D = 9 < 10⁻⁵;
p = 3×10⁻³: 7.1×10⁻³, 4.0×10⁻³, 1.3×10⁻³, 4.5×10⁻⁴; p = 6×10⁻³: 2.7, 2.4, 1.9, 1.3 ×10⁻².  The crossings drift up
with D as usual; the asymptotic circuit-level threshold is **≈ 0.8–0.85%** for this noise model, and **d = 8 behaves
like d = 4** (the 2-adic decoder reduces every level to binary matching, and the 2-torsion sector dominates).  This is
on par with the qubit rotated surface code under the equivalent circuit model (≈ 0.7–1% depending on the idle/reset
details), so **the ring qudit surface code costs nothing in threshold**.  The fit p_L(D) = A (p/p_th)^{(D+1)/2}
with A ≈ 0.053, p_th ≈ 7.6×10⁻³ (all 129 sub-threshold points of d = 4 and 8, both bases; rms 0.2 in log) replaces
the §5f assumptions (0.1, 5×10⁻³): `footprint.HOST_Z4`, `python -m primepower_msd.footprint measured`:
**19  5 832 physical qudits** (D, count) for 10⁻¹⁰ and 23  8 520 for 10⁻¹² at p = 10⁻³ — the whole
cultivation + one-round pipeline is under 10k with margin.  D = 11 (d = 4) and D = 9 (d = 8) points at
5–8×10⁻³ confirm the crossings keep rising towards ≈ 0.8%.

Caveats: the decoder is a heuristic (a level-j graph with weights from the valuation-j faults only; a joint Z_d
decoder would do better, so the threshold is a lower bound); no idle noise during the CSUM steps; the fit's A and
p_th are effective parameters over 3 ≤ D ≤ 9 (the per-pair crossings say p_th is a little higher).

## 5h. The escape step measured: cultivated colour-code state → surface code (`graft_zd.py`, `graft_run.py`)
§5f charged the escape as one host logical error p_L(D); this section simulates it.  **Construction** (Gidney et al.
arXiv:2409.17595, circuits.zip on Zenodo 10.5281/zenodo.13777072, read as data by `graft_extract.py`): the d₁ = 5
colour code (19 qudits, 9 X + 9 Z plaquettes) is grafted into a D = 15 surface code — 229 data qudits = 19 + 105 fresh
|0⟩ + 105 fresh |+⟩, 228 checks (14 colour plaquettes kept, boundary plaquettes extended into the surface region,
interface checks of weight 2/3/5), measured for r₁ = 5 rounds, then the final matchable code (232 checks, the 8
colour-region Z plaquettes split into weight-2 pairs, k = 2 with the second logical commuting with ours) for r₂ = 5
rounds.  **Z_{2^k} transfer**: the check signs are a GF(2) linear system (an X/Z overlap of size 2m needs m sign
disagreements; the colour-code signature τ is fixed; kept/extended/refined checks equal their parents up to a unit;
the logicals commute) — it is consistent, every one of the 415/399 overlaps sums to exactly 0, so the signed graft is a
valid code for **all** d = 2^k (`graft_signs_solve.py`, `data/graft_d5_z4.json`; `tests/test_graft_zd.py` checks
d = 4, 8).  The logicals extend onto one fresh diagonal each (Z_L = Z^τ(19) × Z^v on 14 fresh |0⟩, X_L = X(19) × X^u
on 14 fresh |+⟩).  Syndrome extraction: one ancilla per check, CSUM schedules from the surface-code study (§5g;
weight-6 plaquettes in two conflict-free layers), the §5g noise model, idle noise on data during resets; detectors are
generated automatically (differences, first-round deterministic checks on fresh qudits, parent/child relations,
extensions) and verified deterministic by the frame simulator.  Decoder: §5g's 2-adic matching with the DEM built on
the kept shots only.

**Postselection** (the part that makes a matching decoder work on the colour region).  (i) Gidney's recipe: the
colour-region plaquettes of one colour class per sector are postselected during the r₁ grafted rounds (X: class 2,
Z: class 0 — chosen from the 9 + 1 combinations by hyperedge share and discard, `graft_run select`).  Without it
1–2% of the fault probability is hyperedges and the escape fails at the percent level.  (ii) Faults that fire a
postselected detector are removed from the error model (masking the detector columns instead left a boundary-edge
class with the wrong observable and a first-order failure equal to the "undetectable" share — found by a
single-fault test, `scratchpad/single_fault.py`).  (iii) **Decoder-certified first-order postselection**: at build
time every kept single-fault syndrome (52k patterns) is decoded; the patterns the decoder gets wrong are postselected
as well (`surface_zd.misdecoded_patterns`, `pattern_mask`).  The Z sector has none.  The X sector had 8 (mass
0.50 p per attempt): a Z-type fault on a boundary qudit of the colour code during its last round — that qudit sits in
only two colour plaquettes plus an interface check, so a single postselected colour class misses it, and the
resulting weight-3 hyperedge (observable 0) is matched with the wrong parity.  Before certification the X escape
failed at 1.2×10⁻³ with linear scaling in p (half of it this first-order hole); after it the extra discard is
5×10⁻⁴ per attempt and both sectors scale as p².

**Results** (`runs/surface/graft_variants.jsonl`, `graft_run variants`; 2×10⁵ attempts per sector and point, failures
per kept attempt, the discard is the joint one — both sectors' postselected detectors applied to every attempt, as the
error model assumes).  Three recipes: *one* = one colour class per sector (X: 2, Z: 0); *all* = all three classes in
both sectors; *one+flags* = *one* plus a flag qudit on every check of weight ≥ 5 during the grafted rounds (two
cancelling CSUMs bracketing the data CSUMs, flag measurement postselected, `build(flag_weight=5)`; 12 flags):

| p | one: X, Z, total, discard | all: X, Z, total, discard | one+flags: X, Z, total, discard |
|---|---|---|---|
| 5×10⁻⁴ | 2.0×10⁻⁵ (3), 0 (0), 2×10⁻⁵, 25% | 8×10⁻⁶ (1), 0, 8×10⁻⁶, 39% | 0, 8×10⁻⁶ (1), 8×10⁻⁶, 35% |
| 10⁻³ | 9.7×10⁻⁵ (11), 7.9×10⁻⁵ (9), 1.8×10⁻⁴, 43% | 4.0×10⁻⁵ (3), 2.7×10⁻⁵ (2), 6.6×10⁻⁵, 62% | 1.2×10⁻⁴ (10), 2.4×10⁻⁵ (2), 1.4×10⁻⁴, 58% |
| 2×10⁻³ | 5.8×10⁻⁴ (38), 2.6×10⁻⁴ (17), 8.5×10⁻⁴, 67% | 2.8×10⁻⁴ (8), 1.8×10⁻⁴ (5), 4.6×10⁻⁴, 86% | 5.1×10⁻⁴ (18), 3.7×10⁻⁴ (13), 8.8×10⁻⁴, 82% |

The final graft code's own 10-round memory and a plain D = 15 rotated code are 0 failures in 10⁵ at p ≤ 2×10⁻³
(10⁻⁵ and 2×10⁻⁵ at 3×10⁻³), so every escape failure comes from the five grafted rounds (the colour-code round, phase
0, is noiseless in this model).  Scaling of *one*: ×5 per ×2 in p at both steps, i.e. **fault distance 2** (≈ p^2.4);
*all* is ×5 and ×7 on few events.  **Flags do not help**: the X-sector rate (Z errors on data) is unchanged within
statistics, so the p² failures are not ancilla hooks on the weight-6 plaquettes but pairs of data errors in the colour
region whose plaquette syndrome (two unpostselected colours, a Z_d hyperedge) the matching decoder resolves wrongly;
the flags only add 15 points of discard (`flag_weight` stays as an option).  Postselecting all three classes removes
those pairs' escape routes (×3 lower error) for 19 more points of discard.  (The first version of these runs applied only
the decoded sector's postselection, giving X 3.0×10⁻⁴ + Z 1.1×10⁻⁴ with an understated 33% discard —
`graft_d4.jsonl` rows before the `variants` ones; superseded.)  For comparison Gidney's qubit d = 5 escape
(superdense extraction, full colour-region postselection, decoder-gap postselection) reaches ≈ 10⁻⁹ added error at
99% total discard.

**Decoder-certified second-order postselection** (`Escape.certify_pairs`, `surface_zd.misdecoded_pairs`; variants
*one+pairs*, *one+pairs-cap*, *one+pairs2*).  Since the p² failures are decoding failures, the first-order certification
extends to pairs: every unordered pair of single-fault patterns with at least one member firing a detector of a check
that touches the colour qudits is summed and decoded at build time, and the mis-decoded combined syndromes are
postselected like the first-order ones (cached in `runs/surface/cert_pairs_*.json`; 0.05–0.1 ms per pair).  Three
bookkeeping points that matter: (i) the first-order mis-decoded singles must be *excluded* from the candidates —
paired with anything they are mis-decoded again, and without the exclusion the X-sector "pair mass" read 417 p²
(colour × colour, 13M pairs) and 4 472 p² (colour × everything, 254M pairs, 90 min on 5 workers) when the genuine
numbers are **17.6 p²** (4.6k patterns) and ≈ 20 p²; the Z sector (no first-order bad singles) has 2.8 p² of
colour × colour pairs and 3.2 p² of colour × everything pairs (288M pairs, 97 min on 5 workers) — the surface and
interface regions add nothing at second order.  (ii) A postselected combined syndrome also discards every *single* fault whose syndrome
contains it (`pattern_mask` semantics): 186 of the X patterns carry 75 units of first-order mass (7.5% discard at
10⁻³, +2.5 points measured) for 1.3 units of second-order gain; `degenerate_cap` keeps such patterns in decreasing
order of gain/cost up to a first-order budget (2 units = 0.2% discard).  (iii) Genuine pair mass 17.6 + 2.8 ≈ 20 p²
is only a fifth of the measured p² coefficient (≈ 100 at 10⁻³ for *one*), so after first-order certification the
escape failures are mostly **third order and higher**: fault distance 3 with a prefactor ≈ 6×10⁴ p³ (5 105 and 5 688
colour-region patterns, 13.5k kept patterns in all).  Results (2×10⁵ attempts per sector; *one+pairs* = uncapped,
*one+pairs-cap* = cap 2):

| p | one+pairs-cap: X, Z, total, discard | one+pairs (uncapped) | one | all |
|---|---|---|---|---|
| 5×10⁻⁴ | 1.3×10⁻⁵ (2), 0 (0), 1.3×10⁻⁵, 24% | 1.4×10⁻⁵, 26% | 2×10⁻⁵, 25% | 8×10⁻⁶, 39% |
| 10⁻³ | 4.4×10⁻⁵ (5), 4.4×10⁻⁵ (5), 8.8×10⁻⁵, 43% | 8.2×10⁻⁵, 46% | 1.8×10⁻⁴, 43% | 6.6×10⁻⁵, 62% |
| 2×10⁻³ | 6.1×10⁻⁴ (39), 3.6×10⁻⁴ (23), 9.7×10⁻⁴, 68% | 7.5×10⁻⁴, 71% | 8.5×10⁻⁴, 67% | 4.6×10⁻⁴, 86% |

Colour × everything certification (*one+pairs2*, 10⁻³ only) gives the same numbers within statistics — X 4.6×10⁻⁵
(5) / capped 5.3×10⁻⁵ (6), Z 9×10⁻⁶ (1) / capped 4.4×10⁻⁵ (5) — as its genuine extra pair mass predicts, so the
cheap colour × colour scope is sufficient.  The capped recipe keeps the one-class discard (43% at 10⁻³) and halves the error (1.8×10⁻⁴ → 8.8×10⁻⁵; the X
improvement 11 → 5 events is more than the 18% the genuine pair mass predicts, so part of it is statistics); it is
the best recipe at both p and matches *all* in error with 19 points less discard.  Scaling 5×10⁻⁴ → 10⁻³ is ×6–7
(p^2.8), consistent with a third-order residual.  Beyond this the levers are a decoder that handles the colour-region
hyperedges (Z_d hypergraph / union-find on the colour plaquettes, or decoder-gap postselection as in the qubit escape)
— third-order certification (≈ 10¹¹ triples) is out of reach — or a larger colour code (d = 3 → 5 → 7 cultivation).

**Decoder options beyond matching** (`hyperdecoder.py`, variant *one+hyper*).  (a) *Hypergraph decoder*: BP + LSD
(`ldpc` 2.4) per 2-adic level on the exact single-fault patterns — hyperedges included, nothing decomposed; level j
decodes ((σ ≫ j) & 1) over the distinct level-j bit patterns (6 056 columns per level at d = 4), subtracts the chosen
patterns' exact Z_d values (the key whose values match the syndrome), and offers the next level one "unit-flip" column
per chosen pattern (Z vs Z⁻¹ ambiguity, prior 0.4, set per shot with `update_channel_probs`).  2.7 ms per shot
(matching: 0.1 ms).  It is better at first order — one mis-decoded colour-region pattern (mass 0.03 p) against the
matcher's eight (0.33 p) — but **worse in bulk**: X 2.0×10⁻⁴ (23), Z 7.9×10⁻⁵ (9) per kept attempt at p = 10⁻³ against
4.4×10⁻⁵ / 4.4×10⁻⁵ for the certified matcher (X 4.6×10⁻⁵ / Z 1.3×10⁻⁵ at 5×10⁻⁴; 1.4×10⁻³ / 7.5×10⁻⁴ at 2×10⁻³),
the usual BP+OSD/LSD deficit against MWPM on circuit-level surface-code noise.  Tuning does not close it: on one
12.9k-shot kept sample at 2×10⁻³ (X) the matcher fails 7 times, BP+LSD 13–19 for every setting tried (LSD order 0 /
10 with LSD_CS or LSD_E, min-sum with scaling 0.625 and 100 iterations = best at 13, product-sum 19; `hyper_tune.log`).  (b) *Decoder-gap
postselection*: with the observable as an extra parity row, decoding with the bit forced to 0 and to 1 gives the
complementary gap directly in the BP framework, but forcing the unlikely class defeats BP and falls through to LSD on
the full matrix every time (157 ms per shot, 60× the plain decode), and a pure-BP gap converges for only ≈ 40% of
shots in either class (uninformative for ¾ of them).  A matching-based gap needs a logical representative on a patch
boundary so that every observable-flipping edge is a boundary edge; not built.  Net: the certified 2-adic matcher
(*one+pairs-cap*) stays the decoder of record; a gap would need the matching route or a cheaper complementary decode.

**Complementary-gap postselection with the 2-adic matcher** (`PeelDecoder(gap=True)`, `boundary_potential`; variants
*one+gap*, *one+pairs-cap+gap*).  The qubit trick (logical representative on a patch boundary so that every
observable-flipping edge is a boundary edge) does not transfer directly — signed Z_d hooks flip the observable on
interior edges — so the class is forced differently: a parity potential x on the detectors with o + x_a + x_b ≡ 0 on
every interior edge of the level graph (a BFS; zero inconsistent cycles at both levels) absorbs the interior flips,
the boundary is split into two regular nodes B₀/B₁ (a boundary edge goes to B_{x_a + o}), and a correction's logical
class is the parity of its edge-set degree at B₁.  Decoding the level syndrome with (B₀, B₁) excited as (|σ|, 0) and
(|σ|+1, 1) mod 2 gives the two classes exactly; the lighter wins, the weight difference is the level's gap and the
shot's gap is the minimum over levels; charges are then walked along the edge set the matcher returned (its strings may
end on either boundary).  0.7 ms per shot, identical decisions to the plain matcher on every shot
(`tests/test_surface_zd.py::test_gap_decoder_matches_plain_decoder`).  Failing shots have small gaps: at p = 2×10⁻³
the six failures of a 6.6k-shot sample had gaps 0.3–8 against a median of 21 and a 5th percentile of 12, and in the
2×10⁵-attempt runs (both bases, 5×10⁻⁴–2×10⁻³) discarding the 5% lowest-gap kept shots removes every failure but one
(at 2×10⁻³).  With a fine threshold grid (4×10⁵ attempts per sector at
p = 2×10⁻³, *one+pairs-cap+gap*, `graft_variants_gap_2e-3.log`; total discard without the gap cut 67.7%):

| gap cut | X: fails of 58, p_fail, total discard | Z: fails of 56, p_fail, total discard |
|---|---|---|
| none | 58, 4.5×10⁻⁴, 67.7% | 56, 4.3×10⁻⁴, 67.6% |
| ≥ 3 | 14, 1.1×10⁻⁴, 67.7% | 14, 1.1×10⁻⁴, 67.6% |
| ≥ 5 | 7, 5.4×10⁻⁵, 67.8% | 8, 6.2×10⁻⁵, 67.7% |
| ≥ 8 | 3, 2.4×10⁻⁵, 68.1% | 1, 7.8×10⁻⁶, 67.9% |
| ≥ 10 | 0, < 10⁻⁵, 68.5% | 0, < 10⁻⁵, 68.3% |
| ≥ 15 | 0, —, 71.6% | 0, —, 70.6% |

The failing shots' gaps are 0–9.5 (median 1.5) while the kept shots' median gap is ≈ 18: **a cut at 10 removes every
observed failure for 0.8 points of total discard (2.5% of the kept shots)**, a ≥ 50× reduction of the escape error at
2×10⁻³ for a 1.3% increase in attempts.  At p = 10⁻³ (10⁶ attempts per sector, `graft_variants_gap_1e6.log`; total discard without the cut 43.1%):

| gap cut | X: fails of 38, p_fail, total discard | Z: fails of 14, p_fail, total discard | X + Z |
|---|---|---|---|
| none | 38, 6.7×10⁻⁵, 43.1% | 14, 2.5×10⁻⁵, 43.1% | 9.1×10⁻⁵ |
| ≥ 5 | 4, 7.0×10⁻⁶, 43.2% | 3, 5.3×10⁻⁶, 43.1% | 1.2×10⁻⁵ |
| ≥ 10 | 2, 3.5×10⁻⁶, 43.4% | 0, —, 43.3% | 3.5×10⁻⁶ |
| ≥ 13 | 0, < 4×10⁻⁶ (90% CL), 43.7% | 0, < 4×10⁻⁶, 43.6% | < 8×10⁻⁶ |
| ≥ 15 | 0, —, 44.3% | 0, —, 44.0% | |

(failing shots' gaps 0–12.7, median 1.6; kept shots' median ≈ 24).  **Recipe of record: *one+pairs-cap+gap* with a gap
cut at 13 — escape error < 8×10⁻⁶ per kept attempt (90% CL, no failure in 1.14M kept attempts; ≥ 25× below the
uncertified-gap 9.1×10⁻⁵) for 0.6 points of total discard (43.7%, i.e. 1% more attempts).**  This is what the qubit
escape gets from the complementary gap (Gidney et al.: 6×10⁻¹⁰ → 2×10⁻⁹ end-to-end for the d = 5 state) scaled to our
noisier single-ancilla circuit and smaller statistics; the escape is no longer the pipeline's bottleneck.



**d = 8** (`graft_run variants … 8`, `runs/surface/graft_variants_d8.jsonl`, 2×10⁵ attempts per sector and point, lenore):
the same graft, decoder and recipes run unchanged at d = 8 (Escape(8) builds in 62 min single-threaded, 6 142 edges, no
undecomposed hyperedge mass, 30 postselected detectors) and reproduce the d = 4 numbers:


**d = 16** (`graft_run variants … 16`, `runs/surface/graft_variants_d16.jsonl`, 2×10⁵ attempts per sector and point, local,
2026-11-02; 15–20 min per task once the DEM is built): reproduces d = 8.  One class at 10⁻³: 7 (X) + 8 (Z) failures in
1.08×10⁵ kept per sector, 6.9×10⁻⁵ per kept attempt summed over sectors (d = 8: 7.7×10⁻⁵), discard 0.459 (d = 8: 0.451); at
2×10⁻³ 31 + 28 in 5.9×10⁴ (d = 8: 45 + 24).  Gap recipe: no failure at cut 13 or 15 in 1.08×10⁵ (10⁻³) and 5.8×10⁴
(2×10⁻³) kept per sector; failing gaps ≤ 9.0 (10⁻³) and ≤ 8.1 (2×10⁻³).  Discard grows by under 1 point from d = 8 to d = 16.
**d = 32** (`graft_variants_d32.jsonl`, 2×10⁵ attempts per sector and point, lenore, finished 2026-11-02): reproduces
d = 8 / 16.  One class at 10⁻³: 9 (X) + 5 (Z) in 1.07×10⁵ kept per sector, 1.3×10⁻⁴ summed, discard 0.464; at 2×10⁻³
38 + 27 in 5.7×10⁴ (1.1×10⁻³), discard 0.713.  Gap recipe: no failure at cut 13 or 15 at either p; failing gaps ≤ 7.8
(10⁻³) and ≤ 7.0 (2×10⁻³); discard at cut 15 0.473 (10⁻³), 0.74 (2×10⁻³).  The escape is d-independent through d = 32.
(The pair-certified recipe, *one+pairs\**, is not run beyond d = 8: certification decodes every pair of colour-region
single-fault patterns, ~10¹⁰ pairs per sector at d = 16.)
| p | recipe | X sector | Z sector | total per kept attempt | joint discard | d = 4 (§5h) |
|---|---|---|---|---|---|---|
| 10⁻³ | one class | 1.0×10⁻⁴ (11) | 5.5×10⁻⁵ (6) | 1.6×10⁻⁴ | 0.451 | 1.8×10⁻⁴, 0.43 |
| 10⁻³ | one class + gap cut 13 | 0 / 1.09×10⁵ kept | 0 / 1.09×10⁵ kept | < 4×10⁻⁵ (90% CL) | 0.454 | < 8×10⁻⁶ (1.14M kept), 0.437 |
| 2×10⁻³ | one class | 7.4×10⁻⁴ (45) | 4.0×10⁻⁴ (24) | 1.1×10⁻³ | 0.698 | — |
| 2×10⁻³ | one class + gap cut 13 | 0 / 5.7×10⁴ kept | 0 / 5.7×10⁴ kept | < 8×10⁻⁵ (90% CL) | 0.715 | — |

Every failing shot at 2×10⁻³ had gap ≤ 8.1 (≤ 10.5 at 10⁻³ without the cut), so the cut at 13 removes all of them at the
same discard as d = 4; the escape, like the host code (§5g) and the cultivation check, is d-independent in practice.

**Footprint** (`footprint.ESCAPE_Z4` = *one*, `ESCAPE_Z4_ALL` = *all*, `ESCAPE_Z4_PAIRS` = *one+pairs-cap*,
`ESCAPE_Z4_GAP` = *one+pairs-cap+gap* cut 13 (cut 15 since the margin run below); `python -m primepower_msd.footprint escape` / `escape-all` /
`escape-pairs` / `escape-gap`): without the gap cut the escape error adds to p_cult (9×10⁻⁷) and dominates it — the
cultivated state leaves the escape at 8×10⁻⁵–1.8×10⁻⁴, worse than the 2×10⁻⁵ of a d = 3 cultivation, and the escape,
not the Clifford floor, sets the pipeline; with the gap cut (< 8×10⁻⁶, 43.7% discard) the Clifford floor is again the
binding constraint and one [[7,1,3]] round reaches 10⁻¹² at p = 10⁻³ (D = 23, 8 520 qudits, 59 attempts/output):

| p | target | pipeline | physical qudits | cultivation attempts / output: one+pairs, one, all (§5f free escape) |
|---|---|---|---|---|
| 10⁻³ | 10⁻¹⁰ | one [[7,1,3]] round, D = 19 | 5 832 | 59, 58, 58, 87 (33) |
| 10⁻³ | 10⁻¹² | one round with the gap recipe (two otherwise), D = 23 | 8 520 | 59, 877, 876, 1 311 (33) |
| 5×10⁻⁴ | 10⁻¹⁰ | one [[7,1,3]] round, D = 15 | 3 656 | 30, 30, 30, 37 (22) |
| 5×10⁻⁴ | 10⁻¹² | one [[7,1,3]] round, D = 17 | 4 680 | 30, 30, 30, 37 (22) |

The qudit count is unchanged (the 8 logical qudits are reused round after round); the escape costs a factor 1.4–1.9 in
time through its discard.  Without the gap cut the 10⁻¹² target at p = 10⁻³ needs a second round
(3.15 × (9×10⁻⁵)³ = 2.3×10⁻¹²); with it one round suffices.  The gap recipe at 5×10⁻⁴ (10⁶ attempts per sector, `graft_variants_gap_5e-4.log`): before the cut 6 + 5 failures (1.5×10⁻⁵ per kept attempt, 24.5% discard), after the cut at 13 no failure in 0.75M kept attempts per sector (< 6×10⁻⁶ total, 90% CL) for 0.4 points of discard (24.9%) — `ESCAPE_Z4_GAP[5e-4]`; with the measured CSUM cost (§5i) the 10⁻¹² row moves to D = 19, 5 832 qudits.

**Floor run at 10⁻³** (lenore, 2026-10-31, `graft_variants_gap_floor.jsonl`, 3 replicas × 4×10⁶ attempts per sector): before the gap cut 389 (X) + 270 (Z) failures in 6.83M kept per sector, 9.6×10⁻⁵ per kept attempt — the one+pairs-cap value of 9.1×10⁻⁵ reproduced; at the cut 13 no failure in 6.77M kept per sector (discard 0.436).  Together with the first 1.14M kept: **no failure in 7.9M kept attempts per sector, < 2.9×10⁻⁷ per sector and < 6×10⁻⁷ in total (90% CL)** — the escape adds less error than the cultivation itself (9×10⁻⁷).  The cut has little margin: the failing gaps reach 12.x in one replica (2 failures per 2.25M kept at cut 12, the lowest zero-failure cut per replica is 10, 11, 11, 11, 12, 13), so the true floor at cut 13 is likely within a factor of a few of the bound; a cut at 15 would cost ~1 more point of discard.  Caveats: no idle noise during CSUM layers; the final-code and D = 15 baselines (no postselection, `graft_d4.jsonl`) had no failure in 10⁵ shots at p ≤ 2×10⁻³, so the escape error is the graft's own and postselecting the baselines would add nothing; the 5×10⁻⁴ points rest on 1–3 events each.

**Margin run at 10⁻³** (lenore, 2026-11-02, `graft_variants_gap_margin.log` + `graft_variants_gap_margin.jsonl`, 6 more replicas ×
4×10⁶ attempts per sector; pooled with the floor run, 20.5M kept per sector).  Before the cut 872 (Z) + 1 156 (X) failures
(9.9×10⁻⁵ per kept attempt; the pairs-recipe rate reproduced).  The bound above was optimistic by the margin it warned of: the tail of failing gaps reaches 13–14.

| gap cut | total discard (Z / X) | Z fails | X fails | p_fail per kept attempt (X + Z) |
|---|---|---|---|---|
| 12 | 0.434 / 0.435 | 5 | 19 | 1.2×10⁻⁶ |
| 13 | 0.436 / 0.437 | 1 | 7 | 3.9×10⁻⁷ |
| 14 | 0.438 / 0.438 | 1 | 2 | 1.5×10⁻⁷ |
| 15 | 0.440 / 0.443 | 0 | 1 | 4.9×10⁻⁸ (90% CL < 3×10⁻⁷) |
| 16–18 | 0.443–0.458 | 0 | 1 | 4.9×10⁻⁸ |

Cut 13 already keeps the escape below the cultivation error (3.9×10⁻⁷ vs 9×10⁻⁷).  **The recipe of record now cuts at 15**: the
escape drops to ~5×10⁻⁸ (one event in 20.5M; it survives to cut 18, so it is the decoder's, not the cut's) for 0.5 points of
discard.  `ESCAPE_Z4_GAP` = (3×10⁻⁷, 0.441) at 10⁻³ (90% CL, both sectors), (6×10⁻⁶, 0.250) at 5×10⁻⁴, where no failure was seen at either cut.  No
footprint D moves: 10⁻¹⁰ at 10⁻³ is still D = 21, 7 112 (census), and 1e-12 is still D = 25, 10 056 (sparse census).  Attempts per output go from 58.9 to 59.4.

## 5i. The logical Clifford layer measured: transversal CSUM between surface-code patches (`logical_csum.py`)

The footprint charged every Clifford fault location (the a_C = 0.36 and L_C counts are CSUM gates) one host logical
error p_L(D) — a single patch idling for D rounds — without measuring a logical gate.  `csum_circuit(D, d, rounds, rounds,
p, basis)` builds two identical signed Z_d rotated patches, D rounds of CSUM-extracted memory, one transversal
CSUM(control → target, weight w = 1) layer (two-qudit depolarizing noise on every gate, idle noise on both patches), D
more rounds, and a transversal readout.  Transversal CSUM is a logical CSUM: X-check_c → X-check_c X-check_t^w and
Z-check_t → Z-check_t Z-check_c^{−w}, so the detectors across the layer are m_t^Z(r) − m_t^Z(r−1) − w·m_c^Z(r−1) and
m_c^X(r) − m_c^X(r−1) + w·m_t^X(r−1) (verified deterministic for d = 4 and 8, both bases); the logical action
X̄_c → X̄_c X̄_t^w, Z̄_t → Z̄_t Z̄_c^{−w} makes |0,0⟩ and |+,+⟩ deterministic for both logicals, and a run fails when
*either* logical readout is wrong.  Pre-layer faults on the source patch are copied: an X error on the control before
the layer becomes X_c X_t^w, a weight-4 hyperedge (control pair + target "shadow" pair) that `build_dem` decomposes.

**Two decoders.**  *Independent*: each patch matched on its own with the decomposed DEM.  *Pauli-frame* (`SequentialDecoder`):
X errors flow control → target and Z errors target → control, so the source patch is decoded first with a DEM restricted
to its own detectors, the traversed edges (with their Z_d multipliers, `PeelDecoder.trace`) are mapped through a shadow
table (sink-detector pattern and sink-observable coefficient per unit multiplier of each source edge, from the basis
propagation of the full circuit, majority vote over faults sharing an edge) onto the sink syndrome and observable, and the
sink is then decoded with a DEM of its own noise only.  Every single fault of the D = 3 circuit is corrected in both bases
(`tests/test_logical_csum.py`).  The independent decoder mis-decodes 6 single faults per sector at D = 3: ancilla faults
equivalent to a *source* measurement error at the pre-layer round, which the cross-layer detector turns into one lone
sink defect that per-patch matching cannot pair.  This is a first-order floor, not a distance effect — the independent
CSUM fails at 0.8–2.4% at p = 10⁻³ (growing with D in proportion to the ancilla count, about 24 p at D = 9), 10³× the
memory at D = 9 — and the independent numbers below are a diagnostic, not a protocol.

**Sweep** (`python -m primepower_msd.logical_csum sweep 4 …`, `runs/surface/csum_d4.jsonl`): d = 4, D = 3, 5, 7, 9,
p = 1, 2, 3, 4, 6 ×10⁻³, 10⁵ shots (4×10⁴ at p ≥ 3×10⁻³), both bases; *base* = the same two patches and 2D rounds with
no gate.  Fits of p_fail = A (p/p_th)^{(D+1)/2} over the 17 (16) sub-threshold points per basis (≥ 5 events, p_fail < 0.1):

| experiment | Z basis (X-type errors): A, p_th, rms | X basis (Z-type errors): A, p_th, rms |
|---|---|---|
| two-patch memory, 2D rounds | 0.190, 7.34×10⁻³, 0.14 | 0.193, 6.95×10⁻³, 0.16 |
| logical CSUM, Pauli-frame decoder | 0.219, 7.34×10⁻³, 0.18 | 0.214, 6.90×10⁻³, 0.34 |
| logical CSUM, Pauli-frame decoder, **d = 8** (lenore, `csum_d8.jsonl`) | 0.250, 7.29×10⁻³, 0.15 | 0.216, 6.49×10⁻³, 0.17 |
| logical CSUM, independent decoder | no scaling (fitted "p_th" 1.1×10⁻³, rms 0.9): first-order floor | same |

Failure counts per 10⁵ shots at p = 10⁻³ / 2×10⁻³ (either logical wrong), Pauli-frame CSUM vs memory:

| D | Z: CSUM | Z: memory | X: CSUM | X: memory |
|---|---|---|---|---|
| 3 | 431 / 1 579 | 357 / 1 397 | 436 / 1 731 | 416 / 1 580 |
| 5 | 81 / 508 | 65 / 521 | 93 / 625 | 75 / 571 |
| 7 | 10 / 151 | 7 / 112 | 9 / 162 | 8 / 156 |
| 9 | 2 / 38 | 2 / 29 | 5 / 36 | 2 / 27 |

The logical CSUM costs 1.1–1.2× the idle two-patch memory at the same threshold; the gate's excess over idling is
A ≈ 0.02–0.03, i.e. about half a single-patch D-round memory error (A = 0.053, §5g; the two-patch 2D-round memory is 3.6×
that, as expected).  Control and target fail at comparable rates (the sink slightly more often in the Z basis, the source
in the X basis).  The two-patch sweep gives p_th 7.0–7.3×10⁻³, a little below the single-patch fit of record (7.55×10⁻³).

**Footprint** (`footprint.CSUM_Z4`, `BASE2_Z4`; `python -m primepower_msd.footprint escape-gap csum` / `csum-excess`):
*csum* charges each Clifford location the whole measured two-patch CSUM experiment error (conservative: it includes 2D
rounds of idling of both patches per gate); *csum-excess* charges one p_L(D) as before plus the gate's measured excess
over idling.  At p = 10⁻³ with the gap escape recipe, one [[7,1,3]] round:

| target | Clifford cost model | D | physical qudits |
|---|---|---|---|
| 10⁻¹⁰ | one p_L(D) per CSUM (§5h) | 19 | 5 832 |
| 10⁻¹⁰ | measured excess (*csum-excess*) | 19 | 5 832 |
| 10⁻¹⁰ | whole CSUM experiment (*csum*) | 21 | 7 112 |
| 10⁻¹² | one p_L(D) per CSUM (§5h) | 23 | 8 520 |
| 10⁻¹² | measured excess / whole experiment | 25 | 10 056 |

The measured Clifford layer moves the 10⁻¹⁰ pipeline by at most one distance step and the 10⁻¹² pipeline by one step to
D = 25, which sits at the 10k budget; the sensitivity is the (D+1)/2 power of p_th (7.1 vs 7.55 ×10⁻³ is a factor 2 at
D = 23), so the 10⁻¹² row should be read as 8.5–10k.  Caveats: the model still has no census of idle rounds between
CSUM layers (the *csum* variant over-charges them, *csum-excess* omits them); the d = 8 sweep (same grid, 40 points) gives the same A and a
threshold 0–6% lower, within the scatter between bases; the shadow table is a majority vote over faults sharing an edge (no
ambiguous edges at D = 3; `SequentialDecoder.ambiguous` reports them).

## 5j. Host-code census of the Clifford layer: the depolarizing a_C under-counts (`host_census.py`)

§5e's a_C = 0.357 assumes every Pauli of a failing logical gate is equally likely (depolarizing p_C), and §5f–5i charged
each Clifford location one host logical error p_L(D) or the whole measured CSUM experiment.  Neither is the host's noise:
a surface-code logical CSUM fails mostly by a *single-wire* logical Pauli, a patch that waits fails per round, a Fourier
gate costs a patch rotation (or nothing, if the frame is tracked), and a readout fails like one round.  `host_census.py`
schedules the [[7,1,3]] s = 1 gadget circuit of `meas_circuit.build` into rounds (CSUM = 1, F/FD = τ_F, R/M = 1; 45 CSUMs,
depth 31), injects every single-wire logical Z^a and X^a at every CSUM wire, idle patch-round, Fourier gate, reset and
readout, and counts the undetected bad outputs per unit of per-location error — the *exposure* N of each location class.
The depolarizing fault set reproduces 0.357 (Z-type) and 2.92 (any) exactly.  Z-type exposures per output at d = 4:

| location class | exposure (weight) | undetected bad outputs, any type | Z-type |
|---|---|---|---|
| CSUM wires (90) | 90 | 4.67 (+ 5.67 X-type) | **4.67** |
| idle patch-rounds (130 at τ_F = 0) | 130 | 5.33 (+ 19.0 X-type) | **5.33** |
| Fourier gates (τ_F = 15 rounds each) | 135 | 0 (+ 30 X-type) | 0 |
| resets (8), readouts (7) | 15 | 2.0 (+0) | 0 |

Only Z-type gate and idle faults survive (the check projects X-type faults away, §5e), Fourier rotations and readouts
contribute nothing of Z type, and a τ_F = 15 rotation per Fourier gate adds 14 idle rounds of exposure (the output waits).
The floor per output is therefore

  floor(p, D) = 4.67 q_gate(D) + 5.33 q_round(D),   q_round = Z-type logical error per patch per round = p_base2,X(D) / 4D,
  q_gate = Z-type error per wire per logical CSUM = (p_csum,X − p_base2,X)/2 + q_round  (X-basis fits of §5i),

against a_C p_L(D) = 0.357 × 0.053 (p/7.55×10⁻³)^{(D+1)/2}.  Why the old coefficient is 10× too small: a_C counts one
gate's worth of *depolarizing* failures, 91 bad of 255 Paulis; the host's failures are concentrated on exactly the
single-wire Z-type Paulis that survive, so per unit of gate-failure probability about 2.3 outputs go bad instead of 0.36,
and the idle rounds, never charged before, add as much again.  At p = 10⁻³:

| D | a_C p_L(D) (§5f–h) | *csum* (§5i, conservative) | census, CSUM excess from the fit | census, CSUM = one extra round |
|---|---|---|---|---|
| 19 | 3.1×10⁻¹¹ | 2.4×10⁻¹⁰ | 4.3×10⁻¹⁰ | 1.3×10⁻¹⁰ |
| 21 | 4.2×10⁻¹² | 3.3×10⁻¹¹ | 6.2×10⁻¹¹ | 1.9×10⁻¹¹ |
| 23 | 5.5×10⁻¹³ | 4.7×10⁻¹² | 9.1×10⁻¹² | 2.8×10⁻¹² |
| 25 | 7.3×10⁻¹⁴ | 6.6×10⁻¹³ | 1.3×10⁻¹² | 4.0×10⁻¹³ |
| 27 | 9.7×10⁻¹⁵ | 9.3×10⁻¹⁴ | 2.0×10⁻¹³ | 6.0×10⁻¹⁴ |

The two census columns bracket the one uncertain input: the CSUM's excess over idling.  The fits give the CSUM experiment
1.11× the two-patch memory, i.e. an excess of ≈ 0.4 D patch-rounds per gate (10 rounds at D = 23), but the per-D ratios
(1.10, 1.09, 1.04, 1.33 in the X basis at D = 3–9) cannot distinguish that from a constant one extra round per wire; the
"one extra round" column is the optimistic reading.  **Fixed-D measurement** (`csum_excess.log`, 10⁶ shots per point, 4 workers,
2.4 h): CSUM over no-gate memory, both decoded as in §5i,

| D | basis | p | CSUM / memory failures | ratio | excess in patch-rounds (4D × (ratio − 1)) | one extra round per patch |
|---|---|---|---|---|---|---|
| 5 | X | 2×10⁻³ | 6 156 / 5 604 | 1.099 ± 0.020 | 2.0 ± 0.4 | 2 |
| 5 | Z | 2×10⁻³ | 5 315 / 4 813 | 1.104 ± 0.022 | 2.1 ± 0.4 | 2 |
| 5 | X, Z | 10⁻³ | 787 / 785, 686 / 655 | 1.00 ± 0.05, 1.05 ± 0.06 | 0.1 ± 1.0, 0.9 ± 1.1 | 2 |
| 7 | X | 2×10⁻³ | 1 699 / 1 496 | 1.136 ± 0.040 | 3.8 ± 1.1 | 2 |
| 7 | X | 3×10⁻³ | 8 532 / 7 685 | 1.110 ± 0.017 | 3.1 ± 0.5 | 2 |
| 7 | Z | 2×10⁻³ | 1 429 / 1 303 | 1.097 ± 0.042 | 2.7 ± 1.2 | 2 |
| 7 | Z | 3×10⁻³ | 7 056 / 6 197 | 1.139 ± 0.020 | 3.9 ± 0.6 | 2 |

The excess grows with D — 2.0 ± 0.3 patch-rounds at D = 5, 3.4 ± 0.4 at D = 7 — as ≈ 0.45 D, which is the fitted
0.11 × 4D; a constant one extra round per patch is excluded at D = 7 (3.5σ).  The copying mechanism explains the
scaling: the layer propagates the control's whole *uncorrected* error history (∝ D rounds) into the target as
correlated shadows, which the Pauli-frame decoder removes only when the source decode is right.  So the **"CSUM excess
from the fit" column is the measured one**; the "one round" column is dropped.  (At D = 5 the 10⁻³ points show less
excess than at 2×10⁻³, 0.5 ± 0.7 vs 2.0 ± 0.3 rounds, hinting at a higher power of p in the excess; resolved below.)
**p-dependence** (`csum_excess_p.log`, 10⁶ shots per point, lenore, finished 2026-11-02; ratio, excess in patch-rounds):

| D | p = 1.5×10⁻³ (X / Z) | 3×10⁻³ (X / Z) | 4×10⁻³ (X / Z) |
|---|---|---|---|
| 5 | 1.14 / 1.13 (2.8 / 2.6 ± 0.7) | 1.11 / 1.14 (2.3 / 2.8 ± 0.3) | 1.14 / 1.14 (2.8 / 2.9 ± 0.2) |
| 7 | 1.13 / 1.06 (3.6 / 1.8 ± 2.0) | (above) | 1.14 / 1.15 (4.0 / 4.3 ± 0.3) |
| 9 | — | 1.14 / 1.13 (5.1 / 4.6 ± 1.2) | 1.14 / 1.20 (4.9 / 7.0 ± 0.6) |

The ratio is p-independent from 1.5×10⁻³ to 4×10⁻³ (1.10–1.15; the D = 5, 10⁻³ points are a 2σ low), so the excess is
the same power of p as the memory and grows ≈ 0.55 D patch-rounds (ratio ≈ 1.13 against the fitted 1.11).  That is a 2%
change in the CSUM error per location: no footprint change.
**Footprint** (`python -m primepower_msd.footprint escape-gap census` /
`census-round`), p = 10⁻³, one [[7,1,3]] round, gap escape:

| target | a_C p_L (§5h) | census (fit) | census (one round) |
|---|---|---|---|
| 10⁻¹⁰ | D = 19, 5 832 | D = 21, 7 112 | D = 21, 7 112 |
| 10⁻¹² | D = 23, 8 520 | D = 27, 11 720 | D = 25, 10 056 |
| 10⁻¹⁰ at 5×10⁻⁴ | D = 15, 3 656 | D = 15, 3 656 | D = 15, 3 656 |
| 10⁻¹² at 5×10⁻⁴ | D = 17, 4 680 | D = 19, 5 832 | D = 19, 5 832 |

The 10⁻¹⁰ pipeline at 10⁻³ is 7 112 qudits with the Clifford layer measured and censused (one distance step above the
§5h figure); 10⁻¹² at 10⁻³ needs D = 27 (11 720, the fit column, now measured), over the 10k budget, and 5.8k at 5×10⁻⁴.  Caveats: the census
takes the logical CSUM's failures as single-wire (correlated two-wire logical failures are not weighted; the Pauli-frame
decoder removes most of the copying, §5i); τ_F = 0 (Fourier gates by frame tracking); the exposures are exact for the
dense encoder of the s = 1 gadget (superseded by the sparse census below); d = 4 exposures.

**Sparse-encoder census** (`host_census.gadget(d, sparse=True)`, `scratchpad/sparse_census_scan.py` →
`runs/meas_distill/sparse_census_d{4,8}.jsonl`, 200 / 100 sparse normal bases and encoder networks, 2026-11-01).  The
exposures depend on *which* single-wire faults the check misses, not on the CSUM count: the 29-CSUM networks (the
12-CSUM encoder of §5e plus the check and decoder) sit at gate 3 / rounds 3, while the best networks (41 CSUMs at d = 4,
43 at d = 8, depth 35–37, N = (3,0,0,0,0,1,1) resp. (7,0,0,0,0,1,1)) reach **gate 1.0 / rounds 2.0** (17 with τ_F = 15),
the same at d = 4 and d = 8 (dense d = 8: 3.71 / 6.14).  The inner coefficient is untouched (C₃ = 85 at d = 4, 185 at
d = 8 for every basis).  Floor at p = 10⁻³, D = 21: 1.4×10⁻¹¹ vs 6.2×10⁻¹¹ dense (4.3× lower).  `footprint escape-gap
census-sparse` (`host_census.EXPOSURE_713_SPARSE`):

| target | census, dense | census, best sparse |
|---|---|---|
| 10⁻¹⁰ at 10⁻³ | D = 21, 7 112 | D = 21, 7 112 (err 1.4×10⁻¹¹) |
| 10⁻¹² at 10⁻³ | D = 27, 11 720 | **D = 25, 10 056** (err 3.1×10⁻¹³) |
| 10⁻¹⁵ at 10⁻³ | D = 33, 17 480 | D = 31, 15 432 |
| 10⁻¹² at 5×10⁻⁴ | D = 19, 5 832 | D = 19, 5 832 |

**Two-wire failures and d-scan** (`host_census.faults_pairs`, `scratchpad/census_dscan.py` →
`runs/meas_distill/census_dscan.jsonl`, 2026-11-02).  Correlated logical Z ⊗ Z failures after every CSUM, either uniform
(all Z^a ⊗ Z^b, weight 1/(d−1)²) or the copy pattern a CSUM makes from one earlier wire error (Z_c^−a Z_t^a, weight
1/(d−1)), counted as undetected bad outputs per unit two-wire failure rate:

| d | dense: gate / rounds | two-wire uniform / copy | best sparse: gate / rounds | two-wire uniform / copy |
|---|---|---|---|---|
| 4 | 4.67 / 5.33 | 1.89 / 1.00 | 1.0 / 2.0 | 0 / 0 |
| 8 | 3.71 / 6.14 | 1.49 / 1.00 | 1.0 / 2.0 | 0 / 0 |
| 16 | 3.33 / 6.07 | 1.25 / 1.00 | 1.0 / 2.0 | 0 / 0 |
| 32 | 3.16 / 6.03 | 1.12 / 1.00 | 1.0 / 2.0 | 0 / 0 |

The best sparse gadget (same network seed at every d: 41 CSUMs at d = 4, 43 at d ≥ 8) is immune to correlated two-wire
CSUM failures — every Z ⊗ Z pattern after a CSUM is caught by the check — so the single-wire caveat above does not apply
to the recipe of record, and its exposures 1.0 / 2.0 are d-independent through d = 32.

10⁻¹⁰ stays at D = 21 (at D = 19 the floor alone is 9.9×10⁻¹¹, no room for the other terms); 10⁻¹² comes within 56 qudits of the budget.

## 5k. Idle noise at every tick (`tick_idle`, `zdsim.add_tick_idle`)

Every model above charges idling once per round (data DEP1 at the reset), and the cultivation circuits not at all.  The
circuit-level "idle at every tick" model adds DEP1(f p) on every live qudit a tick leaves alone: f = 1 is the SD6-style
uniform model (idle as noisy as a gate), f = 0.1 SI1000-like.  Surface memory and logical CSUM: `tick_idle = f` in
`surface_zd.memory_circuit` / `logical_csum.csum_circuit` (F / FD ticks, each of the 4 CSUM steps, the ancilla readout);
any other circuit: `zdsim.add_tick_idle(circ, f p)` (ASAP schedule in circuit order, just-in-time resets, k waiting ticks
composed into one DEP1).

**Host** (d = 4, `zd_d4_tick{1,0.1}.jsonl`, 64 points each, local, 2026-11-02):

| f | p_th (Z / X) | A (Z / X) |
|---|---|---|
| 0 (§5g) | 7.8 / 7.2 ×10⁻³ | 0.046 / 0.053 |
| 0.1 | 7.2 / 7.1 ×10⁻³ | 0.049 / 0.059 |
| 1 | 5.4 / 5.5 ×10⁻³ | 0.084 / 0.092 |

**Logical CSUM** (`csum_d4_tick{1,0.1}.jsonl`, D = 5, 7, 2×10⁵ shots, loranne): CSUM / two-patch-memory ratio 1.03–1.33
at f = 1 and 1.06–1.26 at f = 0.1, the same as without idle (§5i); the excess is 0.5–9 patch-rounds at these statistics
(±2), consistent with ≈ 0.45 D; D = 9, f = 1 (lenore, 4×10⁵ shots, p = 10⁻³ / 1.5×10⁻³): ratio 1.17 (95 / 81 fails) / 1.14
(718 / 631).  So the CSUM and census fits scale with the host, which is what
`scratchpad/footprint_tick_proxy.py` assumes.  Footprint at p = 10⁻³ (escape-gap, sparse census):

| target | f = 0 | f = 0.1 | f = 1 |
|---|---|---|---|
| 10⁻¹⁰ | D = 21, 7 112 | D = 21, 7 112 | D = 25, 10 056 |
| 10⁻¹² | D = 25, 10 056 | D = 25, 10 056 | D = 29, 13 512 |
| 10⁻¹⁰ at 5×10⁻⁴ | D = 15, 3 656 | D = 15, 3 656 | D = 17, 4 680 |

— with the cultivation and escape inputs held at their f = 0 values, which the next paragraph shows is not tenable at f = 1.

**Cultivation** (`cult_sweep --tick f --parallel`, d = 4, grow 3 → 5 r1c2, 10⁷ shots per point, local).  The serial
circuits of §CULTIVATION_PLAN (plaquettes one at a time, linear cat fan-outs) carry ~5 200 idle qudit-ticks at D = 5
against ~860 gate locations, so `ColorCult` / `GrowCult(parallel=True)` adds a parallel schedule: block X / Z plaquette
rounds keeping each plaquette's serial leg order, log-depth cat trees, a second cat register (75 qudits instead of 56),
and a serial final round (the block schedule leaves that round's Z hooks undetected; with it the single-fault scan is 0
again, also under tick idle).  That still leaves ~1 600 idle qudit-ticks (≈ 0.7 idle faults per shot at f = 1, p = 10⁻³):

| f | p = 5×10⁻⁴: accept, p_L | p = 10⁻³: accept, p_L |
|---|---|---|
| 0 (serial, §12) | 0.66, 1.5×10⁻⁷ | 0.45, 9.1×10⁻⁷ |
| 0.1 (parallel) | 0.62, 1.1×10⁻⁶ (7 / 6.2M) | 0.39, 1.05×10⁻⁵ (41 / 3.9M) |
| 1 (parallel) | 0.35, 2.3×10⁻⁴ (789 / 3.5M) | 0.12, 1.7×10⁻³ (2 102 / 1.2M) |

Failures are second order (fault distance 2 holds): long waits compose into DEP1 of a few percent.  At f = 1 these
circuits cannot feed the [[7,1,3]] round, and at f = 0.1 the input error rises 7–12×; the cultivation circuit's depth,
not the host, is the f = 1 bottleneck.  (Superseded by the next two paragraphs: that table was taken without hook
sampling and with a just-in-time-prep bug in `add_tick_idle`; both are fixed, and the f = 1 acceptance at 10⁻³ is 0.29.)

**Hook errors and the end-to-end boundary.**  The cultivation `CSUM` primitive used to apply DEP2 once after a whole
fan-out layer, so an ancilla's hook (a fault part-way through its fan-out) was never sampled.  Only cultivation has such
layers (surface, CSUM, escape and [[7,1,3]] gadgets do not).  With `hooks=True` (DEP2 after each pair) every harmful single
fault sits in the final stabiliser round, where the ideal lookup decode of §12 has nothing after it to catch it: the
lookup p_L at 10⁻³ becomes 1.6×10⁻⁴ (serial) / 4.4×10⁻⁴ (parallel), not 9×10⁻⁷.  That boundary is too pessimistic — in the
pipeline the escape comes next.  `scratchpad/cult_escape_e2e.py` injects each accepted shot's bare residual frame
(`FrameSim(gauge=False)`, qudits mapped by the plaquette isomorphism) into the escape after its noiseless colour-state
prep, in both bases, then runs the escape's own noise, postselection, 2-adic matcher and gap cut 15; a shot fails iff
the net logical (a′, b′) is not a power of the cultivated state's check operator, b′ ≢ a′ (mod 4) (calibrated by injecting
X_L, Z_L, X_L Z_L).  The escape's intrinsic error on trivial frames (`ESCAPE_Z4_GAP`, 3×10⁻⁷) adds on top.
Both bases run with the same escape noise (`Escape.run(chunk_seeded=True)`, same seed: the joint keep decision is one
escape's, agreement 0.9998; an earlier pass with independent seeds squared the residual frames' keep rate and is
superseded), and 2×10⁴ trivial frames per run measure the escape's own keep rate a(cut), so the error is per kept
escape output, p_out = fails / ((accepted − residual) a + kept residual).  Results (`cult_escape_e2e.jsonl`, lines with
`keep_disagree`; 4×10⁶ cultivation shots each at p = 10⁻³, d = 4, grow 3 → 5; 2026-11-03):

| cultivation | escape idle | accept | lookup p_L | p_out, no cut | p_out, cut 12 | p_out, cut 15 | cult. attempts / kept escape (cut 15) |
|---|---|---|---|---|---|---|---|
| serial, no hooks (old model) | 0 | 0.452 | 1.1×10⁻⁶ | 1.6×10⁻⁵ | 0 | 0 (< 3.0×10⁻⁶) | 4.1 |
| serial + hooks | 0 | 0.435 | 1.6×10⁻⁴ | 1.1×10⁻⁴ | 0 | 0 (< 3.2×10⁻⁶) | 4.2 |
| parallel + hooks, f = 0 | 0 | 0.435 | 4.4×10⁻⁴ | 7.2×10⁻⁵ | 1.0×10⁻⁶ (1) | 0 (< 3.1×10⁻⁶) | 4.2 |
| parallel + hooks, f = 0.1 | 0 | 0.417 | 4.6×10⁻⁴ | 1.2×10⁻⁴ | 2.2×10⁻⁶ (2) | 0 (< 3.3×10⁻⁶) | 4.5 |
| parallel + hooks, f = 0.1 | 0.1 | 0.417 | 4.6×10⁻⁴ | 1.1×10⁻⁴ | 2.5×10⁻⁶ (2) | 1.3×10⁻⁶ (1) | 5.4 |
| parallel + hooks, f = 1 | 0 | 0.288 | 9.4×10⁻⁴ | 5.3×10⁻⁴ | 6.8×10⁻⁶ (4) | 1.8×10⁻⁶ (1; < 8.3×10⁻⁶) | 7.0 |
| parallel + hooks, f = 1 | 1 | 0.288 | 9.4×10⁻⁴ | 7.0×10⁻⁴ | 0 (< 1.5×10⁻⁵) | 0 (< 2.6×10⁻⁵) | 34.1 (cut 12: 24.9) |
| one check (r1c1), f = 0 | 0 | 0.520 | 4.5×10⁻⁴ | 8.5×10⁻⁵ | 1.0×10⁻⁵ (12) | 1.1×10⁻⁵ (12) | 3.5 |
| one check, f = 1 | 0 | 0.350 | 9.6×10⁻⁴ | 5.7×10⁻⁴ | 1.4×10⁻⁵ (10) | 1.1×10⁻⁵ (8) | 5.7 |
| one check, f = 1 | 1 | 0.350 | 9.6×10⁻⁴ | 7.1×10⁻⁴ | 1.0×10⁻⁵ (2) | 7.2×10⁻⁶ (1) | 28.9 |

(Bounds 95% CL; the trivial frames never fail at cut ≥ 10.)  The escape removes the hook failures of the two-check
recipe: ≤ 3×10⁻⁶ at f ≤ 0.1 and ≈ 2×10⁻⁶ at f = 1 with an idle-free escape — the same order as the old (hook-free)
9×10⁻⁷.  The one-check recipe leaves ≈ 1.1×10⁻⁵ that no cut touches (wrong states the escape cannot see; 12 events
at every cut 8–18), 10× the two-check recipe for 16% fewer attempts.  The escape's own floor with tick idle
(`graft_variants.jsonl`, variants one+gap-tick{0.1,1}, both bases): f = 0.1, 2×10⁵ shots/basis: discard 49% (cut 15: 52%),
0 fails at cut ≥ 10 (< 3×10⁻⁵); f = 1, 1.5×10⁶ shots/basis: discard 81.5% from the postselection itself (cut 12: 83.5%,
cut 15: 86%), 690 fails before the cut, 1 at cut 12, 0 at cut 15 (< 1.4×10⁻⁵ per basis).  So at f = 1 the escape's
postselection, not the cultivation error, sets the cost.  Footprint (`footprint_tick_proxy.py f cult=a,p escape-gap
census-sparse`): the [[7,1,3]] round is Clifford-floor limited up to p_in ≈ 4×10⁻⁵, so the qudits do not move —
**f = 0: 10⁻¹⁰ at 10⁻³ → D = 21, 7 112 qudits, 63 cultivation attempts/output; f = 0.1 (both): D = 21, 81 attempts;
f = 1, idle-free escape: D = 25, 10 056 qudits, 105 attempts; f = 1 everywhere: D = 25, 10 056 qudits, 511 attempts
(cut 15) or 374 (cut 12); 10⁻¹² → D = 29, 13 512.**
**Lighter escape postselection at f = 1** (variant `none+gap-tick1`: no colour-region postselection, first-order
certification and the gap cut only).  Per sector (10⁶ shots/basis, `graft_variants.jsonl`): Z cut 13 / 14 / 15 →
discard 52 / 61 / 68%, p_fail 2.1×10⁻⁵ / 7.7×10⁻⁶ / 0 (< 9×10⁻⁶); the X gap never reaches 14 (max 13.3), cut 12 / 13 →
62 / 67%, 1.8×10⁻⁵ / 6×10⁻⁶.  Joint two-sector keep on trivial frames with common noise
(`scratchpad/escape_joint_cut.py`, 3×10⁵ frames, `runs/surface/escape_joint_cut.jsonl`):

| variant | cut (Z, X) | joint keep | error (per-sector runs, summed) |
|---|---|---|---|
| one+gap-tick1 (recipe) | (12, 12) | 0.149 | ≈ 4×10⁻⁶ |
| one+gap-tick1 | (15, 15) | 0.108 | < 1.4×10⁻⁵ (0 events) |
| none+gap-tick1 | (13, 13) | 0.223 | ≈ 2.7×10⁻⁵ |
| none+gap-tick1 | (14, 13) | 0.191 | ≈ 1.4×10⁻⁵ |

So dropping the colour postselection buys 1.3× keep at (14, 13) for ≈ 3× the escape error, still under the [[7,1,3]]
tolerance (≈ 4×10⁻⁵): f = 1 everywhere would go from 374 to ≈ 290 cultivation attempts/output, same 10 056 qudits.  Not
adopted — the gain is small, its residual-frame (end-to-end) behaviour is untested, and its Z/X keep decisions disagree
on 6% of shots (the certified-pattern discards see the basis-dependent readout), so its joint keep is uncertain at that
level.  The f = 1 cost is set by the escape's colour-region decoding at 1.7× more locations, not by the cut.
**End to end at 5×10⁻⁴** (lenore, 2×10⁷ cultivation shots per config, common-noise escape, `cult_escape_e2e.jsonl`):

| config | escape | cult. accept | p_out no cut | cut 12 | cut 15 | attempts / kept escape (cut 15) |
|---|---|---|---|---|---|---|
| parallel + hooks, f = 0 | idle-free | 0.660 | 2.9×10⁻⁵ | 6.1×10⁻⁷ (6) | 2.0×10⁻⁷ (2) | 2.04 |
| parallel + hooks, f = 0.1 | f = 0.1 | 0.646 | 5.2×10⁻⁵ | 3.3×10⁻⁷ (3) | 1.1×10⁻⁷ (1) | 2.21 |
| parallel + hooks, f = 1 | idle-free | 0.536 | 2.5×10⁻⁴ | 2.7×10⁻⁶ (21) | 9.1×10⁻⁷ (7) | 2.60 |
| parallel + hooks, f = 1 | f = 1 | 0.536 | 2.7×10⁻⁴ | 3.5×10⁻⁶ (15) | 1.7×10⁻⁶ (7) | 4.87 |

The f = 1 escape keeps 40% of trivial frames at cut 15 (vs 14% at 10⁻³), so the escape cost that dominates at 10⁻³
(34 attempts per kept escape) shrinks to 4.9 at 5×10⁻⁴.  Footprint (`footprint_tick_proxy.py f e2e=p,attempts,p_in
escape-gap census-sparse`, the end-to-end numbers replacing cultivation + escape): **10⁻¹⁰ at 5×10⁻⁴: f ≤ 0.1 → D = 15,
3 656 qudits, 31–33 attempts/output; f = 1 → D = 17, 4 680 qudits, 39 attempts (idle-free escape) or 73 (f = 1
everywhere); 10⁻¹² → D = 19, 5 832 (f ≤ 0.1) / D = 21, 7 112 (f = 1).**

**What p_out counts.**  p_out in the end-to-end tables has only the failures of frames that cultivation left nontrivial in
its numerator; the clean frames, the great majority of kept outputs, enter the denominator only (20 000 of them per run, too
few to see their failure rate).  The total input error is p_in = p_out + the escape's own floor on clean frames, taken from
the separate per-sector runs (`graft_variants.jsonl`): negligible for an idle-free escape (5×10⁻⁸), ≈ 4×10⁻⁶ (cut 12) and
< 2.8×10⁻⁵ summed (cut 15, 0 events) for one+gap-tick1 at 10⁻³; not measured at f = 1 below 10⁻³.  The footprints do not
move: rerunning them with p_in = 3×10⁻⁵ at 5×10⁻⁴ and 7×10⁻⁴ gives the same D, since the round is Clifford-floor limited up
to p_in ≈ 4×10⁻⁵.

**End to end at 7×10⁻⁴** (lenore, 10⁷ shots): parallel + hooks, f = 0: accept 0.558, cut 15 2.7×10⁻⁷ (1), 2.73 attempts per
kept escape; f = 1 everywhere: accept 0.418, cut 15 2.0×10⁻⁶ (2), escape keep 0.255, 10.0 attempts.  Attempts per kept
escape at f = 1 everywhere, cut 15: 34.9 / 10.0 / 4.9 at 10⁻³ / 7×10⁻⁴ / 5×10⁻⁴.  Footprint: 10⁻¹⁰ → f = 0: D = 17, 4 680
qudits, 41 attempts/output; f = 1: D = 19, 5 832, 151 attempts; 10⁻¹² → D = 21, 7 112 / D = 25, 10 056.

**Lighter escape end to end at f = 1, 10⁻³** (lenore, 1.6×10⁷ shots, same cultivation frames, joint cuts (Z, X)):

| escape | cut (Z, X) | joint keep | p_out (events) | floor (§ above, per-sector) | attempts per kept | attempts/output |
|---|---|---|---|---|---|---|
| one+gap-tick1 | (12, 12) | 0.150 | 6.4×10⁻⁶ (4) | ≈ 4×10⁻⁶ | 25.4 | 381 |
| one+gap-tick1 | (15, 15) | 0.111 | 6.5×10⁻⁶ (3) | < 2.8×10⁻⁵ | 34.9 | 523 |
| none+gap-tick1 | (14, 13) | 0.188 | 1.3×10⁻⁶ (1) | ≈ 1.4×10⁻⁵ | 20.4 | 305 |
| none+gap-tick1 | (13, 13) | 0.221 | 1.1×10⁻⁶ (1) | ≈ 2.7×10⁻⁵ | 17.3 | 260 |

The residual frames do not punish the lighter escape: its p_out is no worse than the recipe's, so the untested worry above
is removed and the cost would fall by 20% at (14, 13) or 32% at (13, 13), with the same D = 25, 10 056 qudits.  Still not
adopted: its error is dominated by a floor resting on few events, (13, 13) sits at ≈ 2.8×10⁻⁵ against the ≈ 4×10⁻⁵
tolerance, and its Z and X keep decisions disagree on 10% of shots here (one+gap-tick1: 4×10⁻⁵).  The joint keep scores only
shots both basis runs keep, so failures on the other 10% are unsampled and p_out for this variant is uncertain until a
single-readout keep rule is defined and measured.

## 6. Open / next
1. Done (§5d): concatenated order-3 stages; partial discard ruled out.
2. Done: w = 50 outer on 9500 inputs for [[112,50]] (γ 1.372, ≈ 9.6k qudits).  Open: smaller w = 26 outer (< 2080).
3. Done (§5e): circuit-level model; first-order Clifford floor a_C = 0.36 exact for [[7,1,3]] (cannot be made second order without consuming encoded outputs); sparse encoders 1.5–2.6× fewer CSUMs.  Open: a structured (recursive) sparse normal basis for RM(1,m) to remove the heavy tail.
4. Done: the 10⁸-shot growth run (CULTIVATION_PLAN §12, `sweep_grow35_d4_long.json`) has 41 events at 10⁻³ and 10 at 5×10⁻⁴ (p_L 9.1×10⁻⁷ ± 16 %, 1.5×10⁻⁷ ± 32 %); §5d uses these numbers.
5. Done (§5g): ring Z₄/Z₈ surface-code threshold ≈ 0.8% with a 2-adic matching decoder, A ≈ 0.053, p_th ≈ 7.6×10⁻³ → 5.8k physical qudits at p = 10⁻³ for 10⁻¹⁰.  Open: idle noise during CSUM steps; a joint Z_d decoder.
6. Done (§5h): the escape (Gidney's d₁ = 5 → D = 15 graft, signs consistent for all 2^k) simulated with decoder-certified first-order postselection: 1.8×10⁻⁴ per kept attempt at p = 10⁻³ with one colour class postselected (43% joint discard), 6.6×10⁻⁵ with all three (62%); fault distance 2 from colour-region decoding, not from hooks (flags tested, no gain).  Second-order (pair) certification with a cost cap: 9×10⁻⁵ at 43%.  A BP+LSD hypergraph decoder (hyperdecoder.py) is better at first order but 2–3× worse in bulk; BP-based gap is too slow / uninformative.  **Complementary-gap postselection with the 2-adic matcher** (split boundary + parity potential, exact two-class decode, 0.7 ms/shot): cut at 13 → no failure in 7.9M kept attempts per sector at p = 10⁻³ (< 6×10⁻⁷ total, 90% CL; floor run 2026-10-31) for 0.6 points of discard — the escape is no longer the bottleneck; one [[7,1,3]] round reaches 10⁻¹² at p = 10⁻³ (8 520 qudits, 59 attempts/output).  Gap recipe at 5×10⁻⁴: < 6×10⁻⁶ at 24.9% discard.  Open: idle noise during CSUM layers; margin of the cut (failing gaps reach 12.x).
7. Done (§5i): the logical CSUM measured on two surface-code patches with a Pauli-frame (sequential) decoder: 1.1–1.2× the idle two-patch memory, same threshold; independent per-patch decoding has a first-order floor (24 p at D = 9) from source measurement errors leaking into the cross-layer detectors.  Footprint: 10⁻¹⁰ unchanged to +1 step (5 832–7 112 qudits), 10⁻¹² one step up (D = 25, 10 056).  Done: the d = 8 sweep (A 0.25/0.22, p_th 7.3/6.5×10⁻³) matches.  Open: idle census between CSUM layers.
8. Done (§5h, d = 8 paragraph): the escape at d = 8 reproduces d = 4 (one class 1.6×10⁻⁴ at 45% discard; gap cut 13: 0 failures in 1.09×10⁵ kept per sector at p = 10⁻³).  d = 8 cultivation growth 3 → 5 at 2×10⁷ shots (CULTIVATION_PLAN §12): acceptance 0.439 / 0.662 at p = 10⁻³ / 5×10⁻⁴, 3 failures in 8.8M accepted (3.4×10⁻⁷, 90% CL 0.9–8.8×10⁻⁷) and 0 in 13.2M (< 1.7×10⁻⁷), consistent with d = 4 (9.1×10⁻⁷, 1.5×10⁻⁷).  Gap recipe measured at 5×10⁻⁴ (§5h).  Running: 10⁸-attempt floor of the gap recipe at d = 4.
9. Done (§5j): host-code census of the Clifford layer — the depolarizing a_C = 0.357 under-counts the surface-code host by ~10× (single-wire Z-type failures are the surviving type; idle rounds never charged); floor = 4.67 q_gate + 5.33 q_round per output.  Footprint: 10⁻¹⁰ at 10⁻³ → D = 21, 7 112 qudits; 10⁻¹² → D = 27 (11 720).  The CSUM excess measured at fixed D (10⁶ shots): 2.0 ± 0.3 patch-rounds at D = 5, 3.4 ± 0.4 at D = 7, ≈ 0.45 D — the fit's value, not one extra round.  Done (§5j, sparse census): the best of 200 sparse encoder networks has exposures 1.0 / 2.0 (4.3× lower floor), same at d = 8; 10⁻¹² at 10⁻³ → D = 25, 10 056.  Done (§5g): d = 16 host sweep — p_th unchanged, A +15–25%, no footprint change.  Done (§5g): d = 32 host sweep — A, p_th equal to d = 16 (saturated).  Done (§5j): two-wire CSUM failures — the best sparse gadget is immune (0 / 0); census exposures d-independent through d = 32.  Done (§5h): gap-cut margin — cut 13 gives 3.9×10⁻⁷ (8 events in 20.5M kept per sector), cut 15 gives 5×10⁻⁸ for +0.5 points of discard; the recipe now cuts at 15, no footprint change.  Done (§5j): p-dependence of the excess — ratio 1.10–1.15 at every p from 1.5×10⁻³ to 4×10⁻³, ≈ 0.55 D, no footprint change.  Done (§5h): d = 16 and d = 32 escapes reproduce d = 8.  Done (§5k): idle noise at every tick — host costs one distance step at f = 1 (10⁻¹⁰ at 10⁻³ → 10 056), nothing at f = 0.1; the CSUM ratio is unchanged; cultivation scored with hooks against the escape (end to end, common escape noise): ≤ 3×10⁻⁶ at f ≤ 0.1, ≈ 2×10⁻⁶ at f = 1; with tick idle in the escape too, f = 1 still costs D = 25, 10 056 qudits for 10⁻¹⁰ at 10⁻³ but 374–511 cultivation attempts/output (escape discard 81.5–86%); lighter escape postselection buys only 1.3× (not adopted).
