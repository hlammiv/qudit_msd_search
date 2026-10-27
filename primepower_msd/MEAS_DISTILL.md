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

**Footprint** (`footprint.ESCAPE_Z4` = *one*, `ESCAPE_Z4_ALL` = *all*, `ESCAPE_Z4_PAIRS` = *one+pairs-cap*;
`python -m primepower_msd.footprint escape` / `escape-all` / `escape-pairs`): the escape error adds to p_cult (9×10⁻⁷)
and dominates it — the cultivated state leaves the escape at 8×10⁻⁵–1.8×10⁻⁴, worse than the 2×10⁻⁵ of a d = 3
cultivation, so the escape, not the Clifford floor, is now the pipeline's bottleneck:

| p | target | pipeline | physical qudits | cultivation attempts / output: one+pairs, one, all (§5f free escape) |
|---|---|---|---|---|
| 10⁻³ | 10⁻¹⁰ | one [[7,1,3]] round, D = 19 | 5 832 | 58, 58, 87 (33) |
| 10⁻³ | 10⁻¹² | two [[7,1,3]] rounds, D = 23 | 8 520 | 877, 876, 1 311 (33) |
| 5×10⁻⁴ | 10⁻¹⁰ | one [[7,1,3]] round, D = 15 | 3 656 | 30, 30, 37 (22) |
| 5×10⁻⁴ | 10⁻¹² | one [[7,1,3]] round, D = 17 | 4 680 | 30, 30, 37 (22) |

The qudit count is unchanged (the 8 logical qudits are reused round after round); the escape costs a factor 1.4–1.9 in
time through its discard, and at p = 10⁻³ the 10⁻¹² target needs a second round with any of the recipes
(3.15 × (8.8×10⁻⁵)³ = 2.1×10⁻¹²; an escape below 7×10⁻⁵ would remove it).  Caveats: no idle noise during CSUM layers; the final-code and D = 15
baselines were run without postselection; the 5×10⁻⁴ points rest on 1–3 events each.

## 6. Open / next
1. Done (§5d): concatenated order-3 stages; partial discard ruled out.
2. Done: w = 50 outer on 9500 inputs for [[112,50]] (γ 1.372, ≈ 9.6k qudits).  Open: smaller w = 26 outer (< 2080).
3. Done (§5e): circuit-level model; first-order Clifford floor a_C = 0.36 exact for [[7,1,3]] (cannot be made second order without consuming encoded outputs); sparse encoders 1.5–2.6× fewer CSUMs.  Open: a structured (recursive) sparse normal basis for RM(1,m) to remove the heavy tail.
4. Cultivation acceptance and p_L at d = 4 rest on few events (4 failures at 10⁻³); longer lenore sweeps would tighten §5d.
5. Done (§5g): ring Z₄/Z₈ surface-code threshold ≈ 0.8% with a 2-adic matching decoder, A ≈ 0.053, p_th ≈ 7.6×10⁻³ → 5.8k physical qudits at p = 10⁻³ for 10⁻¹⁰.  Open: idle noise during CSUM steps; a joint Z_d decoder.
6. Done (§5h): the escape (Gidney's d₁ = 5 → D = 15 graft, signs consistent for all 2^k) simulated with decoder-certified first-order postselection: 1.8×10⁻⁴ per kept attempt at p = 10⁻³ with one colour class postselected (43% joint discard), 6.6×10⁻⁵ with all three (62%); fault distance 2 from colour-region decoding, not from hooks (flags tested, no gain).  Second-order (pair) certification on the colour region with a first-order cost cap: 8.8×10⁻⁵ at 43% — the recipe of record (genuine pair mass is only 1/5 of the p² coefficient; the residual is third order).  A BP+LSD hypergraph decoder (hyperdecoder.py) is better at first order but 2–3× worse in bulk, tuning included; BP-based gap postselection is too slow / uninformative.  Footprint unchanged in qudits (5 832 at p = 10⁻³ for 10⁻¹⁰), attempts 58 per output instead of 33.  Open: a matching-based decoder gap (logical representative on a boundary); a larger colour code (3 → 5 → 7 cultivation); more statistics at 5×10⁻⁴.
