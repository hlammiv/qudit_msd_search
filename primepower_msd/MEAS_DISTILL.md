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
| 50 / 71 | — | 9500 / 7100, 9940 (round 0 only) |

Minimum size grows ≈ 3–4 w², so w ≲ 50 under 10k inputs.

**Inner codes with D ≥ 7** (`inner_search2.py`: anneal on signed flat collections, objective = number of
≤ 3 + ≤ 3 syndrome collisions, which is 0 iff d_bin(S^⊥) ≥ 7; aligned forbidden flats cap n):
`runs/meas_distill/inner7b/best_inner7.json` — all exactly balanced, normal, Z₄ d > 6, Z₈ d > 5, D ≥ 7 by the 2-adic theorem.

| inner | T/out = 1 + 6n/k | γ (order 7) | outer w | best certified outer |
|---|---|---|---|---|
| [[63,21,7]] signed RM(2,6) | 19.00 | 1.513 | 21 | 1470 inputs |
| **[[64,26,≥7]]** | **15.77** | **1.417** | 26 | **2080 inputs** |
| [[96,40,≥7]] | 15.40 | 1.405 | 40 | 5200 inputs |
| [[112,50,≥7]] | 14.44 | 1.372 | 50 | open |

[[63,21,7]] hosting re-verified at d = 4, 8 (with distance), 16, 32, 64 (balance, normal basis, all logical S
exponents odd): `runs/meas_distill/rm63_alld.log`.

**Recommended protocol: [[64,26,≥7]] inner + (3,26) outer on 2080 inputs (240 checks).**  Order 7, 15.77 T per
output, γ = 1.417, ≈ 2080 + 64 + ancillas ≈ 2.2k qudits, valid for every d = 2^k.

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

  Dominant term = N·(d−1)·n³ (one bad input hidden in each of its 3 checks by a layer-1/layer-2 pair Z^a, Z^{−a} on
  one inner qudit), so p_out ≈ n³ p⁷/(d−1)⁶ (d = 8: ≈ 2 p⁷).
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
| 1e-10 | 69 | 9.8 | 12.5 |
| 1e-15 | 69 | 47 | **12.5** ([[127,113,3]]) |
| 1e-20 | 69 | 47 | 36 ([[64,26,7]]) |

At physical p = 5×10⁻⁴: 1e-15 costs 7.3 (D3) / 8.9 (growth) attempts per output.

Clifford-noise proxy (each location also fails at g·p_C, g = 10): no change for p_C ≤ 10⁻⁷; costs rise from
p_C ≈ 3×10⁻⁷ (for example 3→5 growth at 1e-15: 13.3 → 50 at p_C = 10⁻⁶).  So the logical Cliffords hosting distillation
need p_C ≲ 10⁻⁷.

## 6. Open / next
1. Done (§5d): concatenated order-3 stages; partial discard ruled out.
2. w = 50 outer (≤ 10k) for [[112,50]] (γ 1.372); smaller w = 26 outer (< 2080).
3. Clifford noise: currently a proxy (§5d); an explicit encoder/decoder circuit model would sharpen the p_C ≲ 10⁻⁷ requirement.
4. Cultivation acceptance and p_L at d = 4 rest on few events (4 failures at 10⁻³); longer lenore sweeps would tighten §5d.
