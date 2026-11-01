# Route: IP / CP / SAT for the minimum-distance certificate

Question: can an off-the-shelf solver prove "ker(G0) has no nonzero vector of weight <= w" (w = 5, 6)
for the project's parity-check matrices (15-60 rows over F_p, 400-12 000 columns) faster than the
numba meet-in-the-middle scans in `qmsd/mindist_nb.py`?

**Verdict: no.  Not competitive anywhere in the regime of interest, by many orders of magnitude.**
No solver/encoding could even close weight 3 on the 402-column instance A (MITM: 0.45 s wall,
most of it process + JIT start-up) within 5 minutes, nor weight 3 on a 100-column subset within
200 s.  The solvers finish only for n <= 50-75 at w = 3-4, where MITM is instantaneous.  The
measured growth is at least C(n, w-1) (p-1)^(w-2) (depth-(w-1) enumeration) and in practice worse
(effective exponent ~8 in n between n = 30 and 50 at w = 4), i.e. the solvers never acquire
the birthday advantage of MITM (C(n,3)(p-1)^3 probes against a C(n,3)(p-1)^2 table for w = 6).
Details, table and extrapolations below.  Everything was run with one thread, `ulimit -v 4 GB`.

Files: `ilp_distance.py` (instances, encodings, solver drivers, MITM ground truth; prints a JSON
line), `bench.sh` (sequential sweep driver), `summarize.py` (table from `results.jsonl`),
`results.jsonl` (every run), `cache/*.npy` (the three matrices), `*.log`.

## Instances

| | p | matrix (zero rows dropped) | rank | construction |
|---|---|---|---|---|
| A (negative control) | 7 | 32 x 402 | 21 | class block Gc of plane_sweep k = 270, seed 19 (as `scratchpad/certify_candidate.py`) |
| B (positive control) | 7 | 29 x 400 | 19 | class block of orbit_shift k = 272, seed 3; weight-6 word on class columns [16, 67, 122, 200, 210, 328] |
| C (target scale) | 11 | 30 x 12351 | 30 | `rm_generator(13, 4, 11)` shortened at the 2290 puncture columns of `runs/rm_m4/p11_r13_k2290_seed1.json` |

Subsets are random column subsets (`numpy.default_rng(1)`), zero rows dropped again.  Note that a
C subset with n <= 30 columns is degenerate (the 30 x n block has full column rank, so the kernel
is trivial and every solver proves it in presolve); the smallest meaningful C subset is n = 50.

## Encodings tried

All encodings first row-reduce G0 over F_p (`rref_mod_p`) so only rank-many rows enter the model
(21 for A, 19 for B, 30 for C); the kernel is unchanged.

1. **cpsat-onehot** (OR-tools CP-SAT 9.15).  Booleans x[j][v] = [e_j = v], v = 1..p-1, with
   s_j = sum_v x[j][v] <= 1.  Row t: sum_{j,v} ((R[t,j] v) mod p) x[j][v] = p q_t, q_t an integer in
   [0, floor(w(p-1)/p)] (tight slack bound: residues are in [0, p-1]).  Cardinality sum s_j = w
   (`--exact`, the per-weight semantics of the MITM scans) or 1 <= sum s_j <= w.
   Variant `--lin 0` (no LP relaxation): 5x more branches per second, same outcome.
2. **cpsat-sys** (information-set form).  R = [I | A] after pivoting; e_I = -A e_J.  Only the free
   part e_J is one-hot; per pivot row a residue y_t in [0, p-1] with sum = p q_t + y_t and
   z_t = [y_t != 0]; weight = wt(e_J) + sum z_t.  Fewer free columns (381 for A) but the residue
   freedom destroys what little propagation there was: it could not even close w = 2 on A in 300 s
   (one-hot: 1.75 s).  Dropped after that.
3. **cpsat-elem** (the natural CP model).  w ordered positions: column c_i in [0, n), value v_i in
   [1, p-1], flat index idx_i = (p-1) c_i + v_i - 1, residue r_{t,i} = table_t[idx_i] via
   `AddElement` for each row t, sum_i r_{t,i} = p q_t.  Ordering c_1 < c_2 < ... < c_w removes the
   position permutations; v_1 = 1 removes scalar multiples.  In principle this prunes at depth w-1
   (every row forces the last residue, the element domains intersect to empty), so the tree has
   ~C(n, w-1) (p-1)^(w-2) leaves.  Variant **cpsat-elem-dfs**: a fixed decision strategy
   (c_1, v_1, c_2, ...; `FIXED_SEARCH`) to force exactly that DFS -- slower (560 branches/s on A:
   the decomposed element constraints with 2412-entry tables are expensive to propagate).
4. **sat** (CNF, CaDiCaL 1.9.5 via pysat).  One-hot values with pairwise at-most-one; per row a
   residue automaton over Z_p (state literals st[t][j][r], p(p-1) + p transition clauses per
   column, at-most-one state per step), final state 0; sequential-counter cardinality on s_j.
   A: 64k variables, 540k clauses.  This is the best encoding at tiny n (it shares CaDiCaL's
   inprocessing) but scales the same way.

Symmetry breaking used everywhere: **first nonzero coordinate = 1** (prefix booleans
b_j = OR_{i<=j} s_i, x[j][v >= 2] -> b_{j-1}; in the element model simply v_1 = 1), i.e. the
(p-1)-fold scalar symmetry is removed.  **Support ordering** is implicit in the indicator models
and explicit (c_1 < ... < c_w) in the element model.  **Column-group symmetry**: `--orbit-reps`
adds "support meets a set of orbit representatives" for a known automorphism group acting on
columns; this is only an |orbit|-fold reduction (not exponential), and none of the three test
matrices has a known nontrivial group (random puncture sets), so it was not used in the numbers
below.  Model correctness was validated by restricting B to its known 6-column support: all three
encodings return the (verified) word [1, 6, 2, 2, 2, 6] in < 0.02 s (`--cols 16,67,122,200,210,328`).

## Results on the controls

* **A (must prove infeasible)**: w = 2 proven by cpsat-onehot (1.75 s), cpsat-elem (0.23 s,
  presolve) and CaDiCaL (22 s); cpsat-sys failed (300 s).  **w = 3: every encoding times out at
  300 s** (onehot 491k branches / 48k conflicts; elem 5.9M branches / 384k conflicts; elem-dfs
  169k branches; CaDiCaL was still running after 952 s when killed -- pysat's interrupt() does
  not stop CaDiCaL, fixed afterwards with chunked conflict budgets).  w >= 4 not attempted on the full matrix.  MITM: w = 3 in
  0.45 s wall (JIT load dominated), w = 4 in 1.5 s, w = 5 in 39 s, all single-thread; w = 6 is the
  known 6 min on 8 threads (its 3+3 left table, ~3 GB at n = 300 and ~8 GB at n = 402, does not fit
  under the 4 GB cap used here -- the one place where MITM has a cost the solvers do not).
* **B (must find the weight-6 word)**: cpsat-onehot found nothing in 600 s (1.31M branches),
  CaDiCaL nothing in 422 s (986k conflicts; an earlier run killed at 377 s likewise).  The
  positive control also fails -- single-threaded, no LNS portfolio, but MITM's 3+3 scan finds
  it with early exit in well under its 6-minute full-scan time.

## Scaling table (seconds, one thread; `>T` = time limit T hit; blank = not run)

MITM times below ~0.5 s are process start-up + numba cache load, not scan time.

| inst | n | w | mitm | cpsat-onehot | cpsat-elem | cpsat-elem-dfs | sat-cadical195 | cpsat-sys | cpsat-onehot lin0 |
|---|---|---|---|---|---|---|---|---|---|
| A | 402 | 2 |  | 1.75 | 0.23 |  | 21.99 | >300 |  |
| A | 30 | 3 | 0.47 | 0.34 | 0.10 |  | 0.04 |  |  |
| A | 50 | 3 | 0.39 | >100 | 1.86 |  | 1.87 |  |  |
| A | 75 | 3 | 0.33 | >100 | 5.68 |  | 9.55 |  |  |
| A | 100 | 3 | 0.34 | >200 |  |  |  |  | >200 |
| A | 200 | 3 | 0.39 |  |  |  |  |  |  |
| A | 300 | 3 | 0.45 |  |  |  |  |  |  |
| A | 402 | 3 | 0.45 | >300 | >300 | >300 | >952 | >300 |  |
| A | 30 | 4 | 0.35 | 3.40 | 0.71 |  | 0.27 |  |  |
| A | 50 | 4 | 0.47 | >100 | 41.49 |  | 36.58 |  |  |
| A | 75 | 4 | 0.36 |  | >100 |  | >331 |  |  |
| A | 100 | 4 | 0.31 | >200 |  |  |  |  |  |
| A | 200 | 4 | 0.52 |  |  |  |  |  |  |
| A | 300 | 4 | 0.97 |  |  |  |  |  |  |
| A | 402 | 4 | 1.54 |  |  |  |  |  |  |
| A | 30 | 5 | 0.54 | 70.70 | 4.35 |  | 0.98 |  |  |
| A | 50 | 5 | 0.38 |  |  |  |  |  |  |
| A | 100 | 5 | 2.09 |  |  |  |  |  |  |
| A | 200 | 5 | 3.40 |  |  |  |  |  |  |
| A | 300 | 5 | 12.67 |  |  |  |  |  |  |
| A | 402 | 5 | 39.12 |  |  |  |  |  |  |
| A | 50 | 6 | 0.73 |  |  |  |  |  |  |
| A | 100 | 6 | 7.34 |  |  |  |  |  |  |
| A | 200 | 6 | 30.48 |  |  |  |  |  |  |
| B | 400 | 6 |  | >600 |  |  | >377 |  |  |
| C | 12351 | 2 | 0.04 |  |  |  |  |  |  |
| C | 30 | 3 | 0.50 |  | 0.05 |  | 0.03 |  |  |
| C | 50 | 3 | 0.63 |  | 0.87 |  | 0.89 |  |  |
| C | 75 | 3 | 0.58 |  | 13.86 |  | 19.99 |  |  |
| C | 100 | 3 | 0.36 | >200 |  |  |  |  |  |
| C | 200 | 3 | 0.53 |  |  |  |  |  |  |
| C | 300 | 3 | 0.67 |  |  |  |  |  |  |
| C | 600 | 3 | 1.31 |  |  |  |  |  |  |
| C | 1200 | 3 | 2.85 |  |  |  |  |  |  |
| C | 2400 | 3 | 10.78 |  |  |  |  |  |  |
| C | 12351 | 3 | 348.29 |  |  |  |  |  |  |
| C | 30 | 4 | 0.63 |  | 0.06 |  | 0.03 |  |  |
| C | 50 | 4 | 0.71 |  | 58.13 |  | 19.16 |  |  |
| C | 75 | 4 | 0.66 |  | >100 |  | >136 |  |  |
| C | 100 | 4 | 0.54 | >200 |  |  |  |  |  |
| C | 200 | 4 | 1.71 |  |  |  |  |  |  |
| C | 300 | 4 | 2.16 |  |  |  |  |  |  |
| C | 600 | 4 | 9.42 |  |  |  |  |  |  |
| C | 1200 | 4 | 41.53 |  |  |  |  |  |  |
| C | 30 | 5 | 0.67 |  | 0.07 |  | 0.03 |  |  |
| C | 100 | 5 | 3.06 |  |  |  |  |  |  |
| C | 300 | 5 | 53.13 |  |  |  |  |  |  |
| C | 600 | 5 | 490.48 |  |  |  |  |  |  |

Fits (largest solved points), p = 7 (A): element model w = 3 grows like n^2.7 between n = 50
and 75 (1.86 -> 5.68 s) and w = 4 like n^7.9 between n = 30 and 50 (0.71 -> 41.5 s); CaDiCaL
w = 4 grows n^9.6 (0.27 -> 36.6 s).  p = 11 (C, 30 rows): w = 3 from n = 50 to 75 grows n^6.8
(element, 0.87 -> 13.9 s) and n^7.7 (CaDiCaL, 0.89 -> 20 s); w = 4 at n = 50 takes 58 s / 19 s
and n = 75 is not closed in 100 s by either.  The one-hot CP-SAT model is the worst of all
(cannot close w = 3 at n = 50).  The effective exponents exceed the depth-(w-1) enumeration count
(n^(w-1)) because learning does not generalise across supports: each refuted leaf teaches the
solver essentially nothing about the next one (the row sums are "random" mod p).

## Extrapolation (deliberately optimistic: pure n^w scaling from the largest solved point)

| target | solver estimate | MITM (measured / known) |
|---|---|---|
| A, n = 402, w = 5 | 4.35 s x (402/30)^5 ~ 1.9e6 s ~ 3 weeks (element model) | 39 s (1 thread) |
| A, n = 402, w = 6 | x another ~ n(p-1) ~ 2400: ~ 10^2 years | ~6 min on 8 threads |
| C, n = 12 351, w = 4 | 19 s x (12351/50)^4 ~ 7e10 s ~ 2000 years (CaDiCaL) | n = 1200: 42 s; full n extrapolates to ~1 h of scan (memory-bound 2+2 table under 4 GB) |
| C, n = 12 351, w = 5 | ~ 10^8 years | n = 600: 490 s; full n ~ (12351/600)^3 x 490 s ~ 50 days single-thread (~1 week on 8 threads) |
| C, n = 12 351, w = 6 | ~ 10^13 years | infeasible (~1e18 probes), as stated in the task |

With the measured (super-n^w) exponents the solver numbers are worse still.  Even a 1000x better
solver or encoding would leave the route 3-10 orders of magnitude behind MITM at w = 5, 6.

## Why, and where a solver could still help

The constraints are modular equalities with a free slack per row.  Bounds and LP propagation see
nothing until all but one support column is fixed, and the clause-learning conflict explanations
are as wide as the support itself, so the search degenerates to enumerating (support, values) --
exactly the brute force whose square root MITM takes.  No formulation (indicator, information-set,
ordered-support, CNF automaton) changes this, and the solver's generic presolve / LP machinery adds
a large constant factor on top (CP-SAT ~500-5000 branches/s here).

The route could still be useful for *side problems with strong propagation* that MITM does not
cover: words with heavily restricted structure (e.g. the mixed-word stage, where a few free
coordinates must extend a line word), small-n exact distance with extra constraints, or finding
(not excluding) words in instances expected to have many of them with a portfolio of workers.
For the pure weight-w exclusion on 400-12 000 columns, stay with MITM.

## Ground-truth notes

* MITM under `ulimit -v 4 GB` failed (MemoryError in `_build_left`) for A n = 300 w = 6 (3+3 table)
  and C n = 2400 w = 4 ((p-1)^2 C(n,2) table at p = 11); those cells are absent.
* The full C matrix: w = 2 in 0.04 s, w = 3 in 348 s single-thread (no word, consistent with the
  d >= 4 scan running elsewhere).
