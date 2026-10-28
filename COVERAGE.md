# Best certified code per prime, p <= 97

*2026-10-28. One row per prime: the lowest-gamma code in the catalog (`qmsd/data/catalog/`) with a
certified or closed-form distance, and the structural reason the row cannot be improved where we
have one. gamma = ln(n/k)/ln d.*

| p | best code | gamma | distance evidence | status / why |
|---|---|---|---|---|
| 3 | [[206,37,4]]_3 (paper's puncture set, A_d 572 ours) | 1.239 | MITM exact | gamma < 1 impossible for unweighted RM (flat cap); weighted m <= 7 probed, no window (`WEIGHTED_SUPPORTS.md`) |
| 5 | [[519,106,5]]_5 (paper) | 0.987 | MITM exact (10 s) | the only p <= 13 crossing; Singleton-infeasible at m <= 3 |
| 7 | [[1756,260,7]]_7 (weighted sphere, ours) | 0.982 | three-part proof + weight-7 witness | first gamma < 1 at p = 7; 130 further d >= 7 candidates at k = 265-270 await certification |
| 11 | P11_ROW | | | 2D-capped for unweighted RM; weighted supports dead by line words (r < p-1) |
| 13 | P13_ROW | | | same as p = 11 |
| 17 | [[237,52,6]]_17 (ours) | 0.847 | MITM exact, two-way | first small gamma < 1 search code; 330x smaller than the analytic family |
| 19 | [[293,68,5]]_19 (ours) | 0.908 | MITM exact | |
| 23 | [[18,5,4]]_23 (RS) | 0.924 | MITM exact; eps* = 0.396 | Reed-Solomon m = 1 family from here on: closed-form d = r - k + 2, A_d = C(n,d)(p-1) |
| 29 | [[22,7,4]]_29 | 0.826 | MITM exact; eps* = 0.353 | |
| 31 | [[25,6,5]]_31 | 0.887 | MITM exact; eps* = 0.416 | |
| 37 | [[29,8,5]]_37 | 0.800 | MITM exact; eps* = 0.383 | |
| 41 | [[31,10,5]]_41 | 0.703 | MITM exact; eps* = 0.371 | |
| 43 | [[34,9,6]]_43 | 0.742 | closed-form MDS | MITM tables scale as (p-1)^3: infeasible from here |
| 47 | [[36,11,6]]_47 | 0.662 | closed-form MDS | |
| 53 | [[40,13,6]]_53 | 0.627 | closed-form MDS | |
| 59 | [[45,14,7]]_59 | 0.600 | closed-form MDS | |
| 61 | [[47,14,7]]_61 | 0.622 | closed-form MDS | |
| 67 | [[51,16,7]]_67 | 0.596 | closed-form MDS | |
| 71 | [[53,18,7]]_71 | 0.555 | closed-form MDS | |
| 73 | [[56,17,8]]_73 | 0.573 | closed-form MDS | |
| 79 | [[60,19,8]]_79 | 0.553 | closed-form MDS | |
| 83 | [[62,21,8]]_83 | 0.521 | closed-form MDS | |
| 89 | [[67,22,9]]_89 | 0.507 | closed-form MDS | |
| 97 | [[73,24,9]]_97 | 0.506 | closed-form MDS | |

Notes.

- The RS rows are the minimum-gamma member of the family at each prime (`scratchpad/rs_family.py`,
  bundle `runs/coverage/rs_family_p23_97.json`); every k is admissible and the distance does not
  depend on which columns are punctured. Exact single-round thresholds (arXiv:2610.08012 model) are
  stored for p <= 41 (`qmsd/data/performance/`).
- p = 11 and 13 rows come from the m = 2 flat-spread search on loranne (`runs/coverage/p1[13]_m2_k*.json`),
  the first cataloged codes at those primes. The no-go (`D_NEGATIVE_RESULT.md`) predicts the 2D cap
  [[100,21,4]]_11 (gamma 1.126) and [[143,26,5]]_13 (gamma 1.059) as the best RM constructions.
- The analytic Manhattan family (astronomical n) is cataloged separately and is not a "small code".
