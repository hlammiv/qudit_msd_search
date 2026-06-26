"""Ground-truth numbers transcribed from arXiv:2510.10852 (via the project's verified
notes). These are the deterministic targets
the discovery pipeline must reproduce. Do NOT edit casually -- an error here
corrupts the whole test oracle.

All [[n,k,d]]_p tuples satisfy the build relations:
  - single/maximal-degree RM space: r = r_max = floor((m(p-1)-1)/3)   (3r < m(p-1))
  - dual degree:                     rtilde = m(p-1) - r - 1
  - quantum distance:                d = distance of punctured dual at rtilde
"""

# ---------------------------------------------------------------------------
# Section 3.3 -- single-puncture codes [[p^m - 1, 1, d]]_p (Lemma 2: distance
# is location-independent, d = d_RM(rtilde(r_max,m)) - 1). Strong unit tests for
# the (r_max, rtilde, d_RM) chain end-to-end.  Fields: p, m, n, k, d.
# ---------------------------------------------------------------------------
SINGLE_PUNCTURE_CODES = [
    # p = 3  (n = 3^m - 1, m = 2..8); distances 2,2,5,8,8,17,26
    {"p": 3, "m": 2, "n": 8, "k": 1, "d": 2},
    {"p": 3, "m": 3, "n": 26, "k": 1, "d": 2},
    {"p": 3, "m": 4, "n": 80, "k": 1, "d": 5},
    {"p": 3, "m": 5, "n": 242, "k": 1, "d": 8},
    {"p": 3, "m": 6, "n": 728, "k": 1, "d": 8},
    {"p": 3, "m": 7, "n": 2186, "k": 1, "d": 17},
    {"p": 3, "m": 8, "n": 6560, "k": 1, "d": 26},
    # p = 5  (n = 5^m - 1, m = 1..4); distances 2,3,4,14  (m=1 is Reed-Solomon)
    {"p": 5, "m": 1, "n": 4, "k": 1, "d": 2},
    {"p": 5, "m": 2, "n": 24, "k": 1, "d": 3},
    {"p": 5, "m": 3, "n": 124, "k": 1, "d": 4},
    {"p": 5, "m": 4, "n": 624, "k": 1, "d": 14},
    # p = 7  (m = 1..3); distances 2,4,6
    {"p": 7, "m": 1, "n": 6, "k": 1, "d": 2},
    {"p": 7, "m": 2, "n": 48, "k": 1, "d": 4},
    {"p": 7, "m": 3, "n": 342, "k": 1, "d": 6},
    # p = 11 (m = 1..2); distances 4,7  (m=1 is Reed-Solomon)
    {"p": 11, "m": 1, "n": 10, "k": 1, "d": 4},
    {"p": 11, "m": 2, "n": 120, "k": 1, "d": 7},
]

# ---------------------------------------------------------------------------
# Table 3 -- best triorthogonal codes found by randomized puncture search.
# Fields: p, m, n, k, d, A_d, gamma (paper's value, rounded to 2 dp), and
# has_puncture_json = True iff k > 2 (so explicit columns exist in the oracle).
# gamma = log(n/k)/log(d); the 1-/2-puncture codes are reconstructible from
# Lemma 2 alone (no column list needed).
# ---------------------------------------------------------------------------
TABLE3_CODES = [
    {"p": 3, "m": 4, "n": 80, "k": 1, "d": 5, "A_d": 2080, "gamma": 2.72, "has_puncture_json": False},
    {"p": 3, "m": 4, "n": 79, "k": 2, "d": 4, "A_d": 130, "gamma": 2.65, "has_puncture_json": False},
    {"p": 3, "m": 4, "n": 72, "k": 9, "d": 3, "A_d": 648, "gamma": 1.89, "has_puncture_json": True},
    {"p": 3, "m": 5, "n": 230, "k": 13, "d": 6, "A_d": 572, "gamma": 1.60, "has_puncture_json": True},
    {"p": 3, "m": 5, "n": 215, "k": 28, "d": 5, "A_d": 1104, "gamma": 1.27, "has_puncture_json": True},
    {"p": 3, "m": 5, "n": 206, "k": 37, "d": 4, "A_d": 880, "gamma": 1.24, "has_puncture_json": True},
    {"p": 3, "m": 5, "n": 200, "k": 43, "d": 3, "A_d": 1700, "gamma": 1.40, "has_puncture_json": True},
    {"p": 3, "m": 6, "n": 690, "k": 39, "d": 5, "A_d": 1128, "gamma": 1.79, "has_puncture_json": True},
    {"p": 3, "m": 6, "n": 667, "k": 62, "d": 4, "A_d": 3972, "gamma": 1.71, "has_puncture_json": True},
    {"p": 5, "m": 2, "n": 24, "k": 1, "d": 3, "A_d": 96, "gamma": 2.89, "has_puncture_json": False},
    {"p": 5, "m": 2, "n": 20, "k": 5, "d": 2, "A_d": 760, "gamma": 2.00, "has_puncture_json": True},
    {"p": 5, "m": 3, "n": 124, "k": 1, "d": 4, "A_d": 124, "gamma": 3.48, "has_puncture_json": False},
    {"p": 5, "m": 3, "n": 112, "k": 13, "d": 3, "A_d": 512, "gamma": 1.96, "has_puncture_json": True},
    {"p": 5, "m": 4, "n": 519, "k": 106, "d": 5, "A_d": 2180, "gamma": 0.99, "has_puncture_json": True},
]

# ---------------------------------------------------------------------------
# Table 1 -- asymptotic optimal yield gamma_0(p) and minimizer t_0(p).
# { p: (gamma_0, t_0) }
# ---------------------------------------------------------------------------
TABLE1_ASYMPTOTIC = {
    2: (0.67799, 0.27063),
    3: (0.63215, 0.27369),
    5: (0.55914, 0.27868),
    7: (0.50786, 0.28230),
    11: (0.44108, 0.28720),
    13: (0.41785, 0.28896),
    17: (0.38273, 0.29169),
    19: (0.36901, 0.29278),
    23: (0.34663, 0.29459),
}

# ---------------------------------------------------------------------------
# Table 2 -- smallest Manhattan-puncture codes with gamma < 1.
# Fields: p, m, w, n, k, d.  Invariant n + k = p^m. (p=2 row from [22].)
# These exercise the analytic (integer) engine: n=[m,>w]_p, k=[m,<=w]_p,
# d=Delta_p(m, rtilde, w).
# ---------------------------------------------------------------------------
TABLE2_SMALLEST_SUBLOG = [
    {"p": 2, "m": 58, "w": 14, "n": 288215893050995568, "k": 14483100716176, "d": 21700},
    {"p": 3, "m": 32, "w": 16, "n": 1852445880782154, "k": 574308069687, "d": 3510},
    {"p": 5, "m": 16, "w": 16, "n": 152186472515, "k": 401418110, "d": 429},
    {"p": 7, "m": 13, "w": 20, "n": 96448935471, "k": 440074936, "d": 231},
    {"p": 11, "m": 7, "w": 19, "n": 18874416, "k": 612755, "d": 35},
    {"p": 13, "m": 7, "w": 23, "n": 60848853, "k": 1899664, "d": 35},
    {"p": 17, "m": 4, "w": 17, "n": 77540, "k": 5981, "d": 15},
    {"p": 19, "m": 4, "w": 19, "n": 121470, "k": 8851, "d": 15},
    {"p": 23, "m": 1, "w": 5, "n": 17, "k": 6, "d": 3},
]

# ---------------------------------------------------------------------------
# Worked distillation example (Section 5): the [[519,106,5]]_5 code at
# delta_in = 1e-3 gives delta_out ~ 8e-18 at cost C = n/nbar_T ~ 7.4.
# ---------------------------------------------------------------------------
DISTILLATION_EXAMPLE = {
    "label": "[[519,106,5]]_5",
    "p": 5, "n": 519, "k": 106, "d": 5, "A_d": 2180,
    "delta_in": 1e-3,
    "delta_out_approx": 8e-18,   # order of magnitude
    "cost_C_approx": 7.4,        # +-0.5
}
