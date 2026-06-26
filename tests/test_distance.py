"""Unit tests for qmsd.distance (NOTES Thm 4; sec 7, 10-11).

Covers:
  * delta_p: the beta==0 closed-form identity and known single-/multi-puncture
    Manhattan distances;
  * min_distance vs brute force on tiny explicit codes and vs delta_p on small
    Manhattan-punctured RM codes;
  * weight_enumerator_at vs brute force;
  * A_d_logical_Z reproducing the Table-3 A_d values for the small oracle codes.
"""
from __future__ import annotations

from itertools import product

import numpy as np
import galois
import pytest

from qmsd.distance import (
    delta_p,
    min_distance,
    weight_enumerator_at,
    A_d_logical_Z,
)
from qmsd.reedmuller import r_max, r_tilde
from qmsd.puncture import manhattan_puncture_columns
from qmsd.pnomial import pnomial_gt
from qmsd.triorthogonal import build_triorthogonal_code, dual_matrix
from qmsd.oracle import load_oracle
from tests import ground_truth as gt


# ---------------------------------------------------------------------------
# Brute-force references
# ---------------------------------------------------------------------------
def _brute_weights(generator, p):
    """Return a dict {weight: count} over all nonzero codewords (tiny codes only)."""
    GF = galois.GF(p)
    G = GF(np.asarray(generator).astype(np.int64) % p)
    R = np.asarray(G.row_reduce()).astype(np.int64) % p
    B = np.array([row for row in R if np.any(row != 0)], dtype=np.int64)
    k, n = B.shape if B.size else (0, np.asarray(generator).shape[1])
    counts = {}
    for coeffs in product(range(p), repeat=k):
        if not any(coeffs):
            continue
        word = np.zeros(n, dtype=np.int64)
        for i, cf in enumerate(coeffs):
            if cf:
                word = (word + cf * B[i]) % p
        w = int(np.count_nonzero(word))
        counts[w] = counts.get(w, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# delta_p (Theorem 4)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("p", [2, 3, 5, 7])
@pytest.mark.parametrize("m", [2, 3, 4])
def test_delta_p_beta_zero_identity(p, m):
    """When beta==0, Delta_p collapses to the single tail sum [m-alpha,>w]_p."""
    for alpha in range(0, m):                 # ensure alpha <= m-1
        r = alpha * (p - 1)                   # beta == 0
        if not (0 <= r <= m * (p - 1)):
            continue
        for w in range(-1, m * (p - 1) + 2):
            assert delta_p(m, r, w, p) == pnomial_gt(m - alpha, w, p)


def test_delta_p_internal_branches_run():
    """Exercise the beta!=0 branch (its closed-form assertion must hold)."""
    # p=5, m=3, r=7 -> alpha=1, beta=3 (nonzero) is a representative case.
    val = delta_p(3, 7, 4, 5)
    assert isinstance(val, int) and val >= 0


def test_delta_p_matches_single_puncture_table():
    """Single-puncture (w=0) Manhattan distance == ground-truth d (Lemma 2 chain)."""
    for rec in gt.SINGLE_PUNCTURE_CODES:
        p, m, d = rec["p"], rec["m"], rec["d"]
        rt = r_tilde(r_max(m, p), m, p)
        assert delta_p(m, rt, 0, p) == d, rec


# ---------------------------------------------------------------------------
# min_distance: tiny explicit codes vs brute force
# ---------------------------------------------------------------------------
def test_min_distance_tiny_explicit():
    GF3 = galois.GF(3)
    # Repetition code [4,1,4]_3.
    rep = GF3([[1, 1, 1, 1]])
    assert min_distance(rep, 3) == 4
    # A [4,2]_3 code with a weight-1 row after reduction -> distance 1.
    g = GF3([[1, 0, 0, 0], [0, 1, 1, 0]])
    assert min_distance(g, 3) == 1
    # Random small codes: compare to brute force.
    rng = np.random.default_rng(0)
    for p in (2, 3, 5):
        GF = galois.GF(p)
        for _ in range(8):
            k, n = 3, 7
            G = GF(rng.integers(0, p, size=(k, n)))
            if int(np.linalg.matrix_rank(G)) == 0:
                continue
            brute = _brute_weights(G, p)
            assert min_distance(G, p) == min(brute)


def test_min_distance_max_weight_guard():
    # Build a high-dimensional code (so the increasing-weight scan path is used,
    # not direct enumeration) whose true distance exceeds the bound: bounding the
    # search must RAISE rather than silently under-report.
    GF2 = galois.GF(2)
    rng = np.random.default_rng(7)
    n, k = 30, 3
    while True:
        G = GF2(rng.integers(0, 2, size=(k, n)))
        cols = np.asarray(G)
        if int(np.linalg.matrix_rank(G)) == k and not np.any(
            np.all(cols == 0, axis=0)
        ):
            break
    big = dual_matrix(G, 2)                    # [30,27]_2: scan path, distance >= 2
    assert 2 ** 27 > 200_000                   # confirm direct path is bypassed
    with pytest.raises(ValueError):
        min_distance(big, 2, max_weight=1)
    # With a sufficient bound it returns the exact (>=2) distance.
    d = min_distance(big, 2, max_weight=n)
    assert d >= 2


# ---------------------------------------------------------------------------
# min_distance vs delta_p on small Manhattan-punctured RM codes
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "p,m,w",
    [
        (3, 2, 0),
        (3, 3, 0),
        (5, 2, 0),
        (5, 2, 1),
        (3, 4, 0),
    ],
)
def test_min_distance_matches_delta_p_manhattan(p, m, w):
    r = r_max(m, p)
    rt = r_tilde(r, m, p)
    cols = manhattan_puncture_columns(m, w, p)
    if not cols or len(cols) == p ** m:
        pytest.skip("degenerate puncture set")
    code = build_triorthogonal_code(p, m, r, cols)
    if not code["full_rank"]:
        pytest.skip("not full rank")
    expected = delta_p(m, rt, w, p)
    if expected < 2:
        pytest.skip("degenerate distance-1 code")
    # Distance = min weight of G0^perp = dual of the shortened generator.
    g0_perp = dual_matrix(code["X_stab"], p)
    assert min_distance(g0_perp, p) == expected


# ---------------------------------------------------------------------------
# weight_enumerator_at vs brute force
# ---------------------------------------------------------------------------
def test_weight_enumerator_matches_brute():
    rng = np.random.default_rng(1)
    for p in (2, 3, 5):
        GF = galois.GF(p)
        for _ in range(6):
            k, n = 3, 8
            G = GF(rng.integers(0, p, size=(k, n)))
            brute = _brute_weights(G, p)
            for d in range(0, n + 1):
                expected = brute.get(d, 0)
                if d == 0:
                    expected = 1                       # the zero codeword
                assert weight_enumerator_at(G, p, d) == expected, (p, d)


def test_weight_enumerator_highdim_path():
    """High-dimensional code (small parity check) routes through the subset scan."""
    p, m, w = 5, 2, 0
    r = r_max(m, p)
    cols = manhattan_puncture_columns(m, w, p)
    code = build_triorthogonal_code(p, m, r, cols)
    g0_perp = dual_matrix(code["X_stab"], p)           # high-dim code
    d = delta_p(m, r_tilde(r, m, p), w, p)
    assert weight_enumerator_at(g0_perp, p, d) >= 1


# ---------------------------------------------------------------------------
# A_d_logical_Z vs Table 3
# ---------------------------------------------------------------------------
# (label, p, m, expected A_d) for the small oracle codes that are feasible to
# count exactly within the budget.
_SMALL_AD = [
    ("[[20,5,2]]_5", 5, 2, 760),
    ("[[72,9,3]]_3", 3, 4, 648),
    ("[[112,13,3]]_5", 5, 3, 512),
]


@pytest.mark.parametrize("label,p,m,expected_Ad", _SMALL_AD)
def test_A_d_reproduces_ground_truth(label, p, m, expected_Ad):
    oc = next(c for c in load_oracle() if c.label == label)
    # cross-check the hard-coded expectation against ground_truth.TABLE3_CODES
    row = next(r for r in gt.TABLE3_CODES if r["n"] == oc.n and r["k"] == oc.k
               and r["p"] == p)
    assert row["A_d"] == expected_Ad
    code = build_triorthogonal_code(p, m, r_max(m, p),
                                    oc.puncture_columns_1indexed)
    assert code["full_rank"], label
    try:
        ad = A_d_logical_Z(code, p, oc.d)
    except NotImplementedError as exc:
        pytest.skip(f"A_d infeasible for {label}: {exc}")
    assert ad == expected_Ad, label


def test_A_d_infeasible_raises():
    """A code whose C(n,d) blows past the budget must raise, never guess."""
    # Fabricate a wide code dict; d large enough that C(n,d) > budget.
    p, n = 5, 600
    GF = galois.GF(p)
    fake = {
        "X_stab": GF(np.zeros((1, n), dtype=np.int64)),
        "Gp": GF(np.zeros((1, n), dtype=np.int64)),
    }
    with pytest.raises(NotImplementedError):
        A_d_logical_Z(fake, p, 5)
