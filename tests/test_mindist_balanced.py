"""Tests for the memory-rebalanced MITM minimum distance (qmsd.mindist_balanced).

The contract: ``min_dependent_columns_balanced`` returns the EXACT same value as the
verified ``mindist.min_dependent_columns`` on every small/feasible code, while the d=6
step uses the rebalanced (a=3, left-coefficient-fixed) layout.  We pin:

  1. equality vs mindist on a random fuzz over several primes (covers d = 1..6),
  2. a constructed MDS (Vandermonde) code whose distance is exactly 6 -- forcing the
     rebalanced d=6 path -- both serial and parallel (n_jobs > 1),
  3. the RAM-projection guard refuses (MemoryError) when the projection exceeds the
     allowed fraction, and the projection helpers report the documented sizes.
"""
import itertools
import random

import galois
import numpy as np
import pytest

from qmsd.mindist import min_dependent_columns
from qmsd.mindist_balanced import (
    _weight6_exists_balanced,
    available_ram_bytes,
    left_table_bytes,
    left_table_entries,
    min_dependent_columns_balanced,
)
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.oracle import load_oracle


def _brute_mdc(H, p, dmax):
    GF = galois.GF(p)
    r, n = H.shape
    for d in range(1, dmax + 1):
        for cols in itertools.combinations(range(n), d):
            if np.linalg.matrix_rank(GF(H[:, list(cols)].astype(int) % p)) < d:
                return d
    return None


def _vandermonde(p, r, points):
    """r x len(points) Vandermonde over F_p: rows x^0..x^(r-1). MDS -> min dep cols = r+1."""
    return np.array([[pow(x, i, p) for x in points] for i in range(r)], dtype=np.int64) % p


def test_balanced_matches_mindist_fuzz():
    """Balanced == mindist on random matrices over several primes (d = 1..6)."""
    rng = random.Random(7)
    nprng = np.random.RandomState(11)
    checked = 0
    for _ in range(120):
        p = rng.choice([2, 3, 5, 7])
        r = rng.randint(2, 5)
        n = rng.randint(r + 1, r + 8)
        H = nprng.randint(0, p, size=(r, n)).astype(np.int64)
        if rng.random() < 0.3:
            H[:, rng.randrange(n)] = 0
        if rng.random() < 0.3 and n >= 2:
            a, b = rng.sample(range(n), 2)
            H[:, b] = (rng.randint(1, p - 1) * H[:, a]) % p
        ref = min_dependent_columns(H, p, d_max=6)
        got = min_dependent_columns_balanced(H, p, d_max=6, n_jobs=1)
        assert got == ref, f"p={p} r={r} n={n}: mindist={ref} balanced={got}\n{H}"
        checked += 1
    assert checked == 120


def test_balanced_distance_six_forces_rebalanced_path():
    """A 5xN MDS code has min distance 6 -> exercises the rebalanced d=6 engine."""
    p, r = 7, 5
    H = _vandermonde(p, r, points=[0, 1, 2, 3, 4, 5])  # 6 distinct points in F_7
    # ground truth: MDS, every 5 cols independent, all 6 dependent -> distance 6.
    assert _brute_mdc(H, p, dmax=6) == 6
    assert min_dependent_columns(H, p, d_max=6) == 6
    assert min_dependent_columns_balanced(H, p, d_max=6, n_jobs=1) == 6
    # parallel right stream must agree.
    assert min_dependent_columns_balanced(H, p, d_max=6, n_jobs=4) == 6


def test_balanced_no_weight_six_when_seventh_needed():
    """7xN MDS: every 6 cols independent, so the d<=6 search finds nothing -> ValueError."""
    p, r = 11, 6
    H = _vandermonde(p, r, points=list(range(7)))  # MDS -> min dep cols = 7 > cap
    with pytest.raises(ValueError):
        min_dependent_columns_balanced(H, p, d_max=6, n_jobs=1)


def test_edge_cases():
    # zero column -> distance 1
    assert min_dependent_columns_balanced(np.array([[0, 1, 2], [0, 2, 1]]), 3) == 1
    # proportional columns -> distance 2
    assert min_dependent_columns_balanced(np.array([[1, 2, 1], [1, 2, 0]]), 3) == 2
    # empty parity check -> distance 1
    assert min_dependent_columns_balanced(np.empty((0, 4), dtype=int), 5) == 1


def test_projection_helpers():
    # documented pin numbers: (p-1)^2 * C(237,3) entries, ~32 bytes each.
    from math import comb

    assert left_table_entries(17, 3, 237) == 16 ** 2 * comb(237, 3)
    assert left_table_entries(17, 3, 237) == 256 * comb(237, 3)
    assert left_table_bytes(17, 3, 237) == left_table_entries(17, 3, 237) * 32


def test_ram_guard_refuses():
    """With an absurdly small ram_fraction the d=6 step must refuse (MemoryError), logged."""
    p, r = 7, 5
    H = _vandermonde(p, r, points=[0, 1, 2, 3, 4, 5])
    from qmsd.mindist import _powers

    powers = _powers(p, r)
    assert available_ram_bytes() >= 0
    with pytest.raises(MemoryError):
        _weight6_exists_balanced(H, p, powers, n_jobs=1, ram_fraction=1e-15)


def test_balanced_matches_mindist_small_oracle():
    """On the small oracle codes (feasible d<=6) balanced reproduces the certified distance."""
    for oc in load_oracle():
        built = build_triorthogonal_code(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed)
        G0 = np.asarray(built["X_stab"]) % oc.p
        r, n = G0.shape
        # keep CI fast: only the genuinely small duals (small redundancy r).
        if r > 12 or n > 250:
            continue
        ref = min_dependent_columns(G0, oc.p, d_max=6)
        got = min_dependent_columns_balanced(G0, oc.p, d_max=6, n_jobs=1)
        assert got == ref == oc.d, oc.label
