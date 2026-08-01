"""Tests for the meet-in-the-middle exact minimum-distance routine (qmsd.mindist).

Correctness is pinned two ways: (1) an adversarial fuzz against brute-force column-
dependency search on random matrices, and (2) certifying ALL 10 published oracle codes'
minimum distances -- including the large blocks the naive scan could not finish.
"""
import itertools
import random

import numpy as np
import galois
import pytest

from qmsd.mindist import min_dependent_columns
from qmsd.oracle import load_oracle
from qmsd.triorthogonal import build_triorthogonal_code


def _brute_mdc(H, p, dmax):
    """Ground truth: smallest set of F_p-linearly-dependent columns, by exhaustion."""
    GF = galois.GF(p)
    r, n = H.shape
    for d in range(1, dmax + 1):
        for cols in itertools.combinations(range(n), d):
            if np.linalg.matrix_rank(GF(H[:, list(cols)].astype(int) % p)) < d:
                return d
    return None


def test_mitm_matches_bruteforce_fuzz():
    """MITM == brute force on random matrices over several primes (incl. d=1,2 edge cases)."""
    rng = random.Random(2024)
    nprng = np.random.RandomState(99)
    checked = 0
    for _ in range(150):
        p = rng.choice([2, 3, 5, 7])
        r = rng.randint(2, 5)
        n = rng.randint(r + 1, r + 9)
        H = nprng.randint(0, p, size=(r, n)).astype(np.int64)
        # inject zero / proportional columns to exercise the d=1 and d=2 paths
        if rng.random() < 0.3:
            H[:, rng.randrange(n)] = 0
        if rng.random() < 0.3 and n >= 2:
            a, b = rng.sample(range(n), 2)
            H[:, b] = (rng.randint(1, p - 1) * H[:, a]) % p
        bf = _brute_mdc(H, p, dmax=min(6, n))
        if bf is None:
            continue
        mm = min_dependent_columns(H, p, d_max=6)
        assert mm == bf, f"p={p} r={r} n={n}: brute={bf} mitm={mm}\n{H}"
        checked += 1
    assert checked >= 80


def test_mitm_edge_cases():
    # a zero column -> distance 1
    H = np.array([[0, 1, 2], [0, 2, 1]])
    assert min_dependent_columns(H, 3) == 1
    # two proportional columns, none zero -> distance 2
    H = np.array([[1, 2, 1], [1, 2, 0]])  # col0 and col1 are 2x each other over F_3
    assert min_dependent_columns(H, 3) == 2
    # empty parity check (r=0): every single column is a codeword -> distance 1
    assert min_dependent_columns(np.empty((0, 4), dtype=int), 5) == 1


# All 10 oracle codes: certify the quantum distance d(G0^perp) = min dependent columns
# of G0 (= X_stab). This is the new capability -- the large blocks included.
_ORACLE = list(load_oracle())


@pytest.mark.parametrize("oc", _ORACLE, ids=[oc.label for oc in _ORACLE])
def test_mitm_certifies_oracle_distance(oc):
    built = build_triorthogonal_code(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed)
    G0 = np.asarray(built["X_stab"]) % oc.p
    assert min_dependent_columns(G0, oc.p, d_max=6) == oc.d


def test_hash_encoder_handles_overflow_redundancy():
    # The 64-bit polynomial-hash syndrome match handles redundancy where the OLD p-adic int64
    # radix overflowed (p**r > 2**63). Here p=7, r=25 -> 7**25 = 1.3e21 >> int64: the old
    # encoder raised OverflowError. A planted weight-3 codeword is now found exactly (== brute
    # force), proving the hash + rank-verify path stays certified.
    import itertools
    from qmsd.mindist import _fp_rank
    rng = np.random.default_rng(1)
    p, r, n = 7, 25, 45
    H = rng.integers(0, p, (r, n))
    H[:, 2] = (-(H[:, 0] + H[:, 1])) % p          # weight-3 dependency on columns {0,1,2}
    assert p ** r > 2 ** 63                        # would have overflowed the old int64 radix

    def brute():
        for dd in range(1, 7):
            for c in itertools.combinations(range(n), dd):
                if _fp_rank(H[:, list(c)].T % p, p) < dd:
                    return dd
        return None

    assert min_dependent_columns(H, p, d_max=6) == brute() == 3
