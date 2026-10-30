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

from qmsd.mindist import min_dependent_columns, _weight_d_exists, _powers
from qmsd.mindist_nb import (weight_d_exists_fast, weight6_exists_via_2_4,
                             min_distance_upto6_lowmem, HAVE_NUMBA)
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


def test_numba_kernel_matches_numpy_weight_d_exists():
    """The numba MITM kernel must return the IDENTICAL weight-d existence verdict as the numpy
    reference, per d, over a fuzz battery incl. p=7 (the large-p regime it exists to speed up)."""
    rng = random.Random(7)
    nprng = np.random.RandomState(3)
    checks = 0
    for _ in range(120):
        p = rng.choice([2, 3, 5, 7])
        r = rng.randint(2, 5)
        n = rng.randint(r + 2, r + 8)
        H = (nprng.randint(0, p, size=(r, n))).astype(np.int64) % p
        powers = _powers(p, r)
        for d in range(2, min(6, n) + 1):
            assert weight_d_exists_fast(H, p, d) == _weight_d_exists(H, p, d, powers)
            checks += 1
    assert checks >= 200


def test_2plus4_split_certifies_d6_without_a3_table():
    # The 2+4 split (a=2 left table + b=4 stream) certifies d=6 using only the a=2 table, so d=6
    # is reachable where the standard a=3,b=3 path OOMs. Correctness pinned against brute force.
    import itertools

    def brute_weight6(H, p):
        r, n = H.shape
        for cols in itertools.combinations(range(n), 6):
            M = H[:, list(cols)] % p
            for coef in itertools.product(range(p), repeat=6):
                if all(coef) and not np.any((np.array(coef) @ M.T) % p):
                    return True
        return False

    # hand-built min-distance-6 circuit over F_7: e0..e4 plus col5 = -(sum)
    p = 7
    I = np.eye(5, dtype=np.int64)
    H6 = np.column_stack([I, (-(I.sum(0)) % p).reshape(-1, 1)]) % p
    assert min_distance_upto6_lowmem(H6, p) == 6
    assert weight6_exists_via_2_4(H6, p) is True

    rng = random.Random(11)
    nprng = np.random.RandomState(5)
    for _ in range(40):
        p = rng.choice([2, 3, 5, 7])
        r = rng.randint(3, 6)
        n = rng.randint(6, 9)
        H = nprng.randint(0, p, size=(r, n)).astype(np.int64) % p
        assert weight6_exists_via_2_4(H, p) == brute_weight6(H, p)


def test_numba_is_the_active_path():
    # Guard: in this environment the fast kernel is what min_dependent_columns actually runs.
    assert HAVE_NUMBA


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


def test_weight6_witness_parallel_matches_serial():
    """The thread-parallel 2+4 weight-6 scan agrees with the serial kernel and returns a genuine
    weight-6 dependency, on random F_7 codes with and without weight-6 words."""
    import numpy as np
    from qmsd.mindist_nb import weight6_exists_via_2_4, weight6_witness_parallel, weight_d_exists_fast
    from qmsd.structured_ad import _rank
    rng = np.random.default_rng(11)
    p, seen = 7, {True: 0, False: 0}
    for _ in range(120):
        r, n = int(rng.integers(12, 20)), int(rng.integers(24, 32))
        H = rng.integers(0, p, size=(r, n))
        if any(weight_d_exists_fast(H, p, d) for d in (2, 3, 4, 5)):
            continue
        a = weight6_exists_via_2_4(H, p)
        w = weight6_witness_parallel(H, p, block=4)
        assert (w is not None) == a
        if w is not None:
            assert len(set(w)) == 6 and _rank(H[:, w].T.copy(), p) < 6
        seen[a] += 1
        if seen[True] >= 2 and seen[False] >= 5:
            break
    assert seen[True] >= 1 and seen[False] >= 1


def test_weight6_witness_3x3_and_census_match_brute_force():
    """3+3 witness and per-output census vs exhaustive 6-subsets on small random F_7 codes."""
    import itertools
    import numpy as np
    from qmsd.mindist_nb import weight6_witness_3x3, count_weight6_per_output, weight_d_exists_fast
    from qmsd.structured_ad import _rank, _null_space
    rng = np.random.default_rng(31)
    p, done = 7, 0
    while done < 2:
        r, n = int(rng.integers(9, 11)), int(rng.integers(17, 20))
        H = rng.integers(0, p, size=(r, n))
        if any(weight_d_exists_fast(H, p, d) for d in (2, 3, 4, 5)):
            continue
        H1 = rng.integers(0, p, size=(2, n))
        classes, fails = 0, np.zeros(2, dtype=int)
        for S in itertools.combinations(range(n), 6):
            M = H[:, S]
            if _rank(M.T.copy(), p) < 6:
                v = np.asarray(_null_space(M, p))[0]
                classes += 1
                fails += ((H1[:, S] @ v) % p != 0)
        w = weight6_witness_3x3(H, p)
        assert (w is not None) == (classes > 0)
        if w is not None:
            assert len(set(w)) == 6 and _rank(H[:, w].T.copy(), p) < 6
        L, ncls, per = count_weight6_per_output(H, H1, p)
        assert ncls == classes and (per == fails * (p - 1)).all()
        done += 1


def test_weight7_witness_3x4_matches_brute_force():
    import itertools
    import numpy as np
    from qmsd.mindist_nb import weight7_witness_3x4, weight_d_exists_fast, weight6_exists_via_2_4
    from qmsd.structured_ad import _rank
    rng = np.random.default_rng(41)
    p, seen = 7, {True: 0, False: 0}
    for _ in range(80):
        r, n = int(rng.integers(10, 13)), int(rng.integers(15, 18))
        H = rng.integers(0, p, size=(r, n))
        if any(weight_d_exists_fast(H, p, d) for d in (2, 3, 4, 5)) or weight6_exists_via_2_4(H, p):
            continue
        brute = any(_rank(H[:, list(S)].T.copy(), p) < 7 for S in itertools.combinations(range(n), 7))
        w = weight7_witness_3x4(H, p, block=3)
        assert (w is not None) == brute
        if w is not None:
            assert len(set(w)) == 7 and _rank(H[:, w].T.copy(), p) < 7
        seen[brute] += 1
        if seen[True] >= 1 and seen[False] >= 2:
            break
    assert seen[True] >= 1 and seen[False] >= 1


def test_weight_d_exists_fast_rejects_d_above_six():
    """The 3+3 split cannot certify d >= 7; the kernel must refuse rather than scan weight 6."""
    import pytest
    from qmsd.mindist_nb import weight_d_exists_fast
    H = np.array([[1, 1, 3, 3, 1, 3, 1]], dtype=np.int64)  # 7 columns, one weight-7 dependency over F_7
    with pytest.raises(ValueError):
        weight_d_exists_fast(H, 7, 7)


def test_weight7_census_matches_null_space_enumeration():
    """3+4 weight-7 census against exhaustive enumeration of ker(G0) on small codes with d >= 7
    (nullity 4: 7^4 words), per-output classification included."""
    import itertools
    from qmsd.mindist_nb import count_weight7_per_output, weight_d_exists_fast, weight6_exists_via_2_4
    from qmsd.structured_ad import _null_space
    p = 7
    rng = np.random.default_rng(11)
    checked = 0
    for attempt in range(300):
        r, n = 12, 16
        H = rng.integers(0, p, size=(r, n))
        if np.linalg.matrix_rank(H) < r:
            continue
        if any(weight_d_exists_fast(H, p, d) for d in (2, 3, 4, 5)) or weight6_exists_via_2_4(H, p):
            continue
        H1 = rng.integers(0, p, size=(2, n))
        ns = np.asarray(_null_space(H, p)) % p
        assert ns.shape[0] == n - r
        coefs = np.array(list(itertools.product(range(p), repeat=ns.shape[0])))[1:]
        W = (coefs @ ns) % p
        w7 = W[(W != 0).sum(axis=1) == 7]
        if len(w7) == 0:
            continue
        brute_per_output = ((w7 @ H1.T) % p != 0).sum(axis=0)
        L7, n_classes, per_output = count_weight7_per_output(H, H1, p)
        assert n_classes * (p - 1) == len(w7)
        assert list(per_output) == list(brute_per_output)
        checked += 1
        if checked == 2:
            break
    assert checked >= 1


# ----------------------------------------------------------------------------------------------
# Bucket-indexed right scans (method="index") vs the original searchsorted kernels
# ----------------------------------------------------------------------------------------------

def _random_instances(seed, count, p_choices=(3, 5, 7, 11), r_range=(8, 14), n_range=(14, 24),
                      nullity=None):
    """Random (p, H, H1) instances; ``nullity`` restricts n - r (small nullity -> larger d)."""
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < count:
        p = int(rng.choice(p_choices))
        r = int(rng.integers(r_range[0], r_range[1] + 1))
        n = int(rng.integers(n_range[0], n_range[1] + 1))
        if nullity is not None:
            n = r + int(rng.integers(nullity[0], nullity[1] + 1))
        H = rng.integers(0, p, size=(r, n)).astype(np.int64)
        H1 = rng.integers(0, p, size=(2, n)).astype(np.int64)
        out.append((p, H, H1))
    return out


def _min_distance_reference(H, p):
    """Exact d in 1..7 or 8 (meaning >= 8) using ONLY the original searchsorted kernels."""
    from qmsd.mindist_nb import weight7_witness_3x4
    if (H == 0).all(axis=0).any():
        return 1
    for d in (2, 3, 4, 5, 6):
        if weight_d_exists_fast(H, p, d, method="searchsorted"):
            return d
    return 7 if weight7_witness_3x4(H, p, method="searchsorted") is not None else 8


def _is_dependency(H, p, cols, w):
    from qmsd.structured_ad import _rank
    return cols is not None and len(set(cols)) == w and _rank(H[:, list(cols)].T.copy(), p) < w


def test_left_index_brackets_every_code():
    """bstart[beta] .. bstart[beta+1] is exactly the run of codes with top-K bits == beta."""
    from qmsd.mindist_nb import _left_index
    rng = np.random.default_rng(5)
    for N in (0, 1, 2, 7, 100, 5000):
        codes = np.sort(rng.integers(0, 2 ** 63, size=N).astype(np.uint64) * np.uint64(2)
                        + rng.integers(0, 2, size=N).astype(np.uint64))
        bstart, shift = _left_index(codes)
        K = 64 - int(shift)
        assert bstart.shape[0] == 2 ** K + 1 and bstart[0] == 0 and bstart[-1] == N
        tops = (codes >> np.uint64(shift)).astype(np.int64)
        for beta in range(2 ** K):
            lo, hi = int(bstart[beta]), int(bstart[beta + 1])
            assert lo == np.searchsorted(tops, beta, side="left")
            assert hi == np.searchsorted(tops, beta, side="right")


def test_index_method_matches_searchsorted_weight_d_exists():
    """weight_d_exists_fast: identical verdict for every d (both b=3 splits go through the
    index) on random instances over p in {3, 5, 7, 11}, r 8-14, n 14-24."""
    checks = 0
    for p, H, _ in _random_instances(101, 60):
        for d in (2, 3, 4, 5, 6):
            a = weight_d_exists_fast(H, p, d, method="searchsorted")
            b = weight_d_exists_fast(H, p, d, method="index")
            assert a == b, (p, H.shape, d, a, b)
            checks += 1
            if a:
                break
    assert checks >= 60


def test_index_method_matches_searchsorted_parallel_scans():
    """weight5_exists_parallel, weight6_exists_via_2_4, weight6_witness_parallel,
    weight6_witness_3x3, weight7_witness_3x4: identical booleans, and every witness returned by
    the index path is a genuine dependency of the right weight."""
    from qmsd.mindist_nb import (weight5_exists_parallel, weight6_witness_parallel,
                                 weight6_witness_3x3, weight7_witness_3x4)
    seen = {5: 0, 6: 0, 7: 0, 8: 0}
    for p, H, _ in _random_instances(202, 40, nullity=(2, 6)):
        d = _min_distance_reference(H, p)
        if d < 5:
            continue
        seen[d] += 1
        a5 = weight5_exists_parallel(H, p, method="searchsorted")
        b5 = weight5_exists_parallel(H, p, method="index")
        assert a5 == b5 == (d == 5)
        if d == 5:
            continue
        for fn, kw in ((weight6_exists_via_2_4, {}), (weight6_witness_parallel, {"block": 4}),
                       (weight6_witness_3x3, {})):
            a = fn(H, p, method="searchsorted", **kw)
            b = fn(H, p, method="index", **kw)
            assert (a is not None and a is not False) == (b is not None and b is not False) == (d == 6), fn
            if isinstance(b, list):
                assert _is_dependency(H, p, b, 6)
        if d == 6 or p > 7:
            continue                     # 3+4 at p = 11 is C(n,4) 10^4 probes per instance
        a = weight7_witness_3x4(H, p, block=3, method="searchsorted")
        b = weight7_witness_3x4(H, p, block=3, method="index")
        assert (a is not None) == (b is not None) == (d == 7)
        if b is not None:
            assert _is_dependency(H, p, b, 7)
    assert seen[5] >= 1 and seen[6] >= 1 and seen[7] + seen[8] >= 1, seen


def test_index_method_matches_searchsorted_census_kernels():
    """The per-leading-column census arrays (fails, found) of the 2+3, 3+3 and 3+4 kernels are
    identical element-wise, and the public counters agree (incl. the 10/20/35-fold assertions)."""
    from qmsd.mindist_nb import (_powers_u64, _build_left2_coef, _build_left3_norm_coef,
                                 _left_index, _census_weight5_per_output, _census_weight5_idx,
                                 _census_weight6_3x3, _census_weight6_idx, _census_weight7_3x4,
                                 _census_weight7_idx, count_weight5_per_output,
                                 count_weight6_per_output, count_weight7_per_output)
    done = {5: 0, 6: 0, 7: 0}
    for p, H, H1 in _random_instances(303, 60, nullity=(3, 6)):
        r, n = H.shape
        H = H % p
        powers = _powers_u64(p, r)
        # kernel-level: identical arrays whatever the distance (hits are just verified matches)
        lc, ls, lf = _build_left2_coef(H, p, powers)
        bs, sh = _left_index(lc)
        f_old, g_old = _census_weight5_per_output(H, H1, p, powers, lc, ls, lf)
        f_new, g_new = _census_weight5_idx(H, H1, p, powers, lc, ls, lf, bs, sh, 0, n)
        assert np.array_equal(f_old, f_new) and np.array_equal(g_old, g_new)
        lc, ls, lf = _build_left3_norm_coef(H, p, powers)
        bs, sh = _left_index(lc)
        f_old, g_old = _census_weight6_3x3(H, H1, p, powers, lc, ls, lf)
        f_new, g_new = _census_weight6_idx(H, H1, p, powers, lc, ls, lf, bs, sh, 0, n)
        assert np.array_equal(f_old, f_new) and np.array_equal(g_old, g_new)
        if p <= 7:
            f_old, g_old = _census_weight7_3x4(H, H1, p, powers, lc, ls, lf, 0, n)
            f_new, g_new = _census_weight7_idx(H, H1, p, powers, lc, ls, lf, bs, sh, 0, n)
            assert np.array_equal(f_old, f_new) and np.array_equal(g_old, g_new)
        # public counters on instances whose distance admits them
        d = _min_distance_reference(H, p)
        if d >= 5:
            a = count_weight5_per_output(H, H1, p, method="searchsorted")
            b = count_weight5_per_output(H, H1, p, method="index")
            assert a[:2] == b[:2] and np.array_equal(a[2], b[2])
            done[5] += 1
        if d >= 6:
            a = count_weight6_per_output(H, H1, p, method="searchsorted")
            b = count_weight6_per_output(H, H1, p, method="index")
            assert a[:2] == b[:2] and np.array_equal(a[2], b[2])
            done[6] += 1
        if d >= 7 and p <= 7:
            a = count_weight7_per_output(H, H1, p, method="searchsorted")
            b = count_weight7_per_output(H, H1, p, method="index")
            assert a[:2] == b[:2] and np.array_equal(a[2], b[2])
            done[7] += 1
    assert done[5] >= 3 and done[6] >= 2 and done[7] >= 1, done


def test_index_method_rejects_unknown_method():
    H = np.eye(5, dtype=np.int64)
    with pytest.raises(ValueError):
        weight_d_exists_fast(H, 7, 5, method="bogus")
