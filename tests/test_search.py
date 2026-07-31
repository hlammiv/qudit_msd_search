"""Tests for qmsd.search: the Manhattan sweep, the randomized search, and the top-level driver."""
import pytest

from qmsd.search import manhattan_sweep, random_search, search


def test_manhattan_sweep_returns_valid_sorted_codes():
    codes = manhattan_sweep(3, 6)  # qutrit, p^m = 729 (analytic -> fast, no matrices)
    assert len(codes) >= 2
    # every returned code is well-formed and gamma-sorted ascending
    gammas = [c.gamma for c in codes]
    assert gammas == sorted(gammas)
    for c in codes:
        assert c.d >= 2 and c.k >= 1 and c.n > c.k and c.d_certified


def test_random_search_is_deterministic_and_valid():
    a = random_search(5, 2, trials=80, seed=7)
    b = random_search(5, 2, trials=80, seed=7)
    assert [c.label for c in a] == [c.label for c in b]   # reproducible under fixed seed
    assert len(a) >= 1
    # the search must NEVER emit an invalid or uncertified code
    for c in a:
        assert c.full_rank and c.d_certified and c.d >= 2 and c.n > c.k
        assert c.n + c.k == c.p ** c.m


def test_random_search_target_k_respected():
    codes = random_search(5, 2, trials=120, seed=3, target_k=5)
    assert codes  # finds at least one full-rank [[20,5,d]] code
    assert all(c.k == 5 and c.n == 20 for c in codes)


def test_random_search_parallel_valid_and_deterministic():
    # parallel results are reproducible for a fixed (seed, n_jobs, trials)
    a = random_search(5, 2, trials=200, seed=1, target_k=5, n_jobs=4)
    b = random_search(5, 2, trials=200, seed=1, target_k=5, n_jobs=4)
    assert [c.label for c in a] == [c.label for c in b]
    assert a, "parallel search found no codes"
    # never emits an invalid/uncertified code, even across worker processes
    for c in a:
        assert c.full_rank and c.d_certified and c.d >= 2 and c.n + c.k == c.p ** c.m
    # and it still finds the [[20,5,2]]_5 code (k=5, d=2)
    assert any(c.n == 20 and c.k == 5 and c.d == 2 for c in a)


def test_search_driver_structure():
    # restrict to m=2 so the driver test stays fast (explicit search now reaches large m)
    res = search(5, m_values=[2], trials_per_m=60, seed=0, top=5)
    assert res["p"] == 5
    assert res["scanned"]["n_candidates"] >= 1
    assert 1 <= len(res["best_by_gamma"]) <= 5
    assert 1 <= len(res["best_by_cost"]) <= 5
    # best_by_gamma is sorted ascending in gamma
    g = [c.gamma for c in res["best_by_gamma"]]
    assert g == sorted(g)


def test_capset_climb_reaches_high_distance():
    # The cap-set climb reaches d>=3 where uniform random stalls at d=2 (p=3,m=4,k=9),
    # reconstructing the paper's [[72,9,3]]_3. Deterministic for this seed/budget.
    codes = random_search(3, 4, trials=15, seed=0, target_k=9, sampler="capset_climb")
    assert any(c.d >= 3 for c in codes), "capset_climb did not reach d>=3"
    assert any(c.n == 72 and c.k == 9 and c.d == 3 for c in codes)
    for c in codes:
        assert c.full_rank and c.d_certified and c.n + c.k == 3 ** 4


def test_capset_modes_deterministic_and_parallel_valid():
    a = random_search(3, 4, trials=10, seed=1, target_k=9, sampler="capset_climb")
    b = random_search(3, 4, trials=10, seed=1, target_k=9, sampler="capset_climb")
    assert [c.label for c in a] == [c.label for c in b]   # deterministic per config
    # cap sampler (seed only) and parallel cap search stay valid
    cs = random_search(3, 4, trials=30, seed=0, target_k=9, sampler="capset")
    par = random_search(3, 4, trials=12, seed=1, target_k=9, sampler="capset_climb", n_jobs=3)
    for c in list(cs) + list(par):
        assert c.full_rank and c.d_certified and c.d >= 2 and c.n + c.k == 3 ** 4


def test_arc_climb_valid_deterministic_and_reports_Ad():
    # arc_climb ranks candidates by the (distance, -A_d) surrogate; in the small-dual regime
    # (p=3,m=4) the exact A_d engine is feasible, so found codes carry an exact A_d.
    a = random_search(3, 4, trials=4, seed=0, target_k=9, sampler="arc_climb",
                      climb_steps=8, swap_tries=4)
    b = random_search(3, 4, trials=4, seed=0, target_k=9, sampler="arc_climb",
                      climb_steps=8, swap_tries=4)
    assert [c.label for c in a] == [c.label for c in b]   # deterministic per config
    assert a
    for c in a:
        assert c.full_rank and c.d_certified and c.d >= 2 and c.n + c.k == 3 ** 4
    assert any(c.A_d is not None for c in a)   # the A_d surrogate was actually computed


def test_invalid_sampler_raises():
    with pytest.raises(ValueError):
        random_search(3, 4, trials=2, sampler="bogus")
