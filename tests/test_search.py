"""Tests for qmsd.search: the Manhattan sweep, the randomized search, and the top-level driver."""
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
