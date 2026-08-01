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


def test_plane_spread_reproduces_high_distance_paper_code():
    # plane_spread (cap + no 4 coplanar) reproduces the paper's [[230,13,6]]_3 (d=6) -- the
    # highest-distance qutrit code -- which plain caps never reach (they top out at d=5, and
    # 240k uniform samples never exceed d=3). Every built size-13 plane-spread cap is d=6.
    codes = [c for c in random_search(3, 5, trials=15, seed=0, target_k=13, sampler="plane_spread")
             if c.k == 13]
    assert codes, "plane_spread built no size-13 sets"
    assert any(c.n == 230 and c.k == 13 and c.d == 6 for c in codes), \
        "plane_spread did not reproduce [[230,13,6]] d=6"
    for c in codes:
        assert c.full_rank and c.d_certified and c.d >= 2


def test_near_cap_builds_where_strict_caps_stall():
    # near_cap (cap relaxed by max_triples) builds valid codes at k=43 in AG(5,3), where strict
    # caps stall (near the max-cap = 45). At scale it reproduces [[200,43,3]]_3 (d=3, ~1/2000);
    # here we just confirm it builds valid, full-rank, distance-certified codes.
    codes = [c for c in random_search(3, 5, trials=8, seed=0, target_k=43,
                                      sampler="near_cap", max_triples=45) if c.k == 43]
    assert codes, "near_cap built no size-43 sets"
    for c in codes:
        assert c.full_rank and c.d_certified and c.d >= 2 and c.n + c.k == 3 ** 5


def test_flat_spread_unifies_and_auto_selects_order():
    # The unified flat_spread sampler auto-picks the max feasible arc order for the given k,
    # monotonically stepping cap-order 3 -> 2 -> 1 -> near-cap fallback (0) as k grows. This is
    # the single sampler that subsumes capset (order 1), plane_spread (order 2), near_cap (0).
    from qmsd.sampling import max_feasible_order, all_points
    allpts = all_points(5, 3)
    assert max_feasible_order(5, 3, 6, allpts) == 3       # order-3 arc: the d=7 regime
    assert max_feasible_order(5, 3, 13, allpts) == 2      # order-2 == plane_spread ([[230,13,6]])
    assert max_feasible_order(5, 3, 28, allpts) == 1      # order-1 == cap
    assert max_feasible_order(5, 3, 43, allpts) == 0      # cap stalls -> near-cap fallback
    # and it builds valid, distance-certified codes (k=43 -> near-cap fallback, low d)
    codes = [c for c in random_search(3, 5, trials=6, seed=0, target_k=43, sampler="flat_spread")
             if c.k == 43]
    assert codes and all(c.full_rank and c.d_certified and c.d >= 2 for c in codes)


def test_invalid_sampler_raises():
    with pytest.raises(ValueError):
        random_search(3, 4, trials=2, sampler="bogus")
