import pytest

from primepower_msd.footprint import Host, cheapest, logical_qudits, pipelines
from primepower_msd.msd_planner import STAGES_D4


def _stage(name):
    return next(s for s in STAGES_D4 if s[0] == name)


def test_logical_qudits():
    assert logical_qudits(_stage("O3 [[7,1,3]] s=1")) == 1 + 6 + 1
    assert logical_qudits(_stage("O3 [[15,7,3]]+K8")) == 28 + 8 + 1
    assert logical_qudits(_stage("O7 [[64,26,7]]+2080")) == 2080 + 38 + 1


def test_host_scaling():
    h = Host(p_th=5e-3, A=0.1)
    assert h.p_logical(5e-3, 11) == pytest.approx(0.1)
    assert h.p_logical(1e-3, 5) == pytest.approx(0.1 * 0.2 ** 3)
    assert h.qudits(7) == 98


def test_cheapest_pipeline_is_small_round():
    r = cheapest(1e-3, 1e-10)
    assert r is not None and "[[7,1,3]]" in r["pipeline"] and r["logical"] == 8
    assert r["error"] <= 1e-10 and r["physical"] == 8 * 2 * r["D"] ** 2 + 56
    # a better host code needs a smaller distance and fewer qudits
    r2 = cheapest(1e-3, 1e-10, host=Host(p_th=1e-2))
    assert r2["D"] < r["D"] and r2["physical"] < r["physical"]
    # the Clifford floor a_cliff p_L(D) < target binds: no candidate below the first admissible distance
    assert all(c["D"] >= r["D"] for c in pipelines(1e-3, 1e-10))


def test_tighter_target_costs_more():
    a, b = cheapest(5e-4, 1e-10), cheapest(5e-4, 1e-15)
    assert b["physical"] > a["physical"] and b["error"] <= 1e-15
