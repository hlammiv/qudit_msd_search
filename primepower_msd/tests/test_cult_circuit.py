import pytest

pytest.importorskip("stim")

from primepower_msd.cult_circuit import ColorCult
from primepower_msd.zdsim import check_determinism


@pytest.mark.parametrize("d", [2, 4])
def test_d3_cultivation_detectors_and_fault_distance(d):
    cc = ColorCult(3, d, 1e-3, 1, 1, checks_per_round=2, rounds=1)
    assert check_determinism(cc.c) == []
    r = cc.single_fault_scan()
    assert r["undetected_harmful"] == 0


def test_mismatched_injection_is_caught():
    cc = ColorCult(3, 4, 0.0, sigma=1, inj=3, checks_per_round=1, rounds=1)
    assert check_determinism(cc.c) != []          # wrong state -> the T-check is not deterministic


def test_parallel_schedule_fault_distance_with_tick_idle():
    from primepower_msd.zdsim import add_tick_idle, check_determinism
    cc = ColorCult(3, 4, 1e-3, 1, 1, checks_per_round=2, rounds=1, parallel=True)
    assert check_determinism(cc.c, 64, 0) == []
    assert cc.single_fault_scan()["undetected_harmful"] == 0
    serial = ColorCult(3, 4, 1e-3, 1, 1, checks_per_round=2, rounds=1)
    cc.c = add_tick_idle(cc.c, 1e-3)
    assert check_determinism(cc.c, 64, 0) == []
    assert cc.single_fault_scan()["undetected_harmful"] == 0
    assert cc.c.nq == serial.c.nq + serial.n                  # the second cat register
    assert cc.c.n_meas > 0 and len(cc.c.ops) > len(serial.c.ops)


def test_hook_noise_final_round_boundary():
    """with fan-outs split into noisy legs the final round's hooks reach the ideal decode (fault distance 1 there);
    postselecting the final data syndrome (the next stage detects it) restores fault distance 2, also for the parallel
    schedule with a block final round."""
    from primepower_msd.zdsim import add_tick_idle
    cc = ColorCult(3, 4, 1e-3, 1, 1, checks_per_round=2, rounds=1, hooks=True)
    assert cc.single_fault_scan()["undetected_harmful"] > 0
    assert cc.single_fault_scan(final_postselect=True)["undetected_harmful"] == 0
    cp = ColorCult(3, 4, 1e-3, 1, 1, checks_per_round=2, rounds=1, hooks=True, parallel=True, final_block=True)
    cp.c = add_tick_idle(cp.c, 1e-3)
    assert check_determinism(cp.c, 64, 0) == []
    assert cp.single_fault_scan(final_postselect=True)["undetected_harmful"] == 0
