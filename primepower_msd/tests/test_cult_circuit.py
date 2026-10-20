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
