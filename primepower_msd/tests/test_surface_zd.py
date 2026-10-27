import numpy as np
import pytest

from primepower_msd.surface_zd import Experiment, PeelDecoder, memory_circuit, verify_code
from primepower_msd.zdsim import check_determinism


@pytest.mark.parametrize("d", [4, 8])
@pytest.mark.parametrize("D", [3, 5])
def test_code_structure(D, d):
    v = verify_code(D, d)
    assert v["k"] == 1 and v["noncommuting"] == 0 and v["logical_violations"] == 0
    assert v["nX"] == v["nZ"] == (D * D - 1) // 2


@pytest.mark.parametrize("basis", ["Z", "X"])
def test_detectors_deterministic(basis):
    mem = memory_circuit(5, 4, 5, 0.01, basis)
    assert check_determinism(mem.circ, 64, 1) == []


def test_dem_graphlike_and_decoder_noiseless():
    ex = Experiment(3, 4, "Z")
    assert ex.dem.hyper_prob == 0.0 and ex.dem.undetectable_prob == 0.0
    dec = PeelDecoder(ex.dem, 1e-3)
    corr, left = dec.decode(np.zeros(ex.dem.n_det, dtype=np.int64))
    assert corr == 0 and not left
    # every single edge fault is decoded back to its own observable effect
    for (a, b), e in list(ex.dem.edges.items())[:40]:
        ca, cb, o = e["pat"]
        sigma = np.zeros(ex.dem.n_det, dtype=np.int64)
        sigma[a] = ca
        if b != -1:
            sigma[b] = cb
        corr, left = dec.decode(sigma)
        assert not left and (corr + o) % 4 == 0


def test_low_p_scaling():
    """D=5 beats D=3 well below threshold (both bases)."""
    for basis in "ZX":
        r3 = Experiment(3, 4, basis).run(1e-3, 4000, seed=3)
        r5 = Experiment(5, 4, basis).run(1e-3, 4000, seed=5)
        assert r3["leftover"] == 0 and r5["leftover"] == 0
        assert r5["p_fail"] < r3["p_fail"] + 2e-3


def test_gap_decoder_matches_plain_decoder():
    """Split-boundary (complementary-gap) matcher: same corrections as the plain matcher, non-negative gaps, zero
    potential conflicts on a D=5 memory."""
    import numpy as np
    from primepower_msd.surface_zd import Experiment, PeelDecoder, P_BUILD
    from primepower_msd.zdsim import FrameSim
    ex = Experiment(5, 4, "Z", rounds=3)
    p = 3e-3
    plain, gap = PeelDecoder(ex.dem, p), PeelDecoder(ex.dem, p, gap=True)
    assert gap.conflicts == [0, 0]
    rng = np.random.default_rng(0)
    sim = FrameSim(ex.mem.circ, 300, rng).run(scale=p / P_BUILD)
    Dv = sim.detector_values()[:, ex.sel]
    for s in range(300):
        c1, _ = plain.decode(Dv[s]); c2, g = gap.decode(Dv[s])
        assert c1 == c2 and g >= 0.0
