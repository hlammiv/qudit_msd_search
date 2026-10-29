import numpy as np
import pytest

from primepower_msd.logical_csum import CsumExperiment, csum_circuit
from primepower_msd.surface_zd import P_BUILD
from primepower_msd.zdsim import FrameSim, check_determinism


@pytest.mark.parametrize("d", [4, 8])
@pytest.mark.parametrize("basis", ["Z", "X"])
def test_detectors_deterministic(d, basis):
    mem = csum_circuit(3, d, 2, 2, P_BUILD, basis)
    assert check_determinism(mem.circ, 64, 0) == []


def test_logical_action_z_basis():
    """X on row 0 of the control (logical X_c) just before the layer flips both logical Z readouts."""
    D, d = 3, 4
    mem = csum_circuit(D, d, 2, 2, P_BUILD, "Z")
    row0 = [mem.data[0][j] for j in range(D)]
    mask = np.ones(1, dtype=np.int64)
    P = np.zeros((1, D, 2), dtype=np.int64); P[:, :, 0] = 1
    sim = FrameSim(mem.circ, 1, np.random.default_rng(0), gauge=False).run(scale=0.0, forced={mem.csum_op - 1: [(mask, row0, P)]})
    assert not sim.detector_values().any()
    vals = [sum(cf * sim.meas[0, i] for i, cf in o) % d for o in mem.obs_all]
    assert vals == [1, 1]


def test_experiment_runs_and_decodes():
    ex = CsumExperiment(3, 4, "Z")
    assert all(dem.undetectable_prob == 0 for dem in ex.dems)
    r = ex.run(1e-2, 400, seed=1)
    assert 0 < r["fails"] < 400
    r0 = CsumExperiment(3, 4, "Z", csum=False).run(1e-2, 400, seed=1)
    assert 0 < r0["fails"] < 400


@pytest.mark.parametrize("basis", ["Z", "X"])
def test_sequential_decoder_single_faults(basis):
    """Every single basis fault of the D=3 CSUM circuit is corrected by the Pauli-frame decoder (both observables)."""
    from primepower_msd.surface_zd import propagate_basis
    ex = CsumExperiment(3, 4, basis, decoder="seq")
    dec = ex.make_decoder(1e-3)
    assert ex.seq.n_shadow > 0
    faults, per_loc, V, O = propagate_basis(ex.mem)
    Vs = V[:, ex.sel]
    Oall = []
    for o in ex.mem.obs_all:
        _, _, _, Oo = propagate_basis(__import__("dataclasses").replace(ex.mem, obs=o))
        Oall.append(Oo)
    wrong = 0
    for f in range(len(faults)):
        corr = dec(Vs[f].copy())
        wrong += any((Oall[k][f] + corr[k]) % 4 for k in range(2))
    assert wrong == 0, f"{wrong} of {len(faults)} single faults mis-decoded"


def test_sequential_vs_independent():
    ex = CsumExperiment(3, 4, "Z", decoder="seq")
    r = ex.run(1e-2, 400, seed=1)
    assert 0 < r["fails"] < 400
