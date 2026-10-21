import itertools as it

import numpy as np
import pytest

from primepower_msd.meas_distill import signed_rm
from primepower_msd.meas_sim import Protocol
from primepower_msd.meas_circuit import build, completion, csum_network, _apply_network, forced_verdicts, run
from primepower_msd.zdsim import check_determinism


@pytest.mark.parametrize("d", [4, 8, 16])
def test_encoder_network_and_determinism(d):
    S, tau = signed_rm(1, 4)
    P = Protocol(S, tau, d, np.array(list(it.combinations(range(8), 2))))
    piv, L = completion(np.vstack([P.S, P.N]) % d, d)
    assert np.array_equal(_apply_network(csum_network(L, d), P.n, d), L % d)
    C, _ = build(P, 0, 0, 0)
    assert check_determinism(C) == []
    assert run(P, 0, 0, 0, 100)["accepted"] == 100


def test_circuit_matches_linear_model_713():
    # all weight-3 Z-fault patterns: 85 failing at d = 4 (= meas_sim exact enumeration)
    S, tau = signed_rm(1, 3)
    P = Protocol(S, tau, 4, np.array([[0]]))
    L = P.locations()
    cfgs = [list(zip(c, v)) for c in it.combinations(range(L), 3) for v in it.product(range(1, 4), repeat=3)]
    assert int(forced_verdicts(P, cfgs).sum()) == 85
    cfgs2 = [list(zip(c, v)) for c in it.combinations(range(L), 2) for v in it.product(range(1, 4), repeat=2)]
    assert int(forced_verdicts(P, cfgs2).sum()) == 0
