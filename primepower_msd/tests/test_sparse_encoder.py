import itertools as it

import numpy as np
import pytest

from primepower_msd.meas_distill import _gram, signed_rm
from primepower_msd.meas_sim import Protocol, leading_order
from primepower_msd.meas_circuit import _apply_network, build, completion, csum_network, forced_verdicts, run
from primepower_msd.sparse_encoder import csum_network_ordered, sparse_protocol
from primepower_msd.zdsim import check_determinism


@pytest.mark.parametrize("d", [4, 8])
def test_sparse_normal_basis_and_encoder(d):
    S, tau = signed_rm(1, 4)
    bc = np.array(list(it.combinations(range(8), 2)))
    P, (piv, L, net) = sparse_protocol(S, tau, d, bc, seed=1)
    G, Q = _gram(P.N, tau, d)
    assert not np.any((G - np.diag(np.diag(G))) % d) and np.all(Q % 2)
    assert np.array_equal(_apply_network(net, P.n, d), L % d)
    _, Ld = completion(np.vstack([Protocol(S, tau, d, bc).S, Protocol(S, tau, d, bc).N]) % d, d)
    assert len(net) < len(csum_network(Ld, d))
    C, _ = build(P, 0, 0, 0, (piv, L, net))
    assert check_determinism(C) == []
    assert run(P, 0, 0, 0, 100, encoder=(piv, L, net))["accepted"] == 100


def test_ordered_synthesis_is_exact():
    rng = np.random.default_rng(0)
    for d in (4, 8):
        S, tau = signed_rm(1, 4)
        P = Protocol(S, tau, d, np.array(list(it.combinations(range(8), 2))))
        _, L = completion(np.vstack([P.S, P.N]) % d, d)
        assert np.array_equal(_apply_network(csum_network_ordered(L, d), P.n, d), L % d)


def test_sparse_circuit_matches_linear_713():
    S, tau = signed_rm(1, 3)
    P, enc = sparse_protocol(S, tau, 4, np.array([[0]]), seed=1)
    L = P.locations()
    cfgs = [list(zip(c, v)) for c in it.combinations(range(L), 3) for v in it.product(range(1, 4), repeat=3)]
    assert int(forced_verdicts(P, cfgs, encoder=enc).sum()) == leading_order(P, 3)["C"] == 85
