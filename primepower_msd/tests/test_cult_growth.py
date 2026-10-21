import numpy as np
import pytest

from primepower_msd.cult_growth import GrowCult, embedding, logical_reps
from primepower_msd.zdsim import check_determinism


def test_embedding_corner():
    Hs, Hb, xb, emb = embedding(3, 5)
    big = {frozenset(np.nonzero(r)[0]) for r in Hb}
    for r in Hs:                                   # every D=3 plaquette is (contained in) a D=5 plaquette
        s = frozenset(emb[i] for i in np.nonzero(r)[0])
        assert any(s <= b for b in big)


@pytest.mark.parametrize("d", [2, 4, 8])
def test_growth_deterministic(d):
    g = GrowCult(d=d, p=0.0)
    assert check_determinism(g.c, 128) == []
    assert any(t == "T-check" for t in g.c.det_tags[-len(g.c.det_tags) // 2:])


def test_growth_controls():
    g = GrowCult(d=4, p=0.0)
    g.c.ops = [op for op in g.c.ops if op[0] != "FF"]          # no gauge fixing -> random detectors
    assert check_determinism(g.c, 128)
    g = GrowCult(d=4, p=0.0, sigma=1, inj=3)                   # wrong injected state is caught
    assert check_determinism(g.c, 128)


def test_growth_infeasible_split_rejected():
    Hs, Hb, xb, emb = embedding(3, 5)
    new = [i for i in range(Hb.shape[1]) if i not in emb]
    g = GrowCult(d=4, p=0.0)
    # all new qudits in |0> cannot carry Xbar = X^{all} extension -> some split must be infeasible
    assert logical_reps(Hb, g.MZb, Hs, g.MZs, emb, new, 4) is None


def test_growth_single_faults_d2():
    r = GrowCult(d=2, p=1e-3).single_fault_scan()
    assert r["undetected_harmful"] == 0
