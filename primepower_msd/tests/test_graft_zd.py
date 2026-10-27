import pytest

from primepower_msd.graft_zd import load_graft, schedule


@pytest.mark.parametrize("d", [4, 8])
def test_codes_commute_and_logicals(d):
    G = load_graft(d)
    for name, code in G["codes"].items():
        X = [c for c in code if c.typ == "X"]; Z = [c for c in code if c.typ == "Z"]
        for a in X:
            wa = a.weights()
            for b in Z:
                wb = b.weights()
                assert sum(wa[q] * wb[q] for q in set(wa) & set(wb)) % d == 0, name
        ZL, XL = G["ZL"], G["XL"]
        for a in X:
            wa = a.weights()
            assert sum(wa[q] * ZL[q] for q in set(wa) & set(ZL)) % d == 0, name
        for b in Z:
            wb = b.weights()
            assert sum(wb[q] * XL[q] for q in set(wb) & set(XL)) % d == 0, name
    s = sum(G["ZL"][q] * G["XL"][q] for q in set(G["ZL"]) & set(G["XL"])) % d
    assert s % 2 == 1


def test_sizes():
    G = load_graft(4)
    codes = G["codes"]
    assert len(set(q for c in codes["colour"] for q in c.support)) == 19
    assert len(set(q for c in codes["grafted"] for q in c.support)) == 229
    assert len(G["fresh"]) == 210 and sum(b == "R" for b in G["fresh"].values()) == 105


def test_schedule_conflict_free():
    G = load_graft(4)
    for code in G["codes"].values():
        layers = schedule(code, G["coords"])
        for L in layers:
            qs = [q for _, q, _ in L]
            assert len(qs) == len(set(qs))
        assert sum(len(L) for L in layers) == sum(len(c.legs) for c in code)
