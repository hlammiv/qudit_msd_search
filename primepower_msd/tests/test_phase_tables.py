"""Phase-table hierarchy oracle, decomposition and enumeration (validated against dense oracles)."""

import itertools as it

import numpy as np

from primepower_msd import phase_tables as PT
from primepower_msd.clifford_ring import level_of


def test_single_qudit_level_matches_dense_oracle_all_tables():
    bad = 0
    for rest in it.product(range(16), repeat=3):
        P = np.array((0,) + rest)
        lv = PT.level(P)
        if lv == 0:
            continue
        if lv != level_of(PT.to_unitary(P), 4):
            bad += 1
    assert bad == 0


def test_two_qudit_level_matches_dense_hs_oracle():
    rng = np.random.default_rng(2)
    samples = PT.sample_level3(2, 30, rng) + [PT.normalize(rng.integers(0, 16, size=(4, 4))) for _ in range(15)]
    for P in samples:
        assert PT.level(P) == PT.level_of_dense(PT.to_unitary(P), 2)


def test_carriers_and_strict_level3_count():
    assert len(PT.CARRIERS_D4) == 12
    assert PT.CARRIER_Q in PT.CARRIERS_D4
    s3 = set(PT.strict_level3_single(4))
    assert len(s3) == 48                       # (level<=3 with P(0)=0: 64) - (Clifford: 16)
    assert all(tuple(c) in s3 for c in PT.CARRIERS_D4)
    assert PT.level3_module_size(1) == 1024


def test_decomposition_roundtrip_and_parts():
    rng = np.random.default_rng(3)
    for m in (1, 2, 3):
        P = PT.normalize(rng.integers(0, 16, size=(4,) * m))
        dec = PT.interaction_decomposition(P)
        tot = sum(dec.values()) % 16
        assert np.array_equal(tot, P)
    ccz = PT.ccz_layer(3, [(0, 1, 2, 1)])
    c = PT.classify(ccz)
    assert c["level"] == 3 and c["entangling"] and not c["escapes_A2"]
    assert c["ccz_terms"] == [{"qudits": (0, 1, 2), "coeff": 1, "pure_trilinear": True}]
    assert c["odd_third_difference"] and not c["genuine_T"]
    T = PT.from_qudit_tables([PT.CARRIER_Q, None])
    c = PT.classify(T)
    assert c["genuine_T"] and not c["entangling"] and c["one_body"][0]["escapes"]
    assert PT.level(PT.cz_layer(2, [(0, 1, 1)])) == 2


def test_batch_flags_agree_with_classify():
    rng = np.random.default_rng(4)
    S = PT.sample_level3(2, 80, rng)
    E = np.array([P.reshape(-1) for P in S])
    F = PT.batch_flags(E, 2)
    for i, P in enumerate(S):
        c = PT.classify(P)
        assert (c["level"] <= 3) == F["level3"][i]
        assert (c["level"] <= 2) == F["clifford"][i]
        assert c["entangling"] == F["entangling"][i]
        assert c["escapes_all_sectors"] == F["escapes"][i]
        assert c["one_body_escapes"] == F["one_body_escapes"][i]
        assert c["mixed_part_escapes"] == F["mixed_escapes"][i]


def test_m1_catalog_numbers():
    cat = PT.catalog(1)
    assert cat["counts"]["level3_strict"] == 48
    assert cat["counts"]["escaping"] == 32
    assert cat["counts"]["entangling"] == 0


def test_m2_catalog_has_entangling_escapers():
    cat = PT.catalog(2)
    assert cat["mode"] == "exhaustive"
    assert cat["counts"]["total"] == 131072
    assert cat["counts"]["lemma_violations"] == 65536   # escape carried by the entangling part
