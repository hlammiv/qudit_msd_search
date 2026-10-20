"""Small exhaustive checks of the flat lemma behind the ring-punctured RM distance certificate."""

import itertools as it

import numpy as np

from primepower_msd import rm_puncture_cert as C


def _forms_vanishing_on(G, m):
    fb, span = [], {0}
    for f in range(1, 2 ** (m + 1)):
        if f == (1 << m) or f in span:
            continue
        if all(C.form_values(np.array([g]), f, m)[0] == 0 for g in G):
            fb.append(f)
            span |= {s ^ f for s in span}
    return fb


def test_flat_lemma_m7_r2():
    """Every logical S = A ∪ R (|R| <= 3) of punctured RM(2,7) has A \\ aff(R) in X_{aff R};
    FlatBFS.check agrees with exact_flat_check."""
    m, r = 7, 2
    rng = np.random.default_rng(3)
    allpts = np.array([[(x >> i) & 1 for i in range(m)] for x in range(2 ** m)], dtype=np.uint8)
    Phi = C.phi_matrix(allpts, r)
    enc = [int("".join(map(str, row)), 2) for row in Phi]
    checked = 0
    for _ in range(3):
        while True:
            Pi = sorted(rng.choice(2 ** m, size=10, replace=False).tolist())
            if C.gf2_rank_rows(Phi[Pi]) == len(Pi):
                break
        P = allpts[Pi]
        bfs = C.FlatBFS(P, r)
        basis = {}

        def reduce(v, tr=0):
            while v:
                h = v.bit_length() - 1
                if h not in basis:
                    return v, tr
                bv, bt = basis[h]
                v ^= bv
                tr ^= bt
            return 0, tr
        for j, p in enumerate(Pi):
            v, tr = reduce(enc[p], 1 << j)
            basis[v.bit_length() - 1] = (v, tr)
        comp = [x for x in range(2 ** m) if x not in set(Pi)]
        red = {x: reduce(enc[x]) for x in comp}
        for size in (1, 2, 3):
            for R in it.combinations(comp, size):
                v = tr = 0
                for x in R:
                    v ^= red[x][0]
                    tr ^= red[x][1]
                if v or not tr:
                    continue
                G = C.affine_span(list(R))
                fb = _forms_vanishing_on(G, m)
                ex = C.exact_flat_check(P, fb, r)
                assert bfs.check(fb)[0] == (ex["dim_X"] > 0)
                outside = np.array([(tr >> j & 1) and (Pi[j] not in set(G)) for j in range(len(Pi))], dtype=np.uint8)
                if outside.any():
                    Bm = np.array(ex["basis"], dtype=np.uint8).reshape(-1, len(Pi))
                    assert len(Bm) and C.gf2_rank_rows(np.vstack([Bm, outside])) == C.gf2_rank_rows(Bm)
                checked += 1
    assert checked > 0


def test_sidon_and_witness_helpers():
    pts = [0b0001, 0b0010, 0b0100, 0b1000, 0b0000, 0b1111]
    assert C.is_sidon(pts[:5])
    assert not C.is_sidon([0, 1, 2, 3])
    w = C.flat_witness([0, 1, 2], [0, 1, 2])
    assert w["size"] == 4 and w["occupancy"] == 3 and w["R"] == [3]
