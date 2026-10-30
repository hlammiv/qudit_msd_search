"""Low-weight words of RM_7(20,4) = RM(3,4)^perp as point sets (WEIGHTED_SUPPORTS.md, "Weight-7
words of RM(20,4)"): the batched F_p rank, the line words, the weight-5/6/7 line lemmas on
samples and hand-built configurations, the point-set classifier and the <= 2-column direction
table used by the mixed stages of the p = 7 certificates.  Fast (a few seconds)."""
import itertools

import numpy as np
import pytest

from qmsd.pointsets import (affine_lines, batched_rank_modp, classify_pointset, collinearity_profile,
                            direction_keys, frame_and_flat, full_support_null_vector, lookup_directions,
                            pair_direction_table)
from qmsd.reedmuller import points, rm_generator
from qmsd.structured_ad import _null_space, _rank

p, m = 7, 4
P = np.asarray(points(m, p), dtype=np.int64)
POW = p ** np.arange(m)
R3 = np.asarray(rm_generator(3, m, p), dtype=np.int64) % p
IDX = {tuple(int(v) for v in x): i for i, x in enumerate(P)}


def idx(*coords):
    return [IDX[tuple(c)] for c in coords]


def test_batched_rank_matches_serial_rank():
    rng = np.random.default_rng(3)
    A = rng.integers(0, p, (600, 7, 9))
    A[:100, 6] = A[:100, 5]                      # repeated row
    A[100:200, :, 3] = 0                         # zero column
    A[200:300] = (2 * A[200:300][:, [0, 0, 1, 2, 3, 4, 5]]) % p
    A[300:350, 2:] = 0                           # rank <= 2
    assert (batched_rank_modp(A, p) == np.array([_rank(a, p) for a in A])).all()
    assert (batched_rank_modp(np.swapaxes(A, 1, 2), p) == np.array([_rank(a.T, p) for a in A])).all()
    for q in (3, 5):
        B = rng.integers(0, q, (200, 5, 6)); B[:50, 4] = (B[:50, 0] + B[:50, 1]) % q
        assert (batched_rank_modp(B, q) == np.array([_rank(b, q) for b in B])).all()


def test_affine_lines_are_all_lines_in_parameter_order():
    L = affine_lines(p, m)
    assert L.shape == (137200, 7) and len({tuple(sorted(r)) for r in L}) == 137200
    d = (P[L[:, 1]] - P[L[:, 0]]) % p
    for t in range(2, p):
        assert (((P[L[:, 0]] + t * d) % p) @ POW == L[:, t]).all()
    assert affine_lines(5, 2).shape == (30, 5) and affine_lines(3, 3).shape == (117, 3)


def test_line_words_of_rm20_4():
    """RM(3,4)^perp on a line is RM(2,1): 57 projective words, 21 of weight 5, 14 of weight 6
    (linear and squares), 22 of weight 7 (constants and irreducible quadratics)."""
    line = affine_lines(p, m)[12345]
    ns = np.asarray(_null_space(R3[:, line], p)) % p
    assert ns.shape == (3, 7)
    seen = {}
    for coef in itertools.product(range(p), repeat=3):
        if any(coef):
            u = (np.array(coef) @ ns) % p
            u = tuple(((u * pow(int(u[np.nonzero(u)[0][0]]), -1, p)) % p).tolist())
            seen[u] = sum(1 for v in u if v)
    assert len(seen) == 57 and sorted(np.unique(list(seen.values()), return_counts=True)[1]) == [14, 21, 22]
    assert dict(zip(*np.unique(list(seen.values()), return_counts=True))) == {5: 21, 6: 14, 7: 22}
    u = full_support_null_vector(R3[:, line], p)
    assert u is not None and u.all() and not ((R3[:, line] @ u) % p).any()


def test_weight_5_6_line_lemma_on_a_sample():
    """Every dependent 5- or 6-point set of a random plane / 3-flat is collinear (weight-5/6 words
    of RM(20,4) live on lines), and sets with 5 collinear points + off-line points are dependent
    but never full-support."""
    rng = np.random.default_rng(11)
    frame2, plane = frame_and_flat(m, 2, p); frame3, flat3 = frame_and_flat(m, 3, p)
    for pool, n_pts in ((plane, 5), (plane, 6), (flat3, 5), (flat3, 6)):
        T = np.array([rng.choice(pool, n_pts, replace=False) for _ in range(3000)])
        rk = batched_rank_modp(R3.T[T], p)
        n_dep = 0
        for Tset in T[rk < n_pts]:
            n_dep += 1
            mc = collinearity_profile(P[Tset], p)["max_collinear"]
            assert mc >= 5                                            # dependent => 5 collinear (d = 5)
            assert (full_support_null_vector(R3[:, Tset], p) is not None) == (mc == n_pts)
        assert n_dep > 0 if pool is plane else True
    L = affine_lines(p, m)
    for _ in range(50):
        line = L[rng.integers(len(L))]
        off = [x for x in rng.choice(p ** m, 3, replace=False) if x not in set(line)][:1]
        T = list(rng.choice(line, 5, replace=False)) + off
        assert _rank(R3[:, T], p) == 5 and full_support_null_vector(R3[:, T], p) is None


def test_weight7_lemma_on_hand_built_configurations():
    """Only a full line carries a full-support weight-7 annihilator of the cubics; the natural
    competitors (6 + 1, 5 + 2, two lines 4 + 3 / 3 + 3 + 1, parabola, 7 of an ellipse, a 3-flat
    frame + 3, a span-4 frame + 2) do not."""
    line = idx(*[(t, 0, 0, 0) for t in range(7)])
    assert classify_pointset(P[line], p) == "line7" and full_support_null_vector(R3[:, line], p) is not None
    cases = {
        "line6+1": line[:6] + idx((0, 1, 0, 0)),
        "line5+2": line[:5] + idx((0, 1, 0, 0), (3, 2, 1, 0)),
        "span2:lines[4, 3]:quadric": idx(*[(t, 0, 0, 0) for t in range(4)], *[(t, 1, 0, 0) for t in range(3)]),
        "span2:lines[4, 4]:quadric": idx(*[(t, 0, 0, 0) for t in range(4)], *[(0, t, 0, 0) for t in range(1, 4)]),
        "span2:lines[3, 3]": idx(*[(t, 0, 0, 0) for t in range(3)], *[(t, 1, 0, 0) for t in range(3)], (3, 4, 0, 0)),
        "span2:lines[3, 3, 3]": idx(*[(t, 0, 0, 0) for t in range(3)], *[(t, 1, 0, 0) for t in range(3)], (0, 3, 0, 0)),
        "span2:lines[]:quadric": idx(*[(t, (t * t) % p, 0, 0) for t in range(7)]),                 # parabola
        "span3:lines[]": idx((0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (1, 1, 1, 0), (2, 3, 1, 0), (5, 1, 4, 0)),
        "span4:lines[]": idx((0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 2, 3, 4), (3, 3, 1, 5)),
    }
    for label, T in cases.items():
        assert classify_pointset(P[T], p) == label, (label, classify_pointset(P[T], p))
        assert full_support_null_vector(R3[:, T], p) is None, label
    ell = [IDX[(x, y, 0, 0)] for x in range(p) for y in range(p) if (x * x + y * y - 1) % p == 0]
    assert len(ell) == 8 and _rank(R3[:, ell], p) == 7 and full_support_null_vector(R3[:, ell], p) is not None
    for T in itertools.combinations(ell, 7):                          # 7 of the 8 conic points: independent
        assert _rank(R3[:, list(T)], p) == 7
    assert collinearity_profile(P[cases["span3:lines[]"]], p)["max_coplanar"] == 4


def test_weight7_lemma_on_frame_sets_sample():
    """Frame enumeration (as in scratchpad/rm20_w7_classify.py) on random samples: every
    dependent 7-set of span 2, 3, 4 contains >= 5 collinear points and none has a full-support word."""
    rng = np.random.default_rng(5)
    n_dep = 0
    for j, n_sample in ((2, 20000), (3, 20000), (4, 20000)):
        frame, flat = frame_and_flat(m, j, p)
        others = np.array([x for x in flat if x not in set(frame)])
        E = np.array([rng.choice(others, 7 - (j + 1), replace=False) for _ in range(n_sample)])
        T = np.concatenate([np.broadcast_to(np.array(frame), (n_sample, j + 1)), E], axis=1)
        rk = batched_rank_modp(R3.T[T], p)
        for Tset in T[rk < 7]:
            n_dep += 1
            lab = classify_pointset(P[Tset], p)
            assert lab in ("line5+2", "line6+1"), lab
            assert full_support_null_vector(R3[:, Tset], p) is None
    assert n_dep > 0                                                  # the span-2 sample always has 5-collinear sets


def test_pair_direction_table_lookup():
    rng = np.random.default_rng(9)
    C = rng.integers(0, p, (12, 25))
    tab = pair_direction_table(C, p)
    assert len(tab["keys"]) == 1 + 25 + 6 * 300
    S = np.stack([np.zeros(12, dtype=np.int64), 3 * C[:, 4] % p, (2 * C[:, 7] + 5 * C[:, 19]) % p,
                  (C[:, 0] + C[:, 1] + C[:, 2]) % p, rng.integers(0, p, 12)])
    pos = lookup_directions(tab, S, p)
    assert pos[0] >= 0 and tab["kind"][pos[0]] == 0
    assert tab["kind"][pos[1]] == 1 and tab["i"][pos[1]] == 4
    assert tab["kind"][pos[2]] == 2 and (tab["i"][pos[2]], tab["j"][pos[2]]) == (7, 19)
    assert pos[3] == -1 and pos[4] == -1
    beta = tab["beta"][pos[2]]
    assert (direction_keys((C[:, 7] + beta * C[:, 19])[None, :] % p, p) == direction_keys(S[2:3], p)).all()
