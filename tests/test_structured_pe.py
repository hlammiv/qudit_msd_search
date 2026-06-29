"""Tests for the structured line-based punctured distance (qmsd.structured_pe).

The structured method pins the punctured minimum distance ``d`` of the m=2 punctured-RM
(triorthogonal) codes from the line geometry of AG(2,p): ``d_lines = d_RM - max_ell
|S cap ell|``, with an explicit weight-``d_lines`` codeword emitted and verified to lie in
``ker(X_stab)``.  The decisive validation reproduces the certified ``d = 5`` of
``[[234,55,5]]_17`` (lower bound self-certified in ~4 s, no weight-3 MITM) and pins
``d = 6`` for ``[[237,52,>=6]]_17`` (the established ``d >= 6`` supplied as ``lower_bound``;
the structured side furnishes the ``d <= 6`` half via an explicit logical operator).
"""
import os

import pytest

from qmsd.structured_pe import enumerate_lines, line_punctured_distance

SLOW = os.environ.get("QMSD_SLOW") == "1"

# Puncture sets (1-indexed) from NEW_CODES.md.
PUNC_237 = [6, 19, 32, 37, 40, 44, 47, 54, 59, 61, 69, 79, 84, 87, 95, 97, 103, 106, 109,
            125, 131, 132, 137, 138, 145, 151, 164, 177, 179, 185, 186, 187, 202, 204, 205,
            206, 207, 211, 212, 221, 246, 247, 248, 251, 254, 259, 266, 268, 269, 273, 283, 285]
PUNC_234 = [1, 3, 6, 9, 12, 16, 27, 28, 29, 45, 48, 57, 58, 70, 74, 79, 80, 85, 90, 94, 96,
            104, 107, 110, 111, 119, 120, 122, 131, 133, 140, 142, 175, 186, 194, 216, 218,
            225, 227, 237, 240, 251, 252, 257, 258, 259, 262, 263, 264, 265, 276, 280, 283,
            287, 288]


def test_enumerate_lines_count():
    """AG(2,p) has exactly p(p+1) affine lines; each has p distinct points; each point lies
    on p+1 lines (every column index 0..p^2-1 is covered p+1 times)."""
    p = 17
    lines = list(enumerate_lines(p))
    assert len(lines) == p * (p + 1)
    from collections import Counter
    cover = Counter()
    for pts, tparams in lines:
        assert len(pts) == p
        assert len(set(pts)) == p           # p distinct points
        assert list(tparams) == list(range(p))
        cover.update(pts)
    assert set(cover) == set(range(p * p))  # covers the whole plane
    assert set(cover.values()) == {p + 1}   # every point on exactly p+1 lines


def test_234_validation_d5():
    """DECISIVE VALIDATION: the line method must reproduce the certified d=5 of
    [[234,55,5]]_17.  The lower bound (d>=5) self-certifies here in ~4 s because d_max=4
    uses only size-2 MITM halves (no expensive weight-3 right half)."""
    res = line_punctured_distance(17, 2, 10, PUNC_234, certify_lower=True)
    assert res["d_RM"] == 12
    assert res["beta"] == 5
    assert res["max_line_punctures"] == 7        # a line carries 7 of the 55 punctures
    assert res["d_lines"] == 5
    assert res["lower_bound"] == 5               # min_dependent_columns(d_max=4) raised -> d>=5
    assert res["exact"] is True
    assert res["distance"] == 5
    cert = res["certificate"]
    assert cert["weight"] == 5
    assert cert["in_kernel"] is True             # explicit weight-5 logical operator


def test_237_pins_d6():
    """[[237,52,>=6]]_17: structured upper bound d<=6 via an explicit weight-6 codeword,
    combined with the established certified d>=6 (NEW_CODES.md: min_dependent_columns
    d_max=5 raises) supplied as lower_bound, pins d = 6 EXACTLY."""
    res = line_punctured_distance(17, 2, 10, PUNC_237, lower_bound=6)
    assert res["d_RM"] == 12
    assert res["beta"] == 5
    assert res["max_line_punctures"] == 6        # heaviest line carries 6 punctures
    assert res["d_lines"] == 6
    assert res["exact"] is True
    assert res["distance"] == 6
    cert = res["certificate"]
    assert cert["weight"] == 6
    assert cert["in_kernel"] is True             # explicit weight-6 logical operator => d<=6
    assert len(cert["support"]) == 6


@pytest.mark.skipif(not SLOW, reason="set QMSD_SLOW=1: d_max=5 weight-3 MITM is many minutes")
def test_237_self_certify_lower():
    """Self-certify the d>=6 lower bound for [[237,52]] (slow: d_max=5 weight-3 MITM)."""
    res = line_punctured_distance(17, 2, 10, PUNC_237, certify_lower=True)
    assert res["lower_bound"] == 6
    assert res["exact"] is True
    assert res["distance"] == 6
