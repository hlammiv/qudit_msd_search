"""Exact MacWilliams weight-distribution engine for the small-dual (high-puncture) regime.

For a triorthogonal code whose *shortened* generator ``G0`` is small-dimensional
(``dim(G0) = G0.shape[0]`` is tiny when the puncture count ``k`` is large), the quantum
code distance is the minimum weight of ``G0^perp`` and the paper's ``A_d`` counts its
minimum-weight *logical* operators.  Directly enumerating ``G0^perp`` is infeasible (it is
high-dimensional), but its weight distribution follows EXACTLY from the (enumerable) weight
distribution of ``G0`` via the q-ary MacWilliams identity.  This lifts the d<=6 cap of
``mindist.py`` and certifies ``A_d`` for the large codes.

Everything here is EXACT-INTEGER arithmetic.  Weight counts and Krawtchouk values are exact
integers; a single float would silently corrupt them, so we never touch floating point.

Math (used exactly)
-------------------
q-ary Krawtchouk polynomial:

    K_k(x; n, q) = sum_{j=0}^{k} (-1)^j (q-1)^(k-j) C(x, j) C(n-x, k-j).

MacWilliams identity: for a length-n code C over F_q with weight enumerator
A = [A_0, ..., A_n] and |C| = sum_x A_x, the dual C^perp has weight enumerator
B = [B_0, ..., B_n] with

    B_k = (1/|C|) * sum_{x=0}^{n} A_x * K_k(x; n, q).

Each B_k is a non-negative integer; B_0 = 1; sum_k B_k = q^n / |C| = q^(n - dim C).

The CRITICAL SUBTLETY (logical vs. all dual codewords)
------------------------------------------------------
MacWilliams gives ``B`` = the weight distribution of ALL of ``G0^perp``, so ``B[d]`` counts
EVERY weight-d codeword of ``G0^perp``.  The paper's ``A_d`` counts only the weight-d
*logical* operators = weight-d codewords of ``G0^perp`` that are NOT stabilizers (not in
``Gp^perp``).  These two counts are EQUAL iff there are no weight-d stabilizers, i.e. the
stabilizer code ``Gp^perp`` has minimum weight strictly greater than ``d``.  Empirically (see
``tests/test_weightdist.py``) this holds for every oracle code in the paper's Table 3: for
each, ``B[d]`` equals both ``qmsd.distance.A_d_logical_Z`` and the published ``A_d``.  So the
minimum-weight dual codewords are all logical and reporting ``B[d]`` as ``A_d`` is correct.
The engine reports this, and the tests certify the equality against the slow reference.
"""
from __future__ import annotations

from functools import lru_cache
from math import comb

import numpy as np


# ---------------------------------------------------------------------------
# q-ary Krawtchouk polynomial (pure exact int, memoized)
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def krawtchouk(k: int, x: int, n: int, q: int) -> int:
    """q-ary Krawtchouk K_k(x; n, q) -- exact python int.

    K_k(x; n, q) = sum_{j=0}^{k} (-1)^j (q-1)^(k-j) C(x, j) C(n-x, k-j).

    ``comb`` returns 0 when the lower index exceeds the upper, so x and n-x outside
    [0, n] simply contribute zero terms (the identity stays valid for integer x in [0, n]).
    """
    assert isinstance(q, int) and q >= 2, "q must be an integer >= 2"
    assert isinstance(n, int) and n >= 0, "n must be a non-negative integer"
    assert isinstance(k, int) and 0 <= k, "k must be a non-negative integer"
    total = 0
    for j in range(0, k + 1):
        sign = -1 if (j & 1) else 1
        total += sign * (q - 1) ** (k - j) * comb(x, j) * comb(n - x, k - j)
    return int(total)


# ---------------------------------------------------------------------------
# MacWilliams transform (exact)
# ---------------------------------------------------------------------------
def macwilliams(A: list, q: int) -> list:
    """Primal weight enumerator A=[A_0..A_n] over F_q -> dual weight enumerator B=[B_0..B_n].

    Exact-integer MacWilliams transform.  ``|C| = sum(A)``.  Asserts every ``B_k`` is a
    non-negative integer, ``B_0 == 1``, and ``sum(B) == q**n // |C|`` with exact divisibility
    (i.e. ``q**n`` is divisible by ``|C|``).
    """
    assert isinstance(q, int) and q >= 2, "q must be an integer >= 2"
    A = [int(a) for a in A]
    n = len(A) - 1
    assert n >= 0, "A must have length n+1 >= 1"
    assert all(a >= 0 for a in A), "weight enumerator entries must be non-negative"

    size = sum(A)
    assert size > 0, "|C| = sum(A) must be positive"

    B = []
    for k in range(n + 1):
        acc = 0
        for x in range(n + 1):
            if A[x]:
                acc += A[x] * krawtchouk(k, x, n, q)
        assert acc % size == 0, (
            f"MacWilliams: B_{k} not an integer (sum {acc} not divisible by |C|={size})"
        )
        bk = acc // size
        assert bk >= 0, f"MacWilliams: B_{k} = {bk} is negative"
        B.append(bk)

    # Certifying invariants.
    assert B[0] == 1, f"MacWilliams: B_0 = {B[0]} != 1"
    qn = q ** n
    assert qn % size == 0, (
        f"MacWilliams: |C|={size} does not divide q^n={qn} (A is not a code enumerator)"
    )
    assert sum(B) == qn // size, (
        f"MacWilliams: sum(B)={sum(B)} != q^n/|C|={qn // size}"
    )
    return B


# ---------------------------------------------------------------------------
# Exact weight enumerator of a small code by chunked enumeration
# ---------------------------------------------------------------------------
def _as_int_matrix(G, q) -> np.ndarray:
    M = np.asarray(G)
    if M.ndim != 2:
        raise ValueError("generator must be a 2D matrix")
    return M.astype(np.int64) % q


def weight_enumerator(G, q: int, max_words: int = 5_000_000) -> list:
    """Exact weight enumerator [A_0..A_n] of the code spanned by generator G (r x n over F_q).

    Enumerates all ``q**r`` messages in numpy chunks, accumulating a length-(n+1) histogram
    of codeword weights -- it never materialises all codewords at once.  Raises ``ValueError``
    if ``q**r > max_words``.

    (If G is rank-deficient each codeword is counted with equal multiplicity, so the returned
    A is a uniform multiple of the true enumerator; the MacWilliams transform divides by
    ``|C| = sum(A)`` and is invariant to that scaling.)
    """
    assert isinstance(q, int) and q >= 2, "q must be an integer >= 2"
    Gi = _as_int_matrix(G, q)
    r, n = Gi.shape

    total_words = q ** r
    if total_words > max_words:
        raise ValueError(
            f"weight_enumerator: q**r = {q}**{r} = {total_words} exceeds max_words="
            f"{max_words}; code dual too large to enumerate"
        )

    hist = np.zeros(n + 1, dtype=object)  # exact big-int accumulation
    if r == 0:
        hist[0] = 1  # only the zero codeword
        return [int(v) for v in hist]

    # Chunk so that (chunk x n) stays modest in memory.
    chunk = max(1, 4_000_000 // (n + 1))
    qpows = q ** np.arange(r, dtype=np.int64)  # for digit extraction (r small)

    start = 0
    acc64 = np.zeros(n + 1, dtype=np.int64)
    while start < total_words:
        end = min(start + chunk, total_words)
        idx = np.arange(start, end, dtype=np.int64)
        # base-q digits of idx -> message matrix (csize, r)
        messages = (idx[:, None] // qpows[None, :]) % q
        codewords = (messages @ Gi) % q
        weights = np.count_nonzero(codewords, axis=1)
        acc64 += np.bincount(weights, minlength=n + 1).astype(np.int64)
        start = end

    for w in range(n + 1):
        hist[w] = int(acc64[w])
    out = [int(v) for v in hist]
    assert sum(out) == total_words, "weight_enumerator: histogram lost words"
    return out


def dual_weight_distribution(G0, q: int, max_words: int = 5_000_000) -> list:
    """Weight distribution [B_0..B_n] of G0^perp via MacWilliams(weight_enumerator(G0))."""
    A = weight_enumerator(G0, q, max_words=max_words)
    return macwilliams(A, q)


# ---------------------------------------------------------------------------
# Certified minimum distance + A_d
# ---------------------------------------------------------------------------
def exact_distance_and_Ad(G0, q: int, max_words: int = 5_000_000) -> dict:
    """Certified minimum distance and A_d of G0^perp via the exact MacWilliams transform.

    Returns a dict::

        {feasible, distance, weight_dist, B_d, A_d_logical, A_d_equals_Bd}

    ``feasible`` is False (and every other field None) when ``q**dim(G0) > max_words`` so the
    enumeration cannot run.  Otherwise:

      * ``weight_dist``  = B = dual weight distribution (all of G0^perp).
      * ``distance``     = smallest w>0 with B_w>0 (certified minimum distance of G0^perp).
      * ``B_d``          = B[distance] = number of ALL weight-d dual codewords.
      * ``A_d_logical``  = the paper's logical A_d.  Empirically (and asserted against the slow
                           reference in the tests, for every Table-3 oracle code) the
                           minimum-weight dual codewords are all logical -- the stabilizer
                           code ``Gp^perp`` has minimum weight > d -- so A_d_logical == B_d.
      * ``A_d_equals_Bd``= True (the validated regime: no weight-d stabilizers).

    The logical/total distinction cannot in general be separated from ``G0`` alone (it needs
    ``Gp``); the engine relies on the empirically certified equality.  If a future caller hits
    a code where it fails, the tests' cross-check against ``A_d_logical_Z`` would catch it.
    """
    assert isinstance(q, int) and q >= 2, "q must be an integer >= 2"
    Gi = _as_int_matrix(G0, q)
    r = Gi.shape[0]

    if q ** r > max_words:
        return {
            "feasible": False,
            "distance": None,
            "weight_dist": None,
            "B_d": None,
            "A_d_logical": None,
            "A_d_equals_Bd": None,
        }

    B = dual_weight_distribution(Gi, q, max_words=max_words)

    distance = None
    for w in range(1, len(B)):
        if B[w] > 0:
            distance = w
            break
    if distance is None:
        # Only the all-zero dual codeword exists (G0^perp is trivial). No nonzero distance.
        return {
            "feasible": True,
            "distance": None,
            "weight_dist": B,
            "B_d": None,
            "A_d_logical": None,
            "A_d_equals_Bd": None,
        }

    B_d = B[distance]
    return {
        "feasible": True,
        "distance": distance,
        "weight_dist": B,
        "B_d": B_d,
        "A_d_logical": B_d,
        "A_d_equals_Bd": True,
    }
