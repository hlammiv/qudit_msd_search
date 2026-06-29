"""Correctness core for the distributed qutrit weight-distribution engine.

This module owns the EXACT, certified parts of the pipeline:

  1.  int64 overflow proof for the histogram accumulators (see ``assert_int64_safe``).
  2.  the one-time arbitrary-precision q-ary MacWilliams / Krawtchouk transform that maps
      the (int64) primal weight enumerator ``A`` of ``G0`` to the dual enumerator
      ``B`` of ``G0^perp`` -- ``B`` entries grow to ~3^(n-K) ~ 3^1913, hence Python big ints.
  3.  extraction of the certified minimum distance ``d`` and ``A_d`` from ``B``.
  4.  the certifying invariant checks (B integer, B_0 = 1, sum_w B_w = q^(n-K), all B_w >= 0).
  5.  helpers to merge per-block int64 partial histograms with a per-block checksum.

Everything is exact integer arithmetic.  A single float would silently corrupt a bignum,
so floating point is never touched.

------------------------------------------------------------------------------------------
WHY int64 CANNOT OVERFLOW (proof)
------------------------------------------------------------------------------------------
Let G0 be an [n, K] code over F_q, q = 3, with K = #rows enumerated.  Enumeration runs over
all q^K messages m and histograms wt(m @ G0).  For any bin w,

        A_w = #{ m in F_q^K : wt(m @ G0) = w }  <=  q^K          (every bin <= #messages),

and the TOTAL over all bins is exactly sum_w A_w = q^K.  A partial histogram over a subset
of the message blocks has every bin and its total bounded by q^K as well.  The global
enumerator is the exact vector sum of the partials and again has total q^K.  Hence every
int64 quantity that ever appears -- a partial bin, the merged bin, the running total -- is
bounded by

        q^K = 3^K  <=  3^29 = 68,630,377,364,883  ~ 6.86e13.

int64 holds up to 2^63 - 1 = 9,223,372,036,854,775,807 ~ 9.22e18, a margin of ~1.3e5x over
the largest representable count at the K<=29 feasibility ceiling.  No int64 add in the hot
enumeration or the partial merge can overflow.  (Even at K = 39, 3^39 ~ 4.05e18 < 2^63, so
int64 is safe well beyond the compute-feasible range; we still assert K<=29 as the design
contract.)  The MacWilliams transform that follows uses big ints, NOT int64.
"""
from __future__ import annotations

from math import comb
from typing import Iterable, Sequence

# ---------------------------------------------------------------------------
# int64 overflow proof, enforced.
# ---------------------------------------------------------------------------
INT64_MAX = 2**63 - 1
DESIGN_K_CEILING = 29  # feasibility ceiling on ~52 threads (3^29 codewords).


def assert_int64_safe(K: int, q: int = 3) -> None:
    """Assert that every histogram count for an [n, K] code over F_q fits in a signed int64.

    Proof obligation: the largest count is q^K (a single bin can hold every message, and the
    histogram total is exactly q^K).  We require q^K <= INT64_MAX.  For the design contract we
    also require K <= DESIGN_K_CEILING so the *compute* (not just the arithmetic) is sane.
    """
    assert isinstance(K, int) and K >= 0, "K must be a non-negative int"
    assert isinstance(q, int) and q >= 2, "q must be an int >= 2"
    assert K <= DESIGN_K_CEILING, (
        f"K={K} exceeds the design ceiling {DESIGN_K_CEILING} (3^29 ~ 6.9e13 codewords)"
    )
    assert q**K <= INT64_MAX, (
        f"q^K = {q}^{K} = {q**K} exceeds INT64_MAX={INT64_MAX}; int64 histogram would overflow"
    )


# ---------------------------------------------------------------------------
# q-ary Krawtchouk: exact big-int, via the three-term recurrence in the degree.
# ---------------------------------------------------------------------------
def krawtchouk_direct(k: int, x: int, n: int, q: int) -> int:
    """Reference q-ary Krawtchouk K_k(x; n, q) by its defining sum (exact int).

        K_k(x; n, q) = sum_{j=0}^{k} (-1)^j (q-1)^(k-j) C(x, j) C(n-x, k-j).

    Identical in value to ``qmsd.weightdist.krawtchouk``; used to validate the recurrence.
    """
    total = 0
    for j in range(k + 1):
        sign = -1 if (j & 1) else 1
        total += sign * (q - 1) ** (k - j) * comb(x, j) * comb(n - x, k - j)
    return total


def krawtchouk_column(x: int, n: int, q: int, kmax: int) -> list:
    """Return [K_0(x), K_1(x), ..., K_kmax(x)] for fixed evaluation point x, exact big ints.

    Uses the standard q-ary Krawtchouk three-term recurrence in the degree k:

        (k+1) K_{k+1}(x) = [ (q-1)(n-k) + k - q*x ] K_k(x) - (q-1)(n-k+1) K_{k-1}(x),

    with K_0(x) = 1 and K_1(x) = (q-1)*n - q*x.  Each step is O(1) big-int ops, so the column
    costs O(kmax).  Building the full table for the transform is O(n * kmax) big-int ops.
    """
    assert 0 <= x <= n and kmax >= 0
    K = [0] * (kmax + 1)
    K[0] = 1
    if kmax >= 1:
        K[1] = (q - 1) * n - q * x
    for k in range(1, kmax):
        num = ((q - 1) * (n - k) + k - q * x) * K[k] - (q - 1) * (n - k + 1) * K[k - 1]
        assert num % (k + 1) == 0, "Krawtchouk recurrence produced a non-integer (bug)"
        K[k + 1] = num // (k + 1)
    return K


# ---------------------------------------------------------------------------
# MacWilliams transform (exact big int) + certifying invariant checks.
# ---------------------------------------------------------------------------
def macwilliams(A: Sequence[int], q: int, kmax: int | None = None) -> list:
    """Primal enumerator A=[A_0..A_n] over F_q -> dual enumerator B=[B_0..B_kmax], exact big int.

        B_k = (1/|C|) * sum_x A_x * K_k(x; n, q),   |C| = sum_x A_x.

    Computed by transposing the loop over the recurrence: accumulate, for each evaluation
    point x with A_x != 0, the column [K_0(x)..K_kmax(x)] scaled by A_x.  This evaluates every
    Krawtchouk by recurrence (no per-(k,x) re-summation) and is the efficient one-time path.

    When ``kmax is None`` the full B (length n+1) is returned and the full set of certifying
    invariants is asserted (B_0=1, sum_k B_k = q^n/|C|, every B_k a non-negative integer, and
    |C| divides q^n).  When ``kmax`` is given only B_0..B_kmax are returned (the fast path for
    extracting (d, A_d)); the per-entry integer/non-negativity/B_0 checks still run, but the
    global sum invariant is skipped (it needs the full vector).
    """
    assert isinstance(q, int) and q >= 2, "q must be an int >= 2"
    A = [int(a) for a in A]
    n = len(A) - 1
    assert n >= 0, "A must have length n+1 >= 1"
    assert all(a >= 0 for a in A), "weight enumerator entries must be non-negative"
    size = sum(A)
    assert size > 0, "|C| = sum(A) must be positive"

    full = kmax is None
    kk = n if full else int(kmax)
    assert 0 <= kk <= n, f"kmax={kmax} out of range [0, {n}]"

    acc = [0] * (kk + 1)  # acc[k] = sum_x A_x * K_k(x)   (big ints)
    for x in range(n + 1):
        ax = A[x]
        if not ax:
            continue
        col = krawtchouk_column(x, n, q, kk)
        for k in range(kk + 1):
            acc[k] += ax * col[k]

    B = []
    for k in range(kk + 1):
        assert acc[k] % size == 0, (
            f"MacWilliams: B_{k} not an integer (sum {acc[k]} not divisible by |C|={size})"
        )
        bk = acc[k] // size
        assert bk >= 0, f"MacWilliams: B_{k} = {bk} is negative"
        B.append(bk)

    assert B[0] == 1, f"MacWilliams: B_0 = {B[0]} != 1"
    if full:
        qn = q**n
        assert qn % size == 0, (
            f"MacWilliams: |C|={size} does not divide q^n={qn} (A is not a code enumerator)"
        )
        assert sum(B) == qn // size, (
            f"MacWilliams: sum(B)={sum(B)} != q^n/|C|={qn // size}"
        )
    return B


def assert_dual_invariants(B: Sequence[int], A: Sequence[int], q: int) -> None:
    """Independently re-check the certifying invariants of a full dual enumerator B.

    B_0 = 1; every B_w a non-negative integer; |C| = sum(A) divides q^n; sum_w B_w = q^n/|C|.
    Use after a full ``macwilliams`` (or after deserializing B) as a standalone audit.
    """
    A = [int(a) for a in A]
    B = [int(b) for b in B]
    n = len(A) - 1
    assert len(B) == n + 1, f"B has length {len(B)}, expected n+1={n + 1}"
    size = sum(A)
    assert size > 0, "|C| must be positive"
    assert all(b >= 0 for b in B), "some B_w is negative"
    assert B[0] == 1, f"B_0 = {B[0]} != 1"
    qn = q**n
    assert qn % size == 0, f"|C|={size} does not divide q^n"
    assert sum(B) == qn // size, f"sum(B)={sum(B)} != q^n/|C|={qn // size}"


# ---------------------------------------------------------------------------
# (d, A_d) extraction.
# ---------------------------------------------------------------------------
def extract_d_and_Ad(A: Sequence[int], q: int, search_kmax: int | None = None) -> dict:
    """Certified minimum distance d of G0^perp and A_d, from the primal enumerator A of G0.

    A_d = B_d where B = MacWilliams(A); d = min{ w > 0 : B_w > 0 }.  By default this scans the
    FULL B (and asserts all global invariants), which also certifies B_1..B_{d-1} = 0.  When
    ``search_kmax`` is given (cheaper for the n~1939 target) only B_0..B_{search_kmax} are
    computed; d must be found within that window or a ValueError is raised so a too-small
    window can never be mistaken for "no low-weight codeword".

    Returns {n, K, q, d, A_d, B (or None), B_window}.
    """
    A = [int(a) for a in A]
    n = len(A) - 1
    size = sum(A)
    # K recovered from |C| = q^K (full-rank G0).  If rank-deficient |C| is still a power of q
    # times a uniform multiplicity; the MacWilliams result is scale-invariant either way.
    K = None
    if size > 0:
        kk = 0
        s = 1
        while s < size:
            s *= q
            kk += 1
        if s == size:
            K = kk

    if search_kmax is None:
        B = macwilliams(A, q)  # full, all invariants checked
        window = B
        full = True
    else:
        B = None
        window = macwilliams(A, q, kmax=search_kmax)
        full = False

    d = None
    for w in range(1, len(window)):
        if window[w] > 0:
            d = w
            break
    if d is None and not full:
        raise ValueError(
            f"extract_d_and_Ad: no nonzero B_w found within window kmax={search_kmax}; "
            f"increase search_kmax (true d may be larger)"
        )
    A_d = window[d] if d is not None else None
    return {"n": n, "K": K, "q": q, "d": d, "A_d": A_d, "B": B, "B_window": window}


# ---------------------------------------------------------------------------
# Partial-histogram merge with per-block checksum.
# ---------------------------------------------------------------------------
def merge_partials(partials: Iterable[Sequence[int]], n: int, expected_total: int | None = None) -> list:
    """Exact int64-safe vector sum of per-block partial histograms (returned as a Python list).

    Each partial is a length-(n+1) non-negative integer vector whose entries are bounded by
    q^K (proved overflow-safe; see ``assert_int64_safe``).  The CHECKSUM contract is that the
    sum of all entries across all merged partials equals ``expected_total`` (the number of
    messages covered = q^K for a complete run).  Raises AssertionError on mismatch.
    """
    acc = [0] * (n + 1)
    running_total = 0
    for p in partials:
        assert len(p) == n + 1, f"partial has length {len(p)}, expected n+1={n + 1}"
        for w in range(n + 1):
            v = int(p[w])
            assert v >= 0, "partial histogram has a negative entry"
            acc[w] += v
            running_total += v
    if expected_total is not None:
        assert running_total == expected_total, (
            f"merge checksum failed: total messages {running_total} != expected {expected_total}"
        )
    return acc


def block_checksum_ok(partial: Sequence[int], n_messages_in_block: int) -> bool:
    """A single block's partial histogram must sum to the number of messages it covered."""
    return sum(int(v) for v in partial) == int(n_messages_in_block)
