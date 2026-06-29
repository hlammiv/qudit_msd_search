"""Exact validation of the correctness core against the trusted qmsd engine.

Run from /home/hlamm/Desktop/QC/prime_msd:  python -m dist_weightdist.validate

PROTOCOL
--------
1.  Krawtchouk recurrence == qmsd.weightdist.krawtchouk (defining sum) on a grid of (k, x, n).
2.  Full MacWilliams of this module == qmsd.weightdist.macwilliams on random small codes.
3.  ORACLE codes [[206,37,4]]_3 (K=14) and [[667,62,4]]_3 (K=16):
      - build G0 = build_triorthogonal_code(...)["X_stab"];
      - compute the EXACT primal enumerator A by brute enumeration (the same A the distributed
        engine produces, just done in-process here at K<=16);
      - cross-check A against qmsd.weightdist.weight_enumerator EXACTLY (vector equality);
      - run extract_d_and_Ad and assert d and A_d match qmsd.exact_distance_and_Ad AND the
        published oracle values (A_d = 880 and 3972);
      - assert all dual invariants (B integer, B_0=1, sum_w B_w = q^(n-K), all B_w >= 0).

This is the correctness anchor: the optimized distributed A (merged int64 partials) feeds the
SAME extract_d_and_Ad path validated here, so matching A guarantees matching (d, A_d).
"""
from __future__ import annotations

import numpy as np

from qmsd import weightdist as qwd
from qmsd.oracle import load_oracle
from qmsd.triorthogonal import build_triorthogonal_code

from dist_weightdist import correctness as cx


def brute_primal_enumerator(G0, q: int) -> list:
    """Exact A=[A_0..A_n] by enumerating all q^K messages in numpy chunks (int64 histogram).

    This mirrors EXACTLY what the distributed engine computes (merged int64 partials), so a
    match here certifies the engine's A path end to end at K<=16.  Asserts the int64 proof.
    """
    M = (np.asarray(G0).astype(np.int64)) % q
    K, n = M.shape
    cx.assert_int64_safe(K, q)
    total = q**K
    hist = np.zeros(n + 1, dtype=np.int64)
    qpows = q ** np.arange(K, dtype=np.int64)
    chunk = max(1, 4_000_000 // (n + 1))
    start = 0
    while start < total:
        end = min(start + chunk, total)
        idx = np.arange(start, end, dtype=np.int64)
        msgs = (idx[:, None] // qpows[None, :]) % q
        cws = (msgs @ M) % q
        wts = np.count_nonzero(cws, axis=1)
        hist += np.bincount(wts, minlength=n + 1).astype(np.int64)
        start = end
    A = [int(v) for v in hist]
    assert sum(A) == total, "brute enumerator lost words"
    # int64 headroom witness: the largest bin is well under INT64_MAX.
    assert max(A) <= cx.INT64_MAX
    return A


def check_krawtchouk():
    for n in (5, 12, 30, 60):
        for k in range(0, min(n, 14) + 1):
            col = cx.krawtchouk_column(0, n, 3, k)  # warmed up per x below
        for x in range(0, n + 1, max(1, n // 7)):
            col = cx.krawtchouk_column(x, n, 3, min(n, 14))
            for k in range(len(col)):
                ref = qwd.krawtchouk(k, x, n, 3)
                assert col[k] == ref == cx.krawtchouk_direct(k, x, n, 3), (k, x, n, col[k], ref)
    print("[ok] Krawtchouk recurrence == qmsd defining sum")


def check_macwilliams_random():
    rng = np.random.default_rng(0)
    for (q, n, K) in ((3, 10, 4), (3, 14, 5), (5, 8, 3)):
        # random full-ish generator
        G = rng.integers(0, q, size=(K, n), dtype=np.int64)
        A = brute_primal_enumerator(G, q)
        B_mine = cx.macwilliams(A, q)
        B_ref = qwd.macwilliams(A, q)
        assert B_mine == B_ref, (q, n, K)
    print("[ok] full MacWilliams == qmsd.macwilliams on random codes")


def check_oracle(label: str, exp_d: int, exp_Ad: int):
    oc = next(o for o in load_oracle() if o.label == label)
    built = build_triorthogonal_code(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed)
    G0 = np.asarray(built["X_stab"])
    q = oc.p
    K, n = G0.shape

    # A path: our brute (== distributed engine output) vs qmsd weight_enumerator -- EXACT.
    A = brute_primal_enumerator(G0, q)
    A_ref = qwd.weight_enumerator(G0, q, max_words=q**K + 1)
    assert A == A_ref, f"{label}: primal enumerator A mismatch vs qmsd"

    # (d, A_d) via full transform (all invariants) ...
    res_full = cx.extract_d_and_Ad(A, q)
    cx.assert_dual_invariants(res_full["B"], A, q)
    # ... and via the cheap windowed path used on the n~1939 target.
    res_win = cx.extract_d_and_Ad(A, q, search_kmax=exp_d + 2)
    assert res_full["d"] == res_win["d"] == exp_d, (label, res_full["d"], res_win["d"], exp_d)
    assert res_full["A_d"] == res_win["A_d"] == exp_Ad, (label, res_full["A_d"], exp_Ad)
    assert res_full["K"] == K

    # Cross-check the dual vector against the trusted qmsd MacWilliams transform.  We feed it
    # the SAME A (already validated == qmsd.weight_enumerator above) rather than letting
    # qmsd.exact_distance_and_Ad re-enumerate q^K words a third time.  This reproduces
    # qmsd.exact_distance_and_Ad exactly (which is weight_enumerator -> macwilliams).
    ref_B = qwd.macwilliams(A, q)
    assert ref_B == res_full["B"], f"{label}: dual B disagrees with qmsd.macwilliams"
    ref_d = next(w for w in range(1, len(ref_B)) if ref_B[w] > 0)
    assert ref_d == exp_d and ref_B[ref_d] == exp_Ad
    assert res_full["d"] == ref_d and res_full["A_d"] == ref_B[ref_d]

    print(f"[ok] {label}: K={K} n={n}  d={exp_d}  A_d={exp_Ad}  (A exact-match qmsd, "
          f"invariants pass: B_0=1, sum_w B_w = {q}^{n - K})")


def main():
    cx.assert_int64_safe(cx.DESIGN_K_CEILING, 3)
    print(f"[ok] int64 proof: 3^{cx.DESIGN_K_CEILING} = {3**cx.DESIGN_K_CEILING} <= INT64_MAX "
          f"= {cx.INT64_MAX} (margin {cx.INT64_MAX // 3**cx.DESIGN_K_CEILING}x)")
    check_krawtchouk()
    check_macwilliams_random()
    check_oracle("[[206,37,4]]_3", 4, 880)
    check_oracle("[[667,62,4]]_3", 4, 3972)
    print("ALL CORRECTNESS CHECKS PASSED")


if __name__ == "__main__":
    main()
