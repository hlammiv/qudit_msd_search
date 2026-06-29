"""Collect .partial files -> exact merged primal enumerator A -> certified (d, A_d, B).

This is the one-time, single-machine tail of the pipeline (cheap vs the 3^K
enumeration).  It:
  1. loads every .partial, verifying each one's internal checksum;
  2. merges them with ``correctness.merge_partials(expected_total=3^K)`` which
     fails LOUDLY on any lost/double-counted message (global checksum, design sec.6.3);
  3. NEVER transforms an A whose sum != 3^K;
  4. runs the FULL bignum q-ary MacWilliams (asserting sum_w B_w = 3^(n-K) and all
     dual invariants) and extracts d = min{w>0: B_w>0}, A_d = B_d.
"""
from __future__ import annotations

import json
from pathlib import Path

from . import correctness as cx
from .partial import read_partial


def collect_partials(paths) -> tuple[list[list[int]], int, int]:
    """Read .partial files; assert consistent (q,K,n); return (hists, n, K)."""
    hists, n, K = [], None, None
    for path in paths:
        p = read_partial(path)
        if n is None:
            n, K = p["n"], p["K"]
        assert p["n"] == n and p["K"] == K and p["q"] == 3, f"{path}: mismatched header"
        hists.append(p["hist"])
    assert n is not None, "no partials provided"
    return hists, n, K


def merge_to_A(paths) -> tuple[list[int], int, int]:
    """Merge .partial files into the exact primal enumerator A; verify global checksum 3^K."""
    hists, n, K = collect_partials(paths)
    cx.assert_int64_safe(K, 3)
    expected_total = 3 ** K
    A = cx.merge_partials(hists, n, expected_total=expected_total)
    # Ingest guard (design sec.3.1): |C| = sum(A) must be a power of 3.
    s = sum(A)
    p, e = 1, 0
    while p < s:
        p *= 3; e += 1
    assert p == s, f"|C| = sum(A) = {s} is not a power of 3 (rank-deficient/malformed G0)"
    return A, n, K


def certify(paths, search_kmax: int | None = None) -> dict:
    """Full pipeline tail: merge partials -> A -> (d, A_d, full B) with all invariants.

    When ``search_kmax`` is given, a fast WINDOWED pre-flight runs first (raises if true
    d exceeds the window, never reporting a spuriously small d); certification is always
    the FULL transform on the merged A.
    """
    A, n, K = merge_to_A(paths)
    if search_kmax is not None:
        win = cx.extract_d_and_Ad(A, 3, search_kmax=search_kmax)  # pre-flight; raises if d>window
    res = cx.extract_d_and_Ad(A, 3)                # full transform, all global invariants
    cx.assert_dual_invariants(res["B"], A, 3)
    res["A"] = A
    res["n"], res["K"] = n, K
    if search_kmax is not None:
        assert win["d"] == res["d"] and win["A_d"] == res["A_d"], "windowed vs full mismatch"
    return res


def write_outputs(res: dict, out_dir: str) -> dict:
    """Write A.txt (one count per line) and result.json (d, A_d, B, header)."""
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    A, B = res["A"], res["B"]
    (out / "A.txt").write_text("\n".join(str(a) for a in A) + "\n")
    (out / "B.txt").write_text("\n".join(str(b) for b in B) + "\n")
    summary = {"n": res["n"], "K": res["K"], "q": 3,
               "d": res["d"], "A_d": res["A_d"],
               "size": sum(A), "sum_B": sum(B)}
    (out / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
