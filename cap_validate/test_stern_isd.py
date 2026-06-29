"""Self-test / correctness anchors for the Stern/ISD finder (Method 1).

Cross-checks Stern against the EXACT engine ``qmsd.weightdist.exact_distance_and_Ad``
(MacWilliams on the small generator) on codes whose true minimum distance is known.
For every such code, Stern MUST find a codeword of the exact minimum weight.
"""
import sys
import numpy as np

import os
sys.path.insert(0, os.getcwd())          # project root (for qmsd)
sys.path.insert(0, os.path.dirname(__file__))  # cap_validate (for stern_isd)
from stern_isd import stern_search, stern_search_multi  # noqa: E402

from qmsd.weightdist import exact_distance_and_Ad  # noqa: E402
from qmsd import reedmuller, triorthogonal  # noqa: E402


def exact_dual_distance(H):
    res = exact_distance_and_Ad(np.asarray(H, dtype=np.int64) % 3, 3, max_words=6_000_000)
    assert res["feasible"], "exact engine infeasible for this size"
    return res["distance"]


def check(name, H, **kw):
    H = np.asarray(H, dtype=np.int8) % 3
    # ensure full row rank for Stern
    r0 = int(np.linalg.matrix_rank(H.astype(np.int64)))
    if r0 != H.shape[0]:
        # drop dependent rows
        from qmsd.field import GFp
        rows = GFp(3)(H % 3).row_reduce()
        H = np.asarray(rows, dtype=np.int8)
        H = H[np.any(H, axis=1)]
    d_exact = exact_dual_distance(H)
    res = stern_search_multi(H, target=None, **kw)
    ok = (res.best_weight == d_exact) and res.verify(H)
    print(
        f"[{name}] H={H.shape}  exact d={d_exact}  stern best={res.best_weight}  "
        f"valid={res.verify(H)}  iters={res.iterations}  {res.elapsed:.2f}s  "
        f"{'PASS' if ok else 'FAIL'}"
    )
    if res.witness is not None:
        supp = np.nonzero(res.witness)[0].tolist()
        print(f"       witness weight {int(np.count_nonzero(res.witness))} support {supp}")
    return ok


def main():
    results = []

    # --- 1. random small full-rank parity checks (varied min distances) ---
    for trial, (rr, nn, seed) in enumerate([(8, 40, 1), (9, 55, 2), (10, 60, 3)]):
        rng = np.random.default_rng(seed)
        H = rng.integers(0, 3, size=(rr, nn)).astype(np.int8)
        results.append(
            check(
                f"random{trial} {rr}x{nn}",
                H,
                p_values=(1, 2, 3),
                l=min(rr, 7),
                iterations=900,
                weight_budget=nn,
                threads=2,
                seed=100 * trial,
            )
        )

    # --- 2. a real small triorthogonal cap-style code (same builder as the target) ---
    # m=4, p=3, r=2 (satisfies 3r<m(p-1): 6<8): RM_3(2,4) punctured -> small G0,
    # whose dual is a punctured RM_3 code -- the exact analogue of the cap target.
    mm, rr = 4, 2
    G = reedmuller.rm_generator(rr, mm, 3)
    Gnp = np.asarray(G) % 3
    n0 = Gnp.shape[1]
    rng = np.random.default_rng(11)
    from qmsd.field import GFp
    npunc = 3
    while True:
        S0 = sorted(rng.choice(n0, size=npunc, replace=False).tolist())
        sub = GFp(3)(Gnp[:, S0])
        if int(np.linalg.matrix_rank(sub)) == npunc:
            break
    S1 = [c + 1 for c in S0]  # 1-indexed
    built = triorthogonal.build_triorthogonal_code(3, mm, rr, S1, G=G)
    G0 = np.asarray(built["X_stab"], dtype=np.int8) % 3
    print(f"triorthogonal small code: G0 {G0.shape}, full_rank={built['full_rank']}")
    results.append(
        check(
            "triortho m=4 r=2",
            G0,
            p_values=(1, 2),
            l=min(G0.shape[0], 8),
            iterations=600,
            weight_budget=G0.shape[1],
            threads=2,
            seed=999,
        )
    )

    # --- 3. p=1 sublist path on a code with a weight-2 codeword (2 dependent cols) ---
    H = np.array(
        [[1, 0, 0, 1, 0, 2],
         [0, 1, 0, 0, 1, 0],
         [0, 0, 1, 0, 0, 1]],
        dtype=np.int8,
    )  # cols 2 and 5 are equal -> weight-2 codeword exists
    results.append(
        check("p=1 weight2", H, p_values=(1,), l=3, iterations=200, weight_budget=6,
              threads=1, seed=5)
    )

    print()
    if all(results):
        print(f"ALL {len(results)} CHECKS PASSED")
        return 0
    print(f"{sum(results)}/{len(results)} passed -- FAILURES PRESENT")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
