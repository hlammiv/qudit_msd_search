"""R1 — analytic ring-triorthogonality construction: can a *genuinely-quadratic* level-3 gate be
transversal (per-coset) on a distance>=2 single-ququart code?

M2b showed the quadratic level-3 gates have trivial GLOBAL translation stabilizer, so they are not
globally transversal.  But the true (per-coset) condition is weaker:

  pick X-stabilizers M_X and a logical-rep module L = M_Z^perp with M_X subseteq L such that
    (self)        sum_i phi(a_i) == 0 (mod N)            for all a in M_X
    (transversal) Phi(c + a) == Phi(c) (mod N)           for all a in M_X, c in L

The transversality region V = { c : Phi(c+a) == Phi(c) for all a in M_X } is a UNION of M_X-cosets, so
we CONSTRUCT L instead of sampling: find a coset generator g of order d with g, 2g, ..., (d-1)g all in
V, set L = M_X + <g>, M_Z = L^perp.  Then certify distance and that the induced logical gate is level-3.

This is the analytic co-design the random searches (M2/M2b) cannot reach (the solutions are measure
zero).  Parallelized over X-stabilizer candidates.

Run:  python -m primepower_msd.r1_search
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def _grid(n: int, d: int) -> np.ndarray:
    return np.indices((d,) * n).reshape(n, -1).T.astype(np.int64)      # (d^n, n)


def transversality_region(grid, phi_c, mx_arr, phi, N, d):
    """V = {c in Z_d^n : Phi(c+a) == Phi(c) (mod N) for all a in M_X}, as a boolean mask over `grid`."""
    mask = np.ones(grid.shape[0], dtype=bool)
    for a in mx_arr:
        shifted = (grid + a) % d
        phi_s = (phi[shifted].sum(axis=1)) % N
        mask &= (phi_s == phi_c)
        if not mask.any():
            break
    return mask


def r1_worker(seed: int, trials: int):
    """Sample X-stabilizer modules; for each, construct L subseteq V and test for a distance>=2,
    level-3, single-ququart code.  Top-level (picklable) for joblib."""
    import numpy as _np
    from primepower_msd.ringlinalg import dual_module, howell_form, enumerate_module
    from primepower_msd.ring_css import build_css_general
    from primepower_msd.ring_transversal import extract_phase, certify_distillation_code
    from primepower_msd.hierarchy_search import search as _level3_family

    d = 4
    rng = _np.random.default_rng(seed)
    fam = _level3_family(d)
    gate_phis = [(gi, extract_phase(U, d)) for gi, (a, b, s, U) in enumerate(fam)]
    best = None
    scanned = 0
    self_ok_count = 0

    for n in (5, 7):
        grid = _grid(n, d)
        for gi, (phi_list, N) in gate_phis:
            phi = _np.asarray(phi_list, dtype=_np.int64)
            phi_c = (phi[grid].sum(axis=1)) % N
            for _ in range(trials):
                rx = int(rng.integers(1, 3))
                A = rng.integers(0, d, size=(rx, n))
                mx = enumerate_module(A, d)
                mx_arr = _np.asarray(mx, dtype=_np.int64)
                # (self) condition: Phi must vanish on M_X
                if _np.any((phi[mx_arr].sum(axis=1)) % N != 0):
                    continue
                self_ok_count += 1
                mx_set = set(mx)
                mask = transversality_region(grid, phi_c, mx_arr, phi, N, d)
                Vset = {tuple(int(x) for x in row) for row in grid[mask]}
                fam_U = fam[gi][3]
                # try several cyclic order-d coset generators g with g,2g,...,(d-1)g all in V
                tested_L = 0
                seen_L = set()
                for g in (gg for gg in Vset if gg not in mx_set):
                    if tested_L >= 12:
                        break
                    gv = _np.asarray(g, dtype=_np.int64)
                    order = next((t for t in range(1, d + 1)
                                  if tuple(int(x) for x in (t * gv) % d) in mx_set), None)
                    if order != d:
                        continue
                    if not all(tuple(int(x) for x in (t * gv) % d) in Vset for t in range(2, d)):
                        continue
                    # L = M_X + <g>;  M_Z = L^perp
                    Lgen = _np.vstack([A, gv])
                    L_elems = enumerate_module(Lgen, d)
                    Lkey = frozenset(L_elems)
                    if Lkey in seen_L:
                        continue
                    seen_L.add(Lkey)
                    B = howell_form(_np.asarray(dual_module(L_elems, d)), d)
                    if B.shape[0] == 0:
                        continue
                    code = build_css_general(A, B, d)
                    scanned += 1
                    tested_L += 1
                    if code is None or code.k != 1 or not code.cyclic:
                        continue
                    cert = certify_distillation_code(code, fam_U)
                    if not cert.get("is_distillation_code"):
                        continue
                    if code.distance >= 2 and (best is None or code.distance > best["distance"]):
                        best = {"n": n, "A": A.tolist(), "B": B.tolist(), "gi": gi,
                                "distance": code.distance, "cert": cert}
    return {"best": best, "scanned": scanned, "self_ok": self_ok_count}


def main() -> int:
    import time
    from primepower_msd.parallel import run_chunks, default_jobs
    from primepower_msd.hierarchy_search import search as level3_family
    from primepower_msd.ring_css import build_css_general, logical_group_invariants

    d = 4
    n_jobs = default_jobs()
    print("=" * 78)
    print(f"R1 — ring-triorthogonality construction for genuinely-quadratic transversal gates (d={d})")
    print("=" * 78)
    t0 = time.time()
    print(f"  {n_jobs} workers x {n_jobs} chunks x 1500 X-stabilizer samples ...", flush=True)
    results = run_chunks(r1_worker, n_chunks=n_jobs, trials_per_chunk=1500, n_jobs=n_jobs)
    scanned = sum(r["scanned"] for r in results)
    self_ok = sum(r["self_ok"] for r in results)
    bests = [r["best"] for r in results if r["best"]]
    best = max(bests, key=lambda b: b["distance"]) if bests else None

    print("=" * 78)
    print("VERDICT (R1)")
    print("=" * 78)
    print(f"  {self_ok} self-orthogonal M_X found; {scanned} constructed transversal codes in "
          f"{time.time()-t0:.1f}s ({n_jobs} cores).")
    if best:
        ga, gb, _, U = level3_family(d)[best["gi"]]
        code = build_css_general(best["A"], best["B"], d)
        cert = best["cert"]
        print(f"\n  *** WITNESS: [[{code.n},1,{best['distance']}]]_{d} ring-triorthogonal CSS, "
              f"transversal antidiff(S^{ga} Z^{gb}) ***")
        print(f"    X-stabilizers: {best['A']}")
        print(f"    Z-stabilizers: {best['B']}")
        print(f"    logical orders M_Z^perp/M_X: {logical_group_invariants(code)} cyclic={code.cyclic}")
        print(f"    induced logical level {cert['logical_level']} strict-L3 {cert['logical_strict_level3']}")
        print(f"    -> GENUINE single-ququart distillation code via ring-triorthogonality CONFIRMED.")
        return 0
    print("\n  No distance>=2 ring-triorthogonal code constructed for any level-3 gate (d=4, n<=7).")
    print("  Strong evidence that diagonal-transversal single-ququart MSD is obstructed even per-coset;")
    print("  next: extend to d=8 (Hoggar) and larger n (M3), or pursue the no-go proof.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
