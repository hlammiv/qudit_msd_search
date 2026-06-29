"""M3 — scalable ring-triorthogonality (R1) construction at LARGE block length (n = 9, 11, 13, ...).

The d=4 obstruction (M2/M2b/R1) is bounded to n<=7 only because those searches enumerate the d^n grid
(transversality region V and brute duals).  M3 removes that ceiling: nothing here touches d^n.

  * X-stabilizers M_X and the logical module L = M_X + <g> are SMALL (|L| <= d*|M_X| <= 64), so
    transversality is checked directly on L (Phi constant on each M_X-coset of L).
  * M_Z = L^perp and M_X^perp are obtained as GENERATORS via the scalable Howell kernel (no d^n).
  * distance >= 2 == "no weight-1 logical operator" — checked over the n*(d-1) weight-1 vectors via
    membership, not by min-distance search.

So n is limited only by O(n) bookkeeping; n = 13, 15 are trivial.  This is the M3 substrate + the
larger-n R1 test the d=4 obstruction left open.

Run:  python -m primepower_msd.m3_search
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

N_VALUES = (9, 11, 13)


def _order_mod(g, mx_set, d):
    """Order of coset g in Z_d^n / M_X (smallest t>=1 with t*g in M_X), or None."""
    gv = np.asarray(g) % d
    for t in range(1, d + 1):
        if tuple(int(x) for x in (t * gv) % d) in mx_set:
            return t
    return None


def m3_worker(seed: int, trials: int):
    """Sample (M_X, high-weight g); construct L=M_X+<g>; test transversal + distance>=2 + level-3."""
    import numpy as _np
    from primepower_msd.ringlinalg import enumerate_module, right_kernel_howell, in_span
    from primepower_msd.ring_transversal import extract_phase
    from primepower_msd.clifford_ring import is_level3, is_clifford
    from primepower_msd.hierarchy_search import search as _level3_family

    d = 4
    rng = _np.random.default_rng(seed)
    fam = _level3_family(d)
    gate_phis = [(gi, *extract_phase(U, d), U) for gi, (a, b, s, U) in enumerate(fam)]
    witnesses = []
    scanned = 0
    self_ok = 0
    transversal = 0

    def Phi(state, phi, N):
        return int(sum(int(phi[int(x)]) for x in state) % N)

    for n in N_VALUES:
        for gi, phi_list, N, U in gate_phis:
            phi = _np.asarray(phi_list, dtype=_np.int64)
            for _ in range(trials):
                rx = int(rng.integers(1, 3))
                A = rng.integers(0, d, size=(rx, n))
                mx = enumerate_module(A, d)
                mx_set = set(mx)
                # (self) Phi vanishes on M_X
                if any(Phi(a, phi, N) != 0 for a in mx):
                    continue
                self_ok += 1
                # sample a high-weight coset generator g (weight>=2 => no trivial weight-1 X-logical)
                wt = int(rng.integers(2, n + 1))
                pos = rng.choice(n, size=wt, replace=False)
                g = _np.zeros(n, dtype=_np.int64)
                g[pos] = rng.integers(1, d, size=wt)
                if tuple(int(x) for x in g) in mx_set or _order_mod(g, mx_set, d) != d:
                    continue
                cosets = [tuple(int(x) for x in (t * g) % d) for t in range(d)]
                # (transversal) Phi constant on each M_X-coset of L  (check only on small L)
                ok = True
                for r in cosets:
                    base = Phi(r, phi, N)
                    for a in mx:
                        st = tuple((r[j] + a[j]) % d for j in range(n))
                        if Phi(st, phi, N) != base:
                            ok = False
                            break
                    if not ok:
                        break
                if not ok:
                    continue
                transversal += 1
                # L = M_X + <g>
                Lgen = _np.vstack([A, g])
                L = enumerate_module(Lgen, d)
                L_set = set(L)
                # X-distance >= 2 : no weight-1 vector in L \ M_X
                if any(sum(1 for x in v if x) == 1 for v in L if v not in mx_set):
                    continue
                # Z-distance >= 2 : no weight-1 v in M_X^perp \ M_Z   (M_Z = L^perp)
                mz_gens = right_kernel_howell(Lgen, d)
                bad = False
                for i in range(n):
                    for w in range(1, d):
                        # v = w e_i in M_X^perp  <=>  w * a_i == 0 (mod d) for all a in M_X
                        if all((w * a[i]) % d == 0 for a in mx):
                            v = [0] * n
                            v[i] = w
                            if not in_span(v, mz_gens, d):      # in M_X^perp but not in M_Z -> weight-1 Z-logical
                                bad = True
                                break
                    if bad:
                        break
                if bad:
                    continue
                scanned += 1
                # induced logical gate must be strict level 3
                lam0 = Phi(cosets[0], phi, N)
                lam = [(Phi(cosets[ell], phi, N) - lam0) % N for ell in range(d)]
                UL = _np.diag([_np.exp(2j * _np.pi * (l / N)) for l in lam])
                if is_level3(UL, d) and not is_clifford(UL, d):
                    if len(witnesses) < 40:
                        witnesses.append({"n": n, "A": A.tolist(), "g": [int(x) for x in g],
                                          "gi": gi, "lam": lam, "N": N})
    return {"witnesses": witnesses, "scanned": scanned, "self_ok": self_ok, "transversal": transversal}


def main() -> int:
    import time
    from primepower_msd.parallel import run_chunks, default_jobs
    from primepower_msd.hierarchy_search import search as level3_family

    d = 4
    n_jobs = default_jobs()
    print("=" * 78)
    print(f"M3 — scalable ring-triorthogonality (R1) at large n {N_VALUES} (d={d}); no d^n enumeration")
    print("=" * 78)
    t0 = time.time()
    print(f"  {n_jobs} workers x {n_jobs} chunks x 3000 (M_X,g) samples ...", flush=True)
    results = run_chunks(m3_worker, n_chunks=n_jobs, trials_per_chunk=3000, n_jobs=n_jobs)
    scanned = sum(r["scanned"] for r in results)
    self_ok = sum(r["self_ok"] for r in results)
    transversal = sum(r["transversal"] for r in results)
    witnesses = [w for r in results for w in r["witnesses"]]
    per_n = {}
    for w in witnesses:
        per_n[w["n"]] = per_n.get(w["n"], 0) + 1

    print("=" * 78)
    print("VERDICT (M3 — larger n)")
    print("=" * 78)
    print(f"  {self_ok} self-orthogonal M_X; {transversal} transversal constructions; "
          f"{scanned} distance>=2 transversal codes in {time.time()-t0:.1f}s ({n_jobs} cores).")
    print(f"  WITNESSES (distance>=2 + strict level-3 induced): {len(witnesses)} found; per-n: {per_n}")
    if witnesses:
        w = min(witnesses, key=lambda x: x["n"])      # smallest-n witness
        ga, gb, _, _ = level3_family(d)[w["gi"]]
        print(f"\n  *** SMALLEST WITNESS at n={w['n']}: [[{w['n']},1,>=2]]_{d} ring-triorthogonal CSS, "
              f"transversal antidiff(S^{ga} Z^{gb}) ***")
        print(f"    X-stabilizers: {w['A']}")
        print(f"    logical X gen g: {w['g']}")
        print(f"    induced logical phases (units 2pi/{w['N']}): {w['lam']}  (strict level-3)")
        print(f"    -> GENUINE single-ququart distillation code (distance>=2) at large n CONFIRMED.")
        return 0
    print(f"\n  No distance>=2 transversal level-3 code found at n in {N_VALUES}.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
