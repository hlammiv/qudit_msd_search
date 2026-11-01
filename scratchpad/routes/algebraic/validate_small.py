"""Brute-force validation of the hypersurface theorem at small parameters.

Theorem: O = Z(g) (ALL zeros of g, deg g = delta) with RM(r,m) of full rank on O gives
d(ker G0) >= d(RM(rtilde + delta, m)).  Block version: O a full-rank subset of Z(g),
B = Z(g) \\ O: every kernel word meeting the complement of Z(g) has weight >= that bound,
so d = min(bound, d_block) with d_block = min dependent columns of G0 restricted to B.
"""
import numpy as np
from qmsd.reedmuller import rm_generator
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.mindist import min_dependent_columns
from hyp_common import random_poly_values, rank_modp, bound


def mdc(H, p):
    H = np.asarray(H, dtype=np.int64) % p
    if H.shape[1] == 0:
        return np.inf
    if H.shape[0] == 0:
        return 1
    try:
        return min_dependent_columns(H, p)
    except ValueError:
        return 7  # ">= 7" (MITM cap)


def run(p, m, r, delta, trials, rng):
    G = rm_generator(r, m, p)
    Gi = np.asarray(G, dtype=np.int64)
    K = Gi.shape[0]
    cap = K - rm_generator(r - delta, m, p).shape[0] if r - delta >= 0 else K
    b = bound(delta, r, m, p)
    print(f"\n=== p={p} m={m} r={r} deg g={delta}: dim RM={K}, cap on k = {cap}, theorem bound d >= {b}", flush=True)
    for t in range(trials):
        g, c = random_poly_values(delta, m, p, rng, homogeneous_top=True)
        Z = np.nonzero(g == 0)[0]
        if len(Z) == 0:
            continue
        h = rank_modp(Gi[:, Z], p)
        if h == len(Z):
            code = build_triorthogonal_code(p, m, r, list(Z + 1), G=G)
            assert code["full_rank"]
            d = mdc(code["G0"], p)
            print(f"  trial {t}: |Z|={len(Z)} h_Z={h} INDEPENDENT -> O=Z(g): [[{code['n']},{code['k']},d]] "
                  f"exact d={d}{' (>=7, MITM cap)' if d == 7 else ''}; bound {b}: {'OK' if d >= min(b, 7) else 'VIOLATION'}",
                  flush=True)
        else:
            rho = int(rng.integers(0, 4))
            perm = rng.permutation(Z)
            O, cur = [], 0
            for z in perm:
                if rank_modp(Gi[:, O + [int(z)]], p) > cur:
                    O.append(int(z)); cur += 1
                if cur == h - rho:
                    break
            B = sorted(set(Z.tolist()) - set(O))
            code = build_triorthogonal_code(p, m, r, [o + 1 for o in O], G=G)
            assert code["full_rank"]
            G0 = np.asarray(code["G0"], dtype=np.int64)
            nonout = np.array(sorted(set(range(p ** m)) - set(O)))
            colB = np.searchsorted(nonout, B)
            d = mdc(G0, p)
            dB = mdc(G0[:, colB], p)
            pred = min(b, dB)
            if pred >= 7:
                ok = d >= 7
            else:
                ok = (d == pred)
            print(f"  trial {t}: |Z|={len(Z)} h_Z={h} dependent -> O={len(O)} pts, block |B|={len(B)} rho={rho}: "
                  f"d={d}, d_block={dB}, predicted min(bound,d_block)={pred}: {'OK' if ok else 'MISMATCH'}", flush=True)


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    run(7, 3, 5, 2, 10, rng)
    run(7, 3, 5, 3, 10, rng)
    run(7, 3, 5, 1, 3, rng)
    run(5, 4, 5, 3, 6, rng)
    run(5, 4, 5, 2, 5, rng)
    run(5, 4, 5, 4, 5, rng)
