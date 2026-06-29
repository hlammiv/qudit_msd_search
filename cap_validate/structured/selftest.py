"""Self-test for Method 2 (structured_enum).

Asserts, with no network and no full 2.9e9 materialisation:

  A. BUILD: G0 is 55x1968, full rank; dual code is RM_3(9,7), d_RM=18.
  B. STRUCTURE: weight-18 codewords are two parallel 2-flats.  For sampled directions
     and coset pairs we BUILD the explicit codeword v (1 on coset A, 2 on coset B) and
     assert (i) it is a genuine RM_3(9,7) codeword (G v = 0), (ii) its support is exactly
     the 18 points of the two cosets.  (Closes "family support <=> genuine codeword".)
  C. COUNT: number of dim-2 directions = 99463 and the full weight-18 support count =
     99463 * C(243,2) = 2,924,510,589 (matches prior work ~2.9e9).
  D. KERNEL == BRUTE: for sampled directions, the additive top-2 (numba max_topj) equals
     the brute max over all C(243,2) explicit coset-pair supports (numba inter_counts).
  E. MIN: rigorous min |supp \\ S| over the whole weight-18 family = 10, and over the
     weight-27 (affine 3-flat) family.

Run:  python -m cap_validate.structured.selftest
"""
from __future__ import annotations

from itertools import combinations
from math import comb

import numpy as np

import qmsd.reedmuller as rm

from . import capcode as cc
from . import enum_families as ef
from . import subspaces as sub
from .kernels import inter_counts, max_topj


def _G_dual():
    """Dual generator RM_3(4,7): v in RM_3(9,7) iff G_dual @ v == 0 (mod 3)."""
    return np.asarray(rm.rm_generator(cc.R, cc.M, cc.P)).astype(np.int64) % cc.P


def _coset_ids(H, pts):
    """coset id in [0,3^codim) for every point row of pts under syndrome map H."""
    codim = H.shape[0]
    syn = (pts.astype(np.int64) @ H.T) % cc.P          # (Npts, codim)
    pow3 = (cc.P ** np.arange(codim)).astype(np.int64)
    return syn @ pow3


def run(verbose=True):
    rep = {}

    # ---- A. build -------------------------------------------------------------
    b = cc.verify_build()
    assert b["G0_shape"] == (55, 1968), b
    assert b["full_rank"] and b["dual_code"] == "RM_3(9,7)" and b["d_rm_min"] == 18
    rep["A_build"] = b

    cap = cc.load_cap()
    Spts = cap["cap_pts"]
    pts = cap["all_pts"].astype(np.int64)
    assert cap["size"] == 219
    Gd = _G_dual()

    # ---- C. counts ------------------------------------------------------------
    n_dir2 = sub.gaussian_binomial(7, 2)
    n_supp18 = n_dir2 * comb(3 ** 5, 2)
    assert n_dir2 == 99463, n_dir2
    assert n_supp18 == 2_924_510_589, n_supp18
    rep["C_counts"] = {"n_dir_dim2": n_dir2, "coset_pairs": comb(243, 2),
                       "n_weight18_supports": n_supp18}

    # ---- build all dim-2 directions once (shared by B,D,E) --------------------
    Hs2 = sub.syndrome_maps(2)                 # (99463, 5, 7)
    assert Hs2.shape == (99463, 5, 7)

    # ---- B. structure: sampled family supports are genuine weight-18 codewords -
    rng = np.random.default_rng(7)
    struct_checks = 0
    for s in rng.choice(Hs2.shape[0], size=12, replace=False):
        H = Hs2[s]
        cid_all = _coset_ids(H, pts)           # (2187,)
        # pick two distinct occupied cosets (all 243 are occupied: 2187/243 = 9 pts each)
        cosets = np.unique(cid_all)
        cA, cB = rng.choice(cosets, size=2, replace=False)
        A = np.where(cid_all == cA)[0]
        B = np.where(cid_all == cB)[0]
        assert len(A) == 9 and len(B) == 9      # each 2-flat has 3^2 = 9 points
        v = np.zeros(cc.N, dtype=np.int64)
        v[A] = 1
        v[B] = 2
        # (i) genuine codeword
        assert int(np.abs((Gd @ v) % cc.P).sum()) == 0, "not a RM_3(9,7) codeword"
        # (ii) support = the 18 points, and they are two parallel 2-flats (same dir V)
        supp = np.where(v != 0)[0]
        assert len(supp) == 18
        assert set(supp.tolist()) == set(A.tolist()) | set(B.tolist())
        struct_checks += 1
    rep["B_structure"] = {"sampled_directions": struct_checks,
                          "each": "v(1 on A,2 on B) is a genuine wt-18 codeword; "
                                  "support = two parallel 2-flats"}

    # ---- D. additive kernel == brute over explicit supports -------------------
    # For a handful of directions, compare max_topj (top-2) to the brute max over all
    # C(243,2) explicit 18-point coset-pair supports.
    mask = np.zeros(cc.N, dtype=np.uint8)
    mask[cap["idx0"]] = 1
    d_dirs = rng.choice(Hs2.shape[0], size=6, replace=False)
    Hs_sample = np.ascontiguousarray(Hs2[d_dirs])
    add_top2 = max_topj(Hs_sample, Spts, 5, 2)        # additive, per direction
    brute_ok = True
    for li, s in enumerate(d_dirs):
        H = Hs2[s]
        cid_all = _coset_ids(H, pts)
        # point-id lists per coset (243 cosets, 9 pts each)
        buckets = [np.where(cid_all == c)[0] for c in range(243)]
        # explicit supports for every coset pair
        pairs = list(combinations(range(243), 2))
        supports = np.empty((len(pairs), 18), dtype=np.int64)
        for pi, (i, jc) in enumerate(pairs):
            supports[pi, :9] = buckets[i]
            supports[pi, 9:] = buckets[jc]
        inter = inter_counts(supports, mask)
        brute_max = int(inter.max())
        if brute_max != int(add_top2[li]):
            brute_ok = False
            break
    assert brute_ok, "additive top-2 != brute over explicit supports"
    rep["D_kernel_vs_brute"] = {"directions_checked": len(d_dirs),
                                "agree": True,
                                "note": "additive coset top-2 == brute max over "
                                        "C(243,2) explicit 18-pt supports"}

    # ---- E. rigorous family minima -------------------------------------------
    w18 = ef.family_min_punct(ef.FAMILIES["w18"], Spts, Hs2)
    assert w18["min_punct"] == 10, w18
    assert w18["n_supports"] == 2_924_510_589
    w27 = ef.family_min_punct(ef.FAMILIES["w27"], Spts)
    rep["E_minima"] = {"w18": w18, "w27": w27}
    rep["rigorous_min_over_enumerated"] = min(w18["min_punct"], w27["min_punct"])

    if verbose:
        import json
        print(json.dumps(rep, indent=2, default=int))
    return rep


if __name__ == "__main__":
    run()
