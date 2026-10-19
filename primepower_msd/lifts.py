"""Constructive code families for the diverse-gate search: Z_4 LIFTS of binary triorthogonal /
triply-even codes (the Haah-Hastings / Bravyi-Haah / Reed-Muller template), ZRM codes, Z-refinement
lattices, and residue-field bookkeeping.

Why lifts.  Random Z_4 codes never pass the q-weight screen S1 (Phase 3); binary triply-even and
triorthogonal matrices are the classical objects whose Z_4 lifts have a chance (their rows have
weights == 0 mod 8 / pairwise overlaps == 0 mod 4 / triple overlaps even, which are the mod-2 shadows of
S1-S3).  Lifting is not unique: the canonical 0/1 lift, lifts with extra 2-torsion rows (2 * binary
vectors) and the quaternary Reed-Muller ZRM(r, m) = <lift01(RM(r-1,m)), 2 RM(r,m)> (Hammons et al.)
all reduce to the same binary residue code.  Z-REFINEMENTS add 2 v rows to M_Z (shrinking the support,
lowering k) and are the knob that turns a non-hosting lift into a hosting one.

Residue bound (Gluesing-Luerssen-Pllaha 1710.09884): a FREE Z_4 code has distance <= its residue code's;
`residue_distance_bound` evaluates the binary bound with qmsd.mindist.min_dependent_columns(H, 2).
"""

from __future__ import annotations

import itertools as it

import numpy as np

from .ringlinalg import howell_form, module_size, right_kernel_howell, in_span
from .symplectic import Z4Code


# --------------------------------------------------------------------------- binary sources
def rm_binary(r: int, m: int) -> np.ndarray:
    """Generator of the binary Reed-Muller code RM(r, m) (rows = monomials of degree <= r), via qmsd."""
    from qmsd.reedmuller import rm_generator
    G = np.array(rm_generator(r, m, 2), dtype=np.int64) % 2
    return G


def punctured_rm_triorthogonal(m: int) -> dict:
    """The [[2^m - 1, 1]] punctured-RM(1,m) triorthogonal CSS code (m=4: the 15-qubit code)
    from qmsd.triorthogonal.build_triorthogonal_code(2, m, 1, [1]).  Returns binary X/Z stabilizers."""
    from qmsd.triorthogonal import build_triorthogonal_code
    b = build_triorthogonal_code(2, m, 1, [1])
    return {"n": b["n"], "k": b["k"], "MX": np.array(b["X_stab"], dtype=np.int64) % 2,
            "MZ": np.array(b["Z_stab"], dtype=np.int64) % 2, "label": f"punctured RM(1,{m}) [[{b['n']},1]]"}


def cube_8_3_2() -> dict:
    """The [[8,3,2]] cube code (Steane / Haah-Hastings m=3 CCZ code): X = all-ones, Z = the 6 faces."""
    n = 8
    pts = list(it.product(range(2), repeat=3))
    faces = []
    for axis in range(3):
        for val in range(2):
            faces.append([1 if p[axis] == val else 0 for p in pts])
    return {"n": n, "k": 3, "MX": np.ones((1, n), dtype=np.int64), "MZ": np.array(faces, dtype=np.int64),
            "label": "cube [[8,3,2]]"}


def random_binary_generalized_triorthogonal(n: int, k: int, rng, n_stab: int = 2, tries: int = 20000,
                                            kind: str = "T") -> dict | None:
    """Random binary generalized-triorthogonal matrix (Haah-Hastings Def. 1 / Eq. 2.7) with k logical
    rows G1 and n_stab stabilizer rows G0: all pairwise overlaps even (incl. self-overlaps of G0 rows =
    even weights), all triple overlaps even EXCEPT the diagonal ones prescribed by `kind`:
      kind="T":   each logical row has odd weight (self-triple 1) -> transversal T gives logical T^{(x)k}.
      kind="CCZ": k = 3, the three logical rows have odd 3-way overlap, even weights.
    Returns {"G1","G0","n","k"} or None.  Bravyi-Haah / HH small codes are found this way in seconds."""
    def triple(a, b, c):
        return int(np.sum(a & b & c)) % 2

    for _ in range(tries):
        G1 = rng.integers(0, 2, size=(k, n)).astype(np.int64)
        G0 = rng.integers(0, 2, size=(n_stab, n)).astype(np.int64)
        G = np.vstack([G1, G0])
        if np.any(~np.any(G, axis=0)):          # a zero column gives a weight-1 logical
            continue
        ok = True
        R = G.shape[0]
        for a in range(R):
            for b in range(a, R):
                if a == b and a < k and kind == "T":
                    if int(G[a].sum()) % 2 != 1:
                        ok = False; break
                    continue
                if int(np.sum(G[a] & G[b])) % 2:
                    ok = False; break
            if not ok:
                break
        if not ok:
            continue
        for a in range(R):
            for b in range(a, R):
                for c in range(b, R):
                    t = triple(G[a], G[b], G[c])
                    want = 0
                    if kind == "T" and a == b == c and a < k:
                        want = 1
                    if kind == "CCZ" and (a, b, c) == (0, 1, 2):
                        want = 1
                    if t != want:
                        ok = False; break
                if not ok:
                    break
            if not ok:
                break
        if ok and binary_rank(G) == R:
            return {"G1": G1, "G0": G0, "n": n, "k": k, "kind": kind}
    return None


def rm_puncture_candidate(m: int, punctures, r: int = 1) -> dict | None:
    """Ring version of the qmsd punctured-RM construction at p = 2: C = RM(r, m) restricted to the
    surviving columns; M_X = the subcode vanishing on the punctured columns (shortened code), logical
    X reps = the quotient; k = rank of G restricted to the punctures (must be full = |punctures|).
    Returns binary {"MX","Xlog","n","k"} ready for lift_css_from_x, or None if not full rank."""
    G = rm_binary(r, m)
    N = 2 ** m
    P = sorted(set(int(c) for c in punctures))
    keep = [c for c in range(N) if c not in P]
    GP = G[:, P] % 2
    k = binary_rank(GP)
    if k != len(P):
        return None
    # F_2 row operations: bring G to a form where the first k rows span the puncture part and the rest vanish on P
    A = np.hstack([GP, G[:, keep]]) % 2
    R = A.copy(); rows, cols = R.shape; piv_cols = []; rnk = 0
    for c in range(len(P)):
        piv = next((i for i in range(rnk, rows) if R[i, c]), None)
        if piv is None:
            continue
        R[[rnk, piv]] = R[[piv, rnk]]
        for i in range(rows):
            if i != rnk and R[i, c]:
                R[i] ^= R[rnk]
        piv_cols.append(c); rnk += 1
    Xlog = R[:rnk, len(P):]
    MX = R[rnk:, len(P):]
    MX = MX[np.any(MX, axis=1)]
    return {"MX": MX % 2, "Xlog": Xlog % 2, "n": len(keep), "k": k, "punctures": P, "m": m, "r": r}


# --------------------------------------------------------------------------- lifting
def lift01(Hbin) -> np.ndarray:
    """Canonical 0/1 lift of a binary matrix to Z_4 (a FREE module: |<H>| = 4^rank_2)."""
    H = np.atleast_2d(np.asarray(Hbin, dtype=np.int64)) % 2
    return H.copy()


def lift_with_torsion(Hbin, Tbin) -> np.ndarray:
    """[lift01(H); 2 * lift01(T)]: adds 2-torsion rows."""
    H = lift01(Hbin)
    T = 2 * lift01(Tbin)
    return np.vstack([H, T]) % 4


def zrm(r: int, m: int) -> np.ndarray:
    """Quaternary Reed-Muller ZRM(r, m) = <lift01(RM(r-1, m)), 2 RM(r, m)> (Hammons et al. 1994)."""
    low = lift01(rm_binary(r - 1, m)) if r >= 1 else np.zeros((0, 2 ** m), dtype=np.int64)
    return np.vstack([low, 2 * lift01(rm_binary(r, m))]) % 4


def binary_rank(Hbin) -> int:
    H = (np.atleast_2d(np.asarray(Hbin, dtype=np.int64)) % 2).astype(np.int64)
    H = H.copy(); rnk = 0; rows, cols = H.shape
    for c in range(cols):
        piv = next((i for i in range(rnk, rows) if H[i, c]), None)
        if piv is None:
            continue
        H[[rnk, piv]] = H[[piv, rnk]]
        for i in range(rows):
            if i != rnk and H[i, c]:
                H[i] ^= H[rnk]
        rnk += 1
    return rnk


def lift_css(MXbin, MZbin, n: int, torsion_x=None, torsion_z=None, label: str = "") -> Z4Code:
    """Lift a binary CSS code (X-stabilizers MXbin, Z-stabilizers MZbin) to Z_4.  Binary commutation
    (MX.MZ^T even) does NOT imply mod-4 commutation; where MX.MZ^T == 2 (mod 4) the pair anticommutes
    over Z_4 and the lift fails -- the caller must then use a different lift (e.g. double the offending
    Z rows: 2 z is always orthogonal to a 0/1 x when x.z is even).  This builder repairs by doubling
    offending Z rows and records which were doubled."""
    MX = lift01(MXbin) if len(MXbin) else np.zeros((0, n), dtype=np.int64)
    MZ = lift01(MZbin) if len(MZbin) else np.zeros((0, n), dtype=np.int64)
    if torsion_x is not None:
        MX = np.vstack([MX, 2 * lift01(torsion_x)]) % 4
    if torsion_z is not None:
        MZ = np.vstack([MZ, 2 * lift01(torsion_z)]) % 4
    doubled = []
    if MX.size and MZ.size:
        bad = np.nonzero(np.any((MX @ MZ.T) % 4 != 0, axis=0))[0]
        for j in bad:
            MZ[j] = (2 * MZ[j]) % 4
            doubled.append(int(j))
        assert np.all((MX @ MZ.T) % 4 == 0)
    code = Z4Code.from_css(MX, MZ, n, 4, label or "lift")
    code.lift_info = {"doubled_z_rows": doubled}
    return code


def lift_css_from_x(MXbin, Xlog_bin, n: int, label: str = "") -> Z4Code:
    """The natural Z_4 lift of a binary CSS code given its X-stabilizers and logical-X representatives:
    M_X = lift01(MXbin), M_Z = (M_X + <lift01(Xlog)>)^perp over Z_4.  Commutation is automatic and the
    logical quotient is free iff the lifted logical X's are independent of M_X mod 2 and 2*Xlog is not in M_X."""
    MX = lift01(MXbin) if len(MXbin) else np.zeros((0, n), dtype=np.int64)
    XL = lift01(Xlog_bin) if len(Xlog_bin) else np.zeros((0, n), dtype=np.int64)
    span = np.vstack([MX, XL]) if (MX.size or XL.size) else np.zeros((1, n), dtype=np.int64)
    MZ = right_kernel_howell(span, 4)
    MZ = np.array([z for z in MZ if np.any(np.array(z) % 4)], dtype=np.int64).reshape(-1, n)
    return Z4Code.from_css(MX, MZ, n, 4, label or "lift_from_x")


# --------------------------------------------------------------------------- refinements
def z_refinements(code: Z4Code, rng, n_samples: int = 50, max_added: int = 3) -> list[Z4Code]:
    """Codes obtained by adding 2*v rows (v in M_X^perp mod 2, i.e. any v with 2 v . M_X == 0) to M_Z.
    Each refinement shrinks the support and lowers k; duplicates (same Howell form) are dropped."""
    n, d = code.n, code.d
    MX = np.array(code.x_gens, dtype=np.int64).reshape(-1, n)
    MZ = np.array(code.z_rows, dtype=np.int64).reshape(-1, n)
    # candidates: 2 v with v in (MX mod 2)^perp over F_2  <=>  v . MX even  <=>  2v . MX == 0 mod 4
    K2 = right_kernel_howell((2 * MX) % 4, 4) if MX.size else [list(r) for r in np.eye(n, dtype=int)]
    cands = [np.array(v, dtype=np.int64) % 4 for v in K2 if np.any(np.array(v) % 4)]
    out, seen = [], set()
    for _ in range(n_samples):
        k_add = rng.randrange(1, max_added + 1)
        extra = []
        for _ in range(k_add):
            # random combination of the kernel generators, doubled
            coeffs = [rng.randrange(4) for _ in cands]
            v = sum(c * cv for c, cv in zip(coeffs, cands)) % 4 if cands else np.zeros(n, dtype=np.int64)
            v2 = (2 * v) % 4
            if np.any(v2):
                extra.append(v2)
        if not extra:
            continue
        newMZ = np.vstack([MZ, np.array(extra)]) if MZ.size else np.array(extra)
        if not np.all((MX @ newMZ.T) % 4 == 0):
            continue
        H = howell_form(np.vstack([np.hstack([MX, np.zeros_like(MX)]), np.hstack([np.zeros_like(newMZ), newMZ])]), 4)
        key = H.tobytes()
        if key in seen:
            continue
        seen.add(key)
        out.append(Z4Code.from_css(MX, newMZ, n, d, label=f"{code.label}+{len(extra)}Z2"))
    return out


# --------------------------------------------------------------------------- residue bookkeeping
def residue_code(G4) -> np.ndarray:
    """Rows mod 2 (the residue code generator)."""
    return (np.atleast_2d(np.asarray(G4, dtype=np.int64)) % 2).astype(np.int64)


def torsion_code(G4, n: int) -> np.ndarray:
    """Generators of {v mod 2 : 2 v in <G4>} (the torsion code)."""
    H = howell_form(np.atleast_2d(np.asarray(G4, dtype=np.int64)) % 4, 4)
    rows = []
    for r in H:
        r = np.asarray(r) % 4
        if not np.any(r):
            continue
        if np.all(r % 2 == 0):
            rows.append((r // 2) % 2)
        else:
            rows.append(r % 2)
    return np.array(rows, dtype=np.int64) if rows else np.zeros((0, n), dtype=np.int64)


def residue_distance_bound(code: Z4Code, d_max: int = 8) -> dict:
    """Binary distances of the residue CSS code (upper bounds for a FREE Z_4 code):
    d_Z(res) = min weight in ker(MX mod 2) \\ <MZ mod 2> approximated by the classical distance of
    ker(MX mod 2) (the usual bound), via qmsd.mindist.min_dependent_columns."""
    from qmsd.mindist import min_dependent_columns
    n = code.n
    out = {}
    for name, gens in (("X", code.x_gens), ("Z", code.z_rows)):
        H = residue_code(gens) if len(gens) else np.zeros((0, n), dtype=np.int64)
        H = H[np.any(H, axis=1)] if H.size else H
        if H.shape[0] == 0:
            out[f"d_{'Z' if name == 'X' else 'X'}_res_ub"] = 1
            continue
        try:
            out[f"d_{'Z' if name == 'X' else 'X'}_res_ub"] = int(min_dependent_columns(H, 2, d_max=d_max))
        except Exception:
            out[f"d_{'Z' if name == 'X' else 'X'}_res_ub"] = None
    vals = [v for v in out.values() if v is not None]
    out["d_res_ub"] = min(vals) if vals else None
    return out


# --------------------------------------------------------------------------- family catalog
def lift_families(m_values=(4, 5), rng=None, random_tri=((16, 2), (20, 4), (14, 3)), include_cube=True) -> list[Z4Code]:
    """E2 candidates (Z4Code list, label carries provenance):
      * punctured-RM(1,m) [[2^m-1, 1]] triorthogonal codes: M_X = lift01(X_stab), logical X = all-ones;
      * the cube [[8,3,2]] code (lift01 of X and Z stabilizers; logicals free);
      * random binary generalized-triorthogonal matrices (T kind, and a CCZ triple) lifted via
        lift_css_from_x (M_Z = (M_X + logicals)^perp over Z_4)."""
    rng = rng if rng is not None else np.random.default_rng(0)
    out = []
    for m in m_values:
        b = punctured_rm_triorthogonal(m)
        out.append(lift_css_from_x(b["MX"], np.ones((1, b["n"]), dtype=np.int64), b["n"], label=b["label"] + " liftX"))
    if include_cube:
        c = cube_8_3_2()
        out.append(lift_css(c["MX"], c["MZ"], c["n"], label=c["label"] + " lift01"))
    for (n, k) in random_tri:
        for kind in ("T", "CCZ") if k == 3 else ("T",):
            tri = random_binary_generalized_triorthogonal(n, k, rng, n_stab=2, kind=kind)
            if tri is None:
                continue
            out.append(lift_css_from_x(tri["G0"], tri["G1"], n, label=f"rand-tri {kind} [[{n},{k}]] liftX"))
    return out


def _binary_dual(G) -> np.ndarray:
    """F_2 dual of the row space of G (via the Z_4 kernel of 2G reduced mod 2)."""
    G = np.atleast_2d(np.asarray(G, dtype=np.int64)) % 2
    n = G.shape[1]
    # x with G x == 0 mod 2  <=>  (2G) (x) == 0 mod 4 restricted to 0/1 x; use F_2 Gaussian elimination
    H = G.copy(); rows, cols = H.shape; piv_cols = []; rnk = 0
    for c in range(cols):
        piv = next((i for i in range(rnk, rows) if H[i, c]), None)
        if piv is None:
            continue
        H[[rnk, piv]] = H[[piv, rnk]]
        for i in range(rows):
            if i != rnk and H[i, c]:
                H[i] ^= H[rnk]
        piv_cols.append(c); rnk += 1
    free = [c for c in range(cols) if c not in piv_cols]
    duals = []
    for f in free:
        v = np.zeros(cols, dtype=np.int64); v[f] = 1
        for i, pc in enumerate(piv_cols):
            v[pc] = H[i, f]
        duals.append(v)
    return np.array(duals, dtype=np.int64) if duals else np.zeros((0, cols), dtype=np.int64)
