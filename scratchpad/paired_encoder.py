"""Explicit encoder for the paired-cap round state  sum_{f in RM_7(6,4)|X} |f|_X>,  X = {x1 = 5, 6} x F_7^3 (686 qudits).
RM(6,4)|X = {(a, a + b) : a in RM(6,3), b in RM(5,3)}  (Plotkin; f(6,y) - f(5,y) has degree <= 5 in y).  Two circuits:
  A (value domain):  fast axis-wise encoder of RM(6,3) on slice 5, of RM(5,3) on slice 6, then CSUM slice5 -> slice6 at
                     every point (343 add gates).
  B (coeff domain):  |+> on the degree<=6 coefficient registers of slice 5 and degree<=5 of slice 6, CSUM the 84 degree<=6
                     coefficient registers 5 -> 6, then the RM(6,3) fast encoder on both slices.
Each circuit is SIMULATED as a sequence of invertible line maps / CSUMs acting on the support basis (state |c> -> |Mc>),
the image of the initial support is compared with rowspace(GX) from weighted_support (rank test), and the two-qudit
add-multiple gates are counted exactly as in fast_encoder.py (Gauss-Jordan of each line matrix; scalings/swaps free).
usage: PYTHONPATH=.:scratchpad python scratchpad/paired_encoder.py OUT.json"""
import sys, json, itertools, numpy as np
from qmsd.reedmuller import points
from qmsd.weighted import weighted_support, product_of_linear
from qmsd.structured_ad import _rank
from fast_encoder import line_matrix, add_gates, stage_gates
p = 7; GL = {s: (add_gates(line_matrix(p, s), p) if s else 0) for s in range(p + 1)}
LM = {s: line_matrix(p, s) for s in range(1, p + 1)}

def fast3(A, r):
    """Apply the axis-wise encoder of RM(r,3) to a batch A of shape (B,7,7,7) (coefficient -> value); returns (A, gates)."""
    g = 0
    for ax in range(3):
        later = list(range(ax + 1, 3))
        for idx in itertools.product(range(p), repeat=2):          # the two other coordinates of the line
            others = [a for a in range(3) if a != ax]; fix = dict(zip(others, idx))
            rem = r - sum(fix[a] for a in later); s = 0 if rem < 0 else min(p, rem + 1)
            if s == 0: continue
            sl = [slice(None)] + [fix.get(a, slice(None)) for a in range(3)]
            A[tuple(sl)] = (A[tuple(sl)] @ LM[s].T) % p; g += GL[s]
    return A, g

def basis_states(r):
    """Initial support basis: one unit coefficient vector per monomial of degree <= r."""
    mons = [e for e in itertools.product(range(p), repeat=3) if sum(e) <= r]
    A = np.zeros((len(mons), p, p, p), dtype=np.int64)
    for i, e in enumerate(mons): A[(i,) + e] = 1
    return A, mons

# initial support of the two slices: span{(a-coeffs, 0)} + span{(0, b-coeffs)}
A6, m6 = basis_states(6); A5, m5 = basis_states(5)
def circuit_A():
    S5 = np.concatenate([A6, np.zeros((len(m5), p, p, p), dtype=np.int64)]); S6 = np.concatenate([np.zeros_like(A6), A5])
    S5, g1 = fast3(S5, 6); S6, g2 = fast3(S6, 5)
    S6 = (S6 + S5) % p; return S5, S6, g1 + g2 + p ** 3
def circuit_B():
    S5 = np.concatenate([A6, np.zeros((len(m5), p, p, p), dtype=np.int64)]); S6 = np.concatenate([np.zeros_like(A6), A5])
    S6 = (S6 + S5) % p; cs = len(m6)                               # CSUM on the 84 degree<=6 coefficient registers
    S5, g1 = fast3(S5, 6); S6, g2 = fast3(S6, 6)
    return S5, S6, g1 + g2 + cs

g = product_of_linear(p, 4, [((1, 0, 0, 0), c) for c in range(5)]); GX, X, mu = weighted_support(p, 4, g, 6)
P4 = np.asarray(points(4, p), dtype=np.int64)[X]
rkG = _rank(GX.copy() % p, p); out = {"rank_GX": int(rkG), "n_X": len(X)}
print(f"X = {len(X)} points, rank RM(6,4)|X = {rkG}; fast3 counts: RM(6,3) {stage_gates(7,3,6)[0]}, RM(5,3) {stage_gates(7,3,5)[0]}", flush=True)
for name, fn in (("A_value_domain", circuit_A), ("B_coeff_domain", circuit_B)):
    S5, S6, gates = fn()
    R = np.array([[ (S5 if x[0] == 5 else S6)[i, x[1], x[2], x[3]] for x in P4] for i in range(S5.shape[0])], dtype=np.int64)
    rkR = _rank(R.copy(), p); rkU = _rank(np.vstack([R, GX % p]).copy(), p)
    ok = rkR == rkG == rkU
    out[name] = {"gates_per_round": int(gates), "rank_image": int(rkR), "rank_union": int(rkU), "verified": bool(ok)}
    print(f"{name}: gates {gates}, image rank {rkR}, union rank {rkU} -> {'VERIFIED = RM(6,4)|X' if ok else 'MISMATCH'}", flush=True)
    assert ok
out["old_model_two_3D_plus_686"] = 2 * stage_gates(7, 3, 6)[0] + 2 * p ** 3; out["grid_model"] = stage_gates(7, 4, 6)[0]
out["best"] = min(out["A_value_domain"]["gates_per_round"], out["B_coeff_domain"]["gates_per_round"])
print(out, flush=True); json.dump(out, open(sys.argv[1], "w"), indent=1)
