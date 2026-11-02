"""Axis-by-axis ("fast") preparation of the Reed-Muller code state for a distillation round.

The round needs the state  sum_{f in RM_p(r,m)} |f(x)>_{x in F_p^m}  (outputs = qudits at O; every other
qudit gets the noisy T and is measured in the X basis, which yields the syndrome and teleports the
result -- no per-output decoding).  Evaluation of a polynomial on the whole grid factorises:
values = (V (x) V (x) ... (x) V) coeffs,  V[j, t] = j^t mod p,  so the state is prepared by putting
|+> on the coefficient registers of monomials of degree <= r (|0> elsewhere) and applying the p x p map
V along every line of axis 1, then axis 2, ..., axis m.  In stage i a line (fixed other coordinates) has
only s possibly-nonzero inputs (exponent t with t + (exponents still to come) <= r), so it needs the map
from s inputs to p outputs; we implement it as the invertible p x p matrix [[V_ss, 0], [V_rest, I]]
(points 0..s-1 first) and count its two-qudit add-multiple gates exactly by Gauss-Jordan elimination
(swaps are wire relabels, scalings single-qudit; both free in this count, as in clifford_cost.py).
Weighted supports (p = 7 codes on X subset F_7^4): prepare on the whole grid and measure the points
outside X in the X basis; the resulting phase is a Pauli correction (no extra two-qudit gates).
usage: fast_encoder.py OUT.json"""
import sys, json, itertools, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.structured_ad import _rank

def V(p): return np.array([[pow(j, t, p) if (j or t == 0) else 0 for t in range(p)] for j in range(p)], dtype=np.int64)

def line_matrix(p, s):
    M = np.eye(p, dtype=np.int64); Vp = V(p); M[:, :s] = Vp[:, :s]; return M % p

def add_gates(M, p):
    """Two-qudit add-multiple gates needed to implement M (count of nonzero eliminations in Gauss-Jordan)."""
    A = M.copy() % p; n = A.shape[0]; g = 0
    for c in range(n):
        piv = next(i for i in range(c, n) if A[i, c])
        if piv != c: A[[c, piv]] = A[[piv, c]]
        A[c] = (A[c] * pow(int(A[c, c]), -1, p)) % p
        for i in range(n):
            if i != c and A[i, c]: A[i] = (A[i] - A[i, c] * A[c]) % p; g += 1
    assert (A == np.eye(n, dtype=np.int64)).all()
    return g

def stage_gates(p, m, r):
    gl = {s: (add_gates(line_matrix(p, s), p) if s > 0 else 0) for s in range(p + 1)}
    total = 0; per_stage = []
    for i in range(m):
        later = m - 1 - i; st = 0
        for es in itertools.product(range(p), repeat=later):
            rem = r - sum(es)
            s = 0 if rem < 0 else min(p, rem + 1)
            st += gl[s]
        st *= p ** i                       # lines are also indexed by the i already-transformed value coordinates
        per_stage.append(st); total += st
    return total, per_stage, gl

def verify(p, m, r, trials=3):
    """(a) axis-wise transform == direct evaluation; (b) the transformed degree<=r coefficient space is RM(r,m)."""
    P = np.asarray(points(m, p), dtype=np.int64); rng = np.random.default_rng(0); Vp = V(p)
    exps = list(itertools.product(range(p), repeat=m))
    for _ in range(trials):
        C = np.zeros((p,) * m, dtype=np.int64)
        for e in exps:
            if sum(e) <= r: C[e] = rng.integers(p)
        A = C.copy()
        for ax in range(m): A = np.moveaxis(np.tensordot(Vp, A, axes=([1], [ax])) % p, 0, ax)
        direct = np.zeros(len(P), dtype=np.int64)
        for e in exps:
            if C[e]: direct = (direct + C[e] * np.prod([np.power(P[:, i], e[i]) % p for i in range(m)], axis=0)) % p
        fast = np.array([A[tuple(x)] for x in P])
        assert (fast == direct).all(), "axis-wise transform disagrees with direct evaluation"
    G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
    rows = []
    for e in exps:
        if sum(e) <= r:
            C = np.zeros((p,) * m, dtype=np.int64); C[e] = 1; A = C
            for ax in range(m): A = np.moveaxis(np.tensordot(Vp, A, axes=([1], [ax])) % p, 0, ax)
            rows.append(np.array([A[tuple(x)] for x in P]))
    R = np.array(rows); rk = _rank(G.copy(), p)
    assert _rank(R.copy(), p) == rk == _rank(np.vstack([G, R]).copy(), p), "transformed space != RM(r,m)"
    return rk

if __name__ == "__main__":
    for (p, m, r) in [(3, 2, 1), (5, 2, 2), (7, 2, 3), (5, 3, 3)]:
        print(f"verified p={p} m={m} r={r}: rank {verify(p, m, r)}", flush=True)
    sysrows = {x["code"]: x for x in json.load(open("runs/coverage/clifford_cost.json"))}
    codes = [("[[1745,271,7]]_7", 7, 4, 7, 271, "weighted: grid 2401, 385 extra X measurements"),
             ("[[2093,308,5]]_7", 7, 4, 7, 308, None), ("[[100,21,4]]_11", 11, 2, 6, 21, None), ("[[13332,1309,22]]_11", 11, 4, 13, 1309, None),
             ("[[237,52,6]]_17", 17, 2, 10, 52, None), ("[[78455,5066,>=68]]_17", 17, 4, 21, 5066, None), ("[[296,65,6]]_19", 19, 2, 12, 65, None),
             ("[[123284,7037,>=57]]_19", 19, 4, 23, 7037, None), ("[[666,175,6]]_29", 29, 2, 18, 175, None), ("[[519,106,5]]_5", 5, 4, 5, 106, None), ("[[206,37,4]]_3", 3, 5, 3, 37, None)]
    out = []
    for name, p, m, r, k, note in codes:
        tot, per, gl = stage_gates(p, m, r)
        old = sysrows.get(name, {}).get("gates_per_output")
        row = {"code": name, "p": p, "m": m, "r": r, "k": k, "fast_gates_per_round": tot, "fast_gates_per_output": tot / k, "stage_gates": per,
               "systematic_gates_per_output": old, "ratio": (old / (tot / k)) if old else None, "line_gates_by_s": gl, "note": note}
        out.append(row)
        print(f"{name:26s} fast {tot:>12,} per round = {tot/k:9.1f}/output | systematic {old if old is None else round(old):>9}/output | {'' if not old else f'{old/(tot/k):6.1f}x cheaper'}", flush=True)
    json.dump(out, open(sys.argv[1], "w"), indent=1)
