"""Fast (G0, H1) for large punctured-RM catalog records (n >~ 5000), avoiding the generic galois path
of performance.record_matrices: numba Gaussian elimination mod p throughout.
H1 rows are the canonical logical rows (u_i G with (u_i G)|_O = e_i), restricted to the survivors."""
import numpy as np, numba as nb
from qmsd.reedmuller import rm_generator

@nb.njit(cache=True)
def _rref(A, p):
    A = A.copy(); R, C = A.shape; inv = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        for y in range(1, p):
            if (x * y) % p == 1: inv[x] = y
    piv = np.full(R, -1, dtype=np.int64); rk = 0
    for c in range(C):
        pr = -1
        for i in range(rk, R):
            if A[i, c] != 0: pr = i; break
        if pr < 0: continue
        if pr != rk:
            for j in range(C): t = A[pr, j]; A[pr, j] = A[rk, j]; A[rk, j] = t
        iv = inv[A[rk, c]]
        for j in range(C): A[rk, j] = (A[rk, j] * iv) % p
        for i in range(R):
            if i != rk and A[i, c] != 0:
                f = A[i, c]
                for j in range(C): A[i, j] = (A[i, j] - f * A[rk, j]) % p
        piv[rk] = c; rk += 1
        if rk == R: break
    return A[:rk], piv[:rk]

def _matmul_mod(A, B, p, chunk=512):
    out = np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)
    for s in range(0, A.shape[1], chunk):
        out = (out + A[:, s:s + chunk] @ B[s:s + chunk]) % p
    return out

def fast_matrices(record):
    p, m, r = record.p, record.m, record.r
    O = [c - 1 for c in record.puncture_columns]
    G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p          # K x N
    K, N = G.shape; k = len(O)
    keep = np.ones(N, dtype=bool); keep[O] = False
    # G0 = rows of rowspace(G) vanishing on O: rref of [G_O | G] on the O block first
    GO = G[:, O]
    M = np.ascontiguousarray(np.hstack([GO, G]))
    R_, piv = _rref(M, p)
    rows0 = [i for i in range(R_.shape[0]) if piv[i] >= k]               # pivots beyond the O block -> vanish on O
    G0 = R_[rows0][:, k:][:, keep]
    # H1: rows of R_ with pivot inside the O block are u_i G with (u_i G)|_O in echelon form on O;
    # since G_O has full column rank k, those k rows restricted to O form an invertible k x k matrix E
    rowsO = [i for i in range(R_.shape[0]) if piv[i] < k]
    assert len(rowsO) == k, "output columns not independent"
    A = R_[rowsO][:, k:]                                                # k x N, A|_O = E
    E = A[:, O]
    Ei = _rref(np.ascontiguousarray(np.hstack([E, np.eye(k, dtype=np.int64)])), p)[0][:, k:]
    H = _matmul_mod(Ei, A, p)
    assert (H[:, O] == np.eye(k, dtype=np.int64)).all()
    return G0, H[:, keep]
