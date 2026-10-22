"""Sparse normal logical bases and cheap CSUM encoders for the check-measurement inner codes.

The dense encoder of meas_circuit (Gram–Schmidt normal basis, generic column-reduction synthesis) spends O(k^2)
CSUMs mixing the k input qudits.  Here:

  sequential_basis : the normal basis is grown one vector at a time; each new vector is the lightest odd-norm
                     solution (random information sets, one or two free columns) of S diag(tau) x = 0 and
                     <x, N_i>_tau = 0 for the vectors chosen so far, independent mod 2 of S and of them, and zero on
                     the S pivot columns.  Valid for the given d (mod-d orthogonality).
  encoder          : pivots chosen so that L (S rows on S pivots, N_j on its pivot, unit ancilla rows) is close
                     to triangular; csum_network_ordered synthesises it in Markowitz order (little fill).

Gate counts (runs/meas_distill/sparse_encoders.json): [[127,113,3]] 10554 -> ~4k CSUMs per encoder at d = 4.
The fault model needs Protocol(..., N=basis) with the same basis (sparse_protocol returns both).
"""
from __future__ import annotations

import itertools as it

import numpy as np

from .meas_distill import _gram, logical_x_basis
from .meas_circuit import _apply_network, _det_mod, _inv_unit, csum_network


def _rank2_rows(rows):
    """F_2 rank of an iterable of python-int bit masks (returns the reduced basis dict pivot->mask)."""
    basis = {}
    for v in rows:
        v = int(v)
        while v:
            h = v.bit_length() - 1
            if h in basis:
                v ^= basis[h]
            else:
                basis[h] = v
                break
    return basis


def _insert(basis, v):
    while v:
        h = v.bit_length() - 1
        if h in basis:
            v ^= basis[h]
        else:
            basis[h] = v
            return True
    return False


def _mask(v):
    return sum(1 << int(i) for i in np.nonzero(np.asarray(v) % 2)[0])


def s_pivots(S, d, rng=None):
    """r columns where S (mod 2) has full rank (random order if rng)."""
    S2 = np.asarray(S) % 2
    cols = list(range(S2.shape[1]))
    if rng is not None:
        rng.shuffle(cols)
    piv, basis = [], {}
    for c in cols:
        if _insert(basis, _mask(S2[:, c])):
            piv.append(c)
        if len(piv) == len(_rank2_rows(_mask(r) for r in S2)):
            break
    return piv


def encoder(Sp, N, piv_S, d, S_rows=None):
    """S_rows: optional original (sparser) S rows to place on the S pivots instead of the reduced Sp.
    Pivots for N (greedy: odd entry, prefer columns where few other rows are nonzero), the matrix L and its
    CSUM network."""
    n = Sp.shape[1]
    used = set(piv_S)
    G = np.vstack([Sp, N]) % d
    # N pivots: reduce N to an echelon-like form only through choice of columns (no row ops on N)
    rem = list(range(len(N)))
    piv_N = [None] * len(N)
    W = (N % 2).copy()                                # mod-2 elimination copy: guarantees G[:, piv] invertible
    while rem:
        best = None
        for j in rem:
            for c in np.nonzero(W[j])[0]:
                if c in used:
                    continue
                cost = sum(1 for j2 in rem if N[j2, c] % d)   # 1 = private column: triangular, no mixing
                if best is None or cost < best[0]:
                    best = (cost, j, int(c))
                    if cost == 1:
                        break
            if best and best[0] == 1:
                break
        _, j, c = best
        piv_N[j] = c; used.add(c); rem.remove(j)
        for j2 in rem:
            if W[j2, c]:
                W[j2] ^= W[j]
    piv = list(piv_S) + piv_N
    Lm = np.zeros((n, n), dtype=np.int64)
    for i, c in enumerate(piv):
        Lm[c] = G[i]
    if S_rows is not None:                             # any basis of span(S) works on the S pivots
        rows = list(range(len(S_rows))); W = np.asarray(S_rows) % 2
        for c in piv_S:
            i = next(i for i in rows if W[i, c]); rows.remove(i)
            for i2 in rows:
                if W[i2, c]:
                    W[i2] ^= W[i]
            Lm[c] = np.asarray(S_rows[i]) % d
    rest = [c for c in range(n) if c not in set(piv)]
    for t in rest:
        Lm[t, t] = 1
    det = _det_mod(Lm, d)
    if det != 1:
        Lm[rest[0], rest[0]] = _inv_unit(det, d)
    return piv, Lm, csum_network_ordered(Lm, d)


def _systematic(H, cols, d):
    """Row-reduce H (mod d) choosing unit pivots in the column order `cols`; returns (Hred, pivots) or None."""
    A = H.copy() % d
    m = len(A)
    piv = []
    r = 0
    for c in cols:
        if r == m:
            break
        rows = np.nonzero(A[r:, c] % 2)[0]
        if len(rows) == 0:
            continue
        i = r + rows[0]
        if i != r:
            A[[r, i]] = A[[i, r]]
        A[r] = A[r] * _inv_unit(A[r, c], d) % d
        nz = np.nonzero(A[:, c])[0]
        for i2 in nz:
            if i2 != r:
                A[i2] = (A[i2] - A[i2, c] * A[r]) % d
        piv.append(c)
        r += 1
    if r < m:
        return None
    return A, piv


def sequential_basis(S, tau, d, rng, iters=40, avoid=(), verbose=False):
    """Normal basis built one vector at a time: the lightest odd-norm vector (found by random information sets,
    Prange style: one free column plus the pivots it forces) in {x : S diag(tau) x = 0, <x, N_i>_tau = 0}
    that is independent mod 2 of S and the vectors chosen so far.  Columns in `avoid` are kept free of support
    (they are taken as pivots first, forcing zero there), so the S-pivot rows of the encoder stay triangular."""
    S = np.asarray(S, dtype=np.int64) % d
    tau = np.asarray(tau, dtype=np.int64)
    n = S.shape[1]
    k = len(logical_x_basis(S, tau, d))
    basis = _rank2_rows(_mask(r) for r in S % 2)
    N = []
    avoid = list(avoid)
    rest = [c for c in range(n) if c not in set(avoid)]
    for t in range(k):
        H = np.vstack([S * tau] + [np.asarray(N, dtype=np.int64).reshape(-1, n) * tau]) % d
        best = None
        for _ in range(iters):
            perm = list(rng.permutation(rest))
            red = _systematic(H, perm, d)
            if red is None:
                continue
            A, piv = red
            pset = set(piv)
            for j in perm:
                if j in pset:
                    continue
                v = np.zeros(n, dtype=np.int64)
                v[j] = 1
                v[piv] = (-A[:, j]) % d
                if any(v[c] for c in avoid):
                    continue
                w = int(np.count_nonzero(v))
                if best is not None and w >= best[0]:
                    continue
                if int(((v * tau) @ v) % 2) == 0:             # need odd norm
                    continue
                if not _insert(dict(basis), _mask(v % 2)):
                    continue
                best = (w, v)
        tries2 = 0
        while best is None and tries2 < 4 * iters:      # fallback: two free columns (x_j1 = 1, x_j2 = c)
            tries2 += 1
            red = _systematic(H, list(rng.permutation(rest)), d)
            if red is None:
                continue
            A, piv = red
            free = [c for c in rest if c not in set(piv)]
            for j1, j2 in it.combinations(free, 2):
                for c in range(1, d):
                    v = np.zeros(n, dtype=np.int64)
                    v[j1] = 1; v[j2] = c
                    v[piv] = (-(A[:, j1] + c * A[:, j2])) % d
                    w = int(np.count_nonzero(v))
                    if best is not None and w >= best[0]:
                        continue
                    if int(((v * tau) @ v) % 2) == 0 or not _insert(dict(basis), _mask(v % 2)):
                        continue
                    best = (w, v)
        if best is None:
            raise RuntimeError(f"no candidate at step {t}")
        _insert(basis, _mask(best[1] % 2))
        N.append(best[1])
        if verbose and (t % 10 == 0 or t == k - 1):
            print(f"  step {t}: weight {best[0]}", flush=True)
    N = np.array(N, dtype=np.int64)
    G, Q = _gram(N, tau, d)
    assert not np.any((G - np.diag(np.diag(G))) % d) and np.all(Q % 2), "not a normal basis"
    return N


def s_reduced(S, piv_S, d):
    """S row-reduced (row operations within S) to the identity on piv_S."""
    Sp = np.asarray(S, dtype=np.int64) % d
    red = _systematic(Sp, list(piv_S), d)
    assert red is not None and red[1] == list(piv_S)
    return red[0]


def markowitz_order(L, d):
    """Processing order for csum_network minimising symbolic fill: repeatedly take the unprocessed row r (odd
    diagonal) with least (#other unprocessed rows nonzero in column r) x (#off-diagonal nonzeros in row r)."""
    B = (np.asarray(L) % d) != 0
    n = len(B)
    unp = np.ones(n, dtype=bool)
    order = []
    diag_odd = (np.diag(np.asarray(L)) % 2) == 1
    for _ in range(n):
        colc = (B & unp[:, None]).sum(0) - B.diagonal()          # other unprocessed rows hitting column r
        rowc = B.sum(1) - B.diagonal()
        cost = np.where(unp & diag_odd, colc * rowc, np.iinfo(np.int64).max)
        if np.all(cost == np.iinfo(np.int64).max):           # no odd diagonal left: fall back to index order
            order += [int(i) for i in np.nonzero(unp)[0]]
            break
        r = int(np.argmin(cost + 1e-3 * rowc))
        rows = np.nonzero(B[:, r] & unp)[0]
        rows = rows[rows != r]
        if len(rows):
            B[np.ix_(rows, B[r])] = True
        unp[r] = False
        order.append(r)
    return order


def csum_network_ordered(L, d):
    """csum_network on the symmetrically permuted matrix (Markowitz order), mapped back."""
    L = np.asarray(L, dtype=np.int64) % d
    perm = markowitz_order(L, d)
    Lp = L[np.ix_(perm, perm)]
    net = csum_network(Lp, d)
    net = [(perm[c], perm[t], w) for c, t, w in net]
    assert np.array_equal(_apply_network(net, len(L), d), L)
    return net


def sparse_protocol(S, tau, d, bit_checks, order=None, seed=0, iters=40, tries=1):
    """(Protocol with the sequential sparse normal basis, encoder (piv, L, net)) -- the cheapest of `tries`."""
    from .meas_sim import Protocol
    best = None
    rng = np.random.default_rng(seed)
    fails = 0
    while tries > 0:
        pS = s_pivots(S, d, rng)
        try:
            N = sequential_basis(S, tau, d, rng, iters=iters, avoid=pS)
        except RuntimeError:                           # greedy dead end (alternating remainder): restart
            fails += 1
            assert fails < 50, "sequential_basis keeps dead-ending"
            continue
        tries -= 1
        Sp = s_reduced(S, pS, d)
        for rows in (None, S):                         # reduced or original S rows: keep the cheaper network
            enc = encoder(Sp, N, pS, d, S_rows=rows)
            if best is None or len(enc[2]) < len(best[1][2]):
                best = (N, enc)
    N, enc = best
    return Protocol(S, tau, d, bit_checks, order=order, N=N), enc
