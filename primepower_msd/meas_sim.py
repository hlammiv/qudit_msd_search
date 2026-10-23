"""Fault model, Monte Carlo and exact leading-order failure count for HHPW check-measurement
distillation over Z_d (d = 2^k) with a signed weakly self-dual inner code and an s-regular outer code.

Protocol (one round): N noisy magic-state inputs; the outer code groups them into M checks of width k.
Check j: encode its k inputs into the k logical slots of the inner [[n, k, D]] code (perfect Cliffords),
apply controlled-(T^tau X T^-tau) transversally (2n noisy T gates), measure the control (outcome must be
0), unencode and measure the X-stabilizer syndrome S (must be 0).  Post-select on every check.

Error model (twirled; Cliffords perfect):
  * input magic state  -> Z^e on the input, e uniform in {1..d-1} w.p. p_in   (C-eigenvalue shift)
  * each T gate        -> Z^a on its qudit, a uniform in {1..d-1} w.p. p_T
Z commutes with T, so a layer-1 fault Z^a on qudit q both shifts the outcome by a (it propagates through
the controlled transversal X^1 onto the control) and stays on the data; a layer-2 fault only stays on the
data.  Linear bookkeeping over Z_d:
  outcome_j = sum_slots u_slot * e_input + sum(z1)                 (u = c_j g_j^-1, a unit)
  syndrome  = S (z1 + z2)
  e_input  += <N_slot, z1 + z2>     (logical Z shift; N = normal X-logical basis)
Output failure: an accepted run with some e != 0.
"""
from __future__ import annotations

import itertools as it
import math

import numpy as np
from numba import njit

from .meas_distill import logical_x_basis, normal_basis
from .ringlinalg import right_kernel_howell


class Protocol:
    def __init__(self, S, tau, d, bit_checks, order=None, N=None):
        S = np.asarray(S, dtype=np.int64) % d
        tau = np.asarray(tau, dtype=np.int64)
        self.S, self.tau, self.d = S, tau, d
        self.n = S.shape[1]
        if N is None:                                  # default: Gram-Schmidt normal basis (dense)
            N = normal_basis(logical_x_basis(S, tau, d), tau, d)
        assert N is not None
        N = np.asarray(N, dtype=np.int64) % d
        self.N = N % d
        self.k = len(N)
        g = ((N * tau) * N).sum(1) % d
        c = (N * tau).sum(1) % d
        self.u = np.array([(int(cj) * pow(int(gj), -1, d)) % d for cj, gj in zip(c, g)], dtype=np.int64)
        bc = np.asarray(bit_checks, dtype=np.int64)
        self.bit_checks = bc
        self.n_in, self.s = bc.shape
        self.M = int(bc.max()) + 1
        rows = [[] for _ in range(self.M)]
        for i in range(self.n_in):
            for c_ in bc[i]:
                rows[c_].append(i)
        assert all(len(r) == self.k for r in rows), "outer row weight must equal inner k"
        self.check_inputs = np.array(rows, dtype=np.int64)               # slot order = sorted input index
        self.order = np.arange(self.M) if order is None else np.asarray(order, dtype=np.int64)
        self.T_per_run = 2 * self.n * self.M

    def locations(self):
        return self.n_in + self.T_per_run


# ----------------------------------------------------------------------------------------------------
# Monte Carlo (geometric skipping over fault locations; exact for the i.i.d. model)

@njit(cache=False)
def _gap(p):
    if p <= 0.0:
        return 1 << 62
    return int(math.floor(math.log(1.0 - np.random.random()) / math.log(1.0 - p))) + 1


@njit(cache=False)
def _mc(S, NL, u, check_inputs, order, d, n_in, p_in, p_T, shots, seed, out):
    np.random.seed(seed)
    n = S.shape[1]
    r = S.shape[0]
    k = NL.shape[0]
    M = order.shape[0]
    e = np.zeros(n_in, np.int64)
    z = np.zeros(n, np.int64)
    acc = 0
    fail_shots = 0
    bad_outputs = 0
    nxt_in = _gap(p_in) - 1
    nxt_T = _gap(p_T) - 1
    per_check = 2 * n
    tot_T = per_check * M
    for shot in range(shots):
        # inputs
        dirty = False
        e[:] = 0
        while nxt_in < n_in:
            e[nxt_in] = 1 + np.random.randint(d - 1)
            dirty = True
            nxt_in += _gap(p_in)
        nxt_in -= n_in
        if nxt_T >= tot_T and not dirty:
            nxt_T -= tot_T
            acc += 1
            continue
        ok = True
        for jj in range(M):
            j = order[jj]
            base = jj * per_check
            outc = 0
            for sl in range(k):
                outc += u[sl] * e[check_inputs[j, sl]]
            anyz = False
            while nxt_T < base + per_check:
                loc = nxt_T - base
                a = 1 + np.random.randint(d - 1)
                if loc < n:
                    z[loc] += a
                    outc += a
                else:
                    z[loc - n] += a
                anyz = True
                nxt_T += _gap(p_T)
            if ok and outc % d != 0:
                ok = False
            if anyz:
                if ok:
                    for rr in range(r):
                        sy = 0
                        for q in range(n):
                            sy += S[rr, q] * z[q]
                        if sy % d != 0:
                            ok = False
                            break
                if ok:
                    for sl in range(k):
                        sh = 0
                        for q in range(n):
                            sh += NL[sl, q] * z[q]
                        i = check_inputs[j, sl]
                        e[i] = (e[i] + sh) % d
                z[:] = 0
        nxt_T -= tot_T
        if ok:
            acc += 1
            nb = 0
            for i in range(n_in):
                if e[i] % d != 0:
                    nb += 1
            if nb:
                fail_shots += 1
                bad_outputs += nb
    out[0] = acc
    out[1] = fail_shots
    out[2] = bad_outputs


def monte_carlo(P: Protocol, p_in, p_T=None, shots=10**6, seed=0):
    p_T = p_in if p_T is None else p_T
    out = np.zeros(3, np.int64)
    _mc(P.S, P.N, P.u, P.check_inputs, P.order, P.d, P.n_in, float(p_in), float(p_T), int(shots), int(seed), out)
    return {"shots": int(shots), "accepted": int(out[0]), "fail_shots": int(out[1]), "bad_outputs": int(out[2])}


# ----------------------------------------------------------------------------------------------------
# Exact leading-order counting

def _closed_sets(P: Protocol, W: int):
    """All input sets v (connected through checks) with |v| + 2 * #unique-checks <= W (binary bound)."""
    bc, ci, s = P.bit_checks, P.check_inputs, P.s
    found = set()

    def rec(V, deg, frozen):
        U, opened = 0, []
        for b in V:
            for c in bc[b]:
                if deg[c] == 1:
                    U += 1
                    if c not in frozen:
                        opened.append(c)
        t = len(V)
        if not opened:
            if t + 2 * U <= W:
                found.add(frozenset(V))
            return
        if t >= W:
            return
        umax = (W - t) // 2
        if len(opened) > s * (W - t) + max(0, umax - len(frozen)):
            return
        best = min(opened)
        if len(frozen) + 1 <= umax:
            rec(V, deg, frozen | {best})
        for b in ci[best]:
            b = int(b)
            if b in V or b < V[0] or any(int(c) in frozen for c in bc[b]):
                continue
            for c in bc[b]:
                deg[c] += 1
            rec(V + [b], deg, frozen)
            for c in bc[b]:
                deg[c] -= 1

    deg = np.zeros(P.M, np.int64)
    for b0 in range(P.n_in):
        for c in bc[b0]:
            deg[c] += 1
        rec([b0], deg, frozenset())
        for c in bc[b0]:
            deg[c] -= 1
    return found


def count_input_failures(P: Protocol, W: int):
    """Weighted number of weight-W fault configurations made of input errors + outcome-flip gadgets
    (layer-1 Z^a and layer-2 Z^-a on the same inner qudit: n choices, value fixed) that are accepted and
    leave an output error.  Returns {set size t: count}."""
    d, n = P.d, P.n
    res = {}
    for v in _closed_sets(P, W):
        v = sorted(v)
        checks = sorted({int(c) for b in v for c in P.bit_checks[b]})
        coef = []
        for c in checks:
            row = list(P.check_inputs[c])
            coef.append([(int(P.u[row.index(b)]) if b in row else 0) for b in v])
        coef = np.array(coef, dtype=np.int64)
        tot = 0
        for vals in it.product(range(1, d), repeat=len(v)):
            K = int(np.count_nonzero((coef @ np.array(vals)) % d))
            if len(v) + 2 * K == W:
                tot += n ** K
        if tot:
            res[len(v)] = res.get(len(v), 0) + tot
    return res


def min_weight_logicals(P: Protocol, w: int):
    """All Z_d vectors z of weight w with S z = 0 and N z != 0 (inner-code logical Z errors), via
    binary supports in ker(S mod 2) (the 2-adic argument: odd part or 2*(binary word))."""
    S2 = P.S % 2
    n = P.n
    a = w // 2
    left = {}
    for A in it.combinations(range(n), a):
        key = S2[:, list(A)].sum(1) % 2
        left.setdefault(key.tobytes(), []).append(A)
    supports = set()
    for B in it.combinations(range(n), w - a):
        key = S2[:, list(B)].sum(1) % 2
        for A in left.get(key.tobytes(), ()):
            if not set(A) & set(B):
                supports.add(tuple(sorted(A + B)))
    out = []
    d = P.d
    for A in supports:
        sub = P.S[:, list(A)] % d
        K = np.array(right_kernel_howell(sub, d), dtype=np.int64).reshape(-1, len(A))
        elems = {tuple([0] * len(A))}
        for g in K:
            elems = {tuple((np.array(x) + m * g) % d) for x in elems for m in range(d)}
        for zA in elems:
            zA = np.array(zA)
            if np.all(zA % d):
                z = np.zeros(n, np.int64); z[list(A)] = zA
                if np.any((P.N @ z) % d):
                    out.append(z)
    return supports, out


def count_inner_failures(P: Protocol, w: int):
    """Weighted number of weight-w configurations: a weight-w inner logical Z error z in one check,
    split over the two T layers with zero outcome shift, whose induced input errors pass every later check."""
    _, zs = min_weight_logicals(P, w)
    pos = np.empty(P.M, np.int64); pos[P.order] = np.arange(P.M)
    d = P.d
    total = 0
    for z in zs:
        sh = (P.N @ z) % d
        supp = np.nonzero(z)[0]
        vals = z[supp]
        splits = sum(1 for m in range(1 << len(supp))
                     if sum(int(vals[i]) for i in range(len(supp)) if m >> i & 1) % d == 0)
        for j in range(P.M):
            e = {int(P.check_inputs[j, sl]): int(sh[sl]) for sl in range(P.k) if sh[sl]}
            later = {int(c) for i in e for c in P.bit_checks[i] if pos[c] > pos[j]}
            ok = True
            for c in later:
                row = P.check_inputs[c]
                tot = sum(int(P.u[sl]) * e.get(int(row[sl]), 0) for sl in range(P.k))
                if tot % d:
                    ok = False
                    break
            if ok:
                total += splits
    return total


def leading_order(P: Protocol, W: int):
    """p_fail ~ C_W * (p/(d-1))^W  (both input and T fault rates = p, uniform nonzero values)."""
    A = count_input_failures(P, W)
    B = count_inner_failures(P, W)
    C = sum(A.values()) + B
    return {"W": W, "input_by_size": A, "inner": B, "C": C,
            "per_output_coef": C / P.n_in / (P.d - 1) ** W}


# numba enumeration of closed input sets (same branching as outer_hyper._dfs, recording instead of stopping)

@njit
def _enum(b0, V, t, deg, frozen, nfrozen, bit_checks, check_bits, W, buf, nbuf):
    s = bit_checks.shape[1]
    best = -1
    U = 0
    O = 0
    for i in range(t):
        for j in range(s):
            c = bit_checks[V[i], j]
            if deg[c] == 1:
                U += 1
                if frozen[c] == 0:
                    O += 1
                    if best == -1 or c < best:
                        best = c
    if best == -1:
        if t >= 2 and t + 2 * U <= W and nbuf[0] < buf.shape[0]:
            for i in range(t):
                buf[nbuf[0], i] = V[i]
            buf[nbuf[0], W] = t
            nbuf[0] += 1
        return
    if t >= W:
        return
    umax = (W - t) // 2
    if O > s * (W - t) + max(0, umax - nfrozen):
        return
    if nfrozen + 1 <= umax:
        frozen[best] = 1
        _enum(b0, V, t, deg, frozen, nfrozen + 1, bit_checks, check_bits, W, buf, nbuf)
        frozen[best] = 0
    w = check_bits.shape[1]
    for kk in range(w):
        b = check_bits[best, kk]
        if b <= b0:
            continue
        inV = False
        for i in range(t):
            if V[i] == b:
                inV = True
        if inV:
            continue
        ok = True
        for j in range(s):
            if frozen[bit_checks[b, j]] == 1:
                ok = False
        if not ok:
            continue
        V[t] = b
        for j in range(s):
            deg[bit_checks[b, j]] += 1
        _enum(b0, V, t + 1, deg, frozen, nfrozen, bit_checks, check_bits, W, buf, nbuf)
        for j in range(s):
            deg[bit_checks[b, j]] -= 1


@njit
def _enum_range(bit_checks, check_bits, M, W, cap, lo, hi):
    """Closed sets whose minimum input lies in [lo, hi); every set is found exactly once from its minimum."""
    buf = -np.ones((cap, W + 1), np.int64)
    nbuf = np.zeros(1, np.int64)
    V = np.zeros(W + 1, np.int64)
    deg = np.zeros(M, np.int64)
    frozen = np.zeros(M, np.int64)
    for b0 in range(lo, hi):
        V[0] = b0
        for j in range(bit_checks.shape[1]):
            deg[bit_checks[b0, j]] += 1
        _enum(b0, V, 1, deg, frozen, 0, bit_checks, check_bits, W, buf, nbuf)
        for j in range(bit_checks.shape[1]):
            deg[bit_checks[b0, j]] -= 1
    return buf[: nbuf[0]], nbuf[0] >= cap


def _enum_all(bit_checks, check_bits, M, W, cap, lo=0, hi=None):
    # plain-Python wrapper: picklable for joblib workers (a jitted dispatcher is not)
    return _enum_range(bit_checks, check_bits, M, W, cap, lo, bit_checks.shape[0] if hi is None else hi)


def closed_sets_fast(P: Protocol, W: int, cap: int = 2_000_000, jobs: int = 1):
    """Same sets as _closed_sets (|v| >= 2 only; singletons are handled analytically).
    jobs > 1 splits the root range over processes (roots are disjoint work; late roots are cheaper)."""
    if jobs <= 1:
        parts = [_enum_all(P.bit_checks, P.check_inputs, P.M, W, cap)]
    else:
        from joblib import Parallel, delayed
        edges = np.linspace(0, P.n_in, 8 * jobs + 1).astype(int)
        parts = Parallel(n_jobs=jobs)(delayed(_enum_all)(P.bit_checks, P.check_inputs, P.M, W, cap, int(lo), int(hi))
                                      for lo, hi in zip(edges[:-1], edges[1:]))
    assert not any(ov for _, ov in parts), "enumeration buffer overflow"
    return {frozenset(int(x) for x in row[: row[W]]) for buf, _ in parts for row in buf}


def count_input_failures_fast(P: Protocol, W: int, jobs: int = 1):
    d, n = P.d, P.n
    res = {}
    if W == 1 + 2 * P.s:                                   # singletons: every check unique and nonzero
        res[1] = P.n_in * (d - 1) * n ** P.s
    for v in closed_sets_fast(P, W, jobs=jobs):
        v = sorted(v)
        checks = sorted({int(c) for b in v for c in P.bit_checks[b]})
        coef = np.array([[(int(P.u[list(P.check_inputs[c]).index(b)]) if b in P.check_inputs[c] else 0) for b in v]
                         for c in checks], dtype=np.int64)
        tot = 0
        for vals in it.product(range(1, d), repeat=len(v)):
            K = int(np.count_nonzero((coef @ np.array(vals)) % d))
            if len(v) + 2 * K == W:
                tot += n ** K
        if tot:
            res[len(v)] = res.get(len(v), 0) + tot
    return res


# Partial discard: a failed check discards only its own inputs (an output is kept iff all its checks pass).
# Per-output order drops to min(D_inner, s + 1), but the discard rate ~ s (k + 2n) p no longer grows with N.

@njit(cache=False)
def _mc_partial(S, NL, u, check_inputs, order, bit_checks, d, n_in, p_in, p_T, shots, seed, out):
    np.random.seed(seed)
    n = S.shape[1]
    r = S.shape[0]
    k = NL.shape[0]
    M = order.shape[0]
    e = np.zeros(n_in, np.int64)
    z = np.zeros(n, np.int64)
    failed = np.zeros(M, np.int64)
    kept_tot = 0
    bad_kept = 0
    shots_with_bad = 0
    nxt_in = _gap(p_in) - 1
    nxt_T = _gap(p_T) - 1
    per_check = 2 * n
    tot_T = per_check * M
    for shot in range(shots):
        dirty = False
        e[:] = 0
        while nxt_in < n_in:
            e[nxt_in] = 1 + np.random.randint(d - 1)
            dirty = True
            nxt_in += _gap(p_in)
        nxt_in -= n_in
        if nxt_T >= tot_T and not dirty:
            nxt_T -= tot_T
            kept_tot += n_in
            continue
        failed[:] = 0
        for jj in range(M):
            j = order[jj]
            base = jj * per_check
            outc = 0
            for sl in range(k):
                outc += u[sl] * e[check_inputs[j, sl]]
            anyz = False
            while nxt_T < base + per_check:
                loc = nxt_T - base
                a = 1 + np.random.randint(d - 1)
                if loc < n:
                    z[loc] += a
                    outc += a
                else:
                    z[loc - n] += a
                anyz = True
                nxt_T += _gap(p_T)
            if outc % d != 0:
                failed[j] = 1
            if anyz:
                for rr in range(r):
                    sy = 0
                    for q in range(n):
                        sy += S[rr, q] * z[q]
                    if sy % d != 0:
                        failed[j] = 1
                        break
                for sl in range(k):
                    sh = 0
                    for q in range(n):
                        sh += NL[sl, q] * z[q]
                    i = check_inputs[j, sl]
                    e[i] = (e[i] + sh) % d
                z[:] = 0
        nxt_T -= tot_T
        nb = 0
        for i in range(n_in):
            keep = True
            for c in range(bit_checks.shape[1]):
                if failed[bit_checks[i, c]]:
                    keep = False
            if keep:
                kept_tot += 1
                if e[i] % d != 0:
                    nb += 1
        bad_kept += nb
        if nb:
            shots_with_bad += 1
    out[0] = kept_tot
    out[1] = bad_kept
    out[2] = shots_with_bad


def monte_carlo_partial(P: Protocol, p_in, p_T=None, shots=10**6, seed=0):
    """Returns kept-output fraction and the error rate per kept output under partial discard."""
    p_T = p_in if p_T is None else p_T
    out = np.zeros(3, np.int64)
    _mc_partial(P.S, P.N, P.u, P.check_inputs, P.order, P.bit_checks, P.d, P.n_in, float(p_in), float(p_T),
                int(shots), int(seed), out)
    kept = int(out[0])
    return {"shots": int(shots), "kept_frac": kept / (shots * P.n_in), "bad_kept": int(out[1]),
            "err_per_kept": int(out[1]) / max(kept, 1), "shots_with_bad": int(out[2])}
