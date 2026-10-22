"""Explicit Clifford circuit of the HHPW check-measurement protocol for Pauli-frame simulation (zdsim).

T gates are replaced by identity: twirled T faults are Z^a, which commute with T, so the Z-fault statistics are
exact; X-type (Clifford) errors are treated conservatively (the real T-check would also detect part of them).
Inputs are |+> (the identity-proxy "magic" state; the transversal check X^1 stabilises the encoded inputs).

Per check (in schedule order):
  encoder  : CSUM network realising x -> x L, with L an invertible n x n matrix whose pivot rows are the
             X-stabilisers S (pivot qudits prepared in |+>) and the normal logical X basis N (pivot qudits = the
             check's k input qudits); remaining rows are unit vectors (ancillas in |0>);
  T layer 1: ZERR(p_T) on all n data qudits;
  check    : control |+>, CSUM(control -> q) for every data qudit, control measured in the X basis (detector);
  T layer 2: ZERR(p_T) on all data;
  decoder  : inverse CSUM network; S-pivots measured in the X basis (detect Z), ancillas in Z (detect X);
             the input qudits now carry their (possibly corrupted) states out of the block.
Clifford noise p_C: DEP2 after every CSUM, DEP1 after every F/FD, XERR after resets and before measurements.
"""
from __future__ import annotations

import numpy as np

from .zdsim import Circuit, FrameSim
from .meas_sim import Protocol


def _inv_unit(a, d):
    return pow(int(a) % d, -1, d)


def completion(G, d):
    """Pivot columns P (one per row of G, unit pivot) and an invertible L (n x n, det 1) with L[P[i]] = G[i]."""
    G = np.asarray(G, dtype=np.int64) % d
    m, n = G.shape
    A = G.copy()
    piv = []
    for i in range(m):
        c = next(c for c in range(n) if c not in piv and A[i, c] % 2 == 1)
        piv.append(c)
        inv = _inv_unit(A[i, c], d)
        for i2 in range(i + 1, m):
            A[i2] = (A[i2] - A[i2, c] * inv * A[i]) % d
    L = np.zeros((n, n), dtype=np.int64)
    rest = [c for c in range(n) if c not in piv]
    for i, c in enumerate(piv):
        L[c] = G[i]
    for t in rest:
        L[t, t] = 1
    det = int(round(_det_mod(L, d)))
    if det != 1:                                  # rescale one ancilla row (an X on |0> is not a stabiliser)
        t = rest[0] if rest else None
        assert t is not None, "no ancilla to fix the determinant"
        L[t, t] = _inv_unit(det, d)
    assert _det_mod(L, d) == 1
    return piv, L


def _det_mod(L, d):
    A = L.copy() % d
    n = len(A)
    det = 1
    for c in range(n):
        r = next((r for r in range(c, n) if A[r, c] % 2 == 1), None)
        if r is None:
            return 0
        if r != c:
            A[[r, c]] = A[[c, r]]; det = -det
        det = det * A[c, c] % d
        inv = _inv_unit(A[c, c], d)
        for r2 in range(c + 1, n):
            A[r2] = (A[r2] - A[r2, c] * inv * A[c]) % d
    return det % d


def csum_network(L, d):
    """Gate list [(c, t, w)] (CSUM: x_t += w x_c, i.e. x -> x E with E = I + w e_c e_t^T) whose product is L
    (det L = 1).  Column reduction: row r is turned into e_r using only columns >= r for the pivot."""
    n = len(L)
    A = L.copy() % d
    ops = []

    def colop(c, t, w):                            # column t += w * column c
        w %= d
        if w:
            A[:, t] = (A[:, t] + w * A[:, c]) % d
            ops.append((c, t, w))

    for r in range(n):
        if A[r, r] % 2 == 0:
            c = next(c for c in range(r + 1, n) if A[r, c] % 2 == 1)
            colop(c, r, 1)
        if A[r, r] != 1 and r < n - 1:
            others = [c for c in range(r + 1, n) if A[r, c] % 2 == 1]
            if others:
                h = others[0]
            else:
                h = r + 1
                colop(r, h, 1)
            colop(h, r, (1 - A[r, r]) * _inv_unit(A[r, h], d))
        assert A[r, r] == 1, "pivot"
        for t in range(n):
            if t != r and A[r, t]:
                colop(r, t, -A[r, t])
    assert np.array_equal(A % d, np.eye(n, dtype=np.int64)), "network synthesis failed"
    # L E_1 ... E_m = I  =>  L = E_m^-1 ... E_1^-1 : gates applied in that order
    return [(c, t, (-w) % d) for (c, t, w) in reversed(ops)]


def _apply_network(net, n, d):
    M = np.eye(n, dtype=np.int64)
    for c, t, w in net:
        M[:, t] = (M[:, t] + w * M[:, c]) % d
    return M


def build(P: Protocol, p_in: float, p_T: float, p_C: float, encoder=None):
    """encoder: optional (piv, L, net) with piv = S-pivots then one pivot per row of P.N, L[piv] = (S basis; P.N)
    (e.g. sparse_encoder.sparse_protocol); default: dense completion of (S; N)."""
    d, n, k, r = P.d, P.n, P.k, P.S.shape[0]
    if encoder is None:
        G = np.vstack([P.S, P.N]) % d
        piv, L = completion(G, d)
        net = csum_network(L, d)
    else:
        piv, L, net = encoder
        assert np.array_equal(L[list(piv[r:])] % d, P.N % d)
    assert np.array_equal(_apply_network(net, n, d), L % d)
    s_piv, n_piv = piv[:r], piv[r:]
    anc_cols = [c for c in range(n) if c not in piv]
    C = Circuit(d)
    inputs = C.q(P.n_in)
    blk = C.q(n - k)                               # reusable block qudits for the S-pivots and ancillas
    ctrl = C.q(1)[0]
    C.R(inputs); C.F(inputs)
    if p_C: C.XERR(inputs, p_C); C.DEP1(inputs, p_C)
    zops = {"inputs": len(C.ops)}
    C.ZERR(inputs, p_in)
    other = [c for c in range(n) if c not in n_piv]  # block columns hosted on `blk`
    for j in P.order:
        col2q = {}
        for sl, c in enumerate(n_piv):
            col2q[c] = inputs[P.check_inputs[j, sl]]
        for i, c in enumerate(other):
            col2q[c] = blk[i]
        data = [col2q[c] for c in range(n)]
        sq = [col2q[c] for c in s_piv]
        aq = [col2q[c] for c in anc_cols]
        C.R(sq + aq)
        if p_C: C.XERR(sq + aq, p_C)
        C.F(sq)
        if p_C: C.DEP1(sq, p_C)
        for c, t, w in net:
            C.CSUM([(col2q[c], col2q[t])], w)
            if p_C: C.DEP2([(col2q[c], col2q[t])], p_C)
        zops[(int(j), 1)] = len(C.ops)
        C.ZERR(data, p_T)
        C.R([ctrl]); C.F([ctrl])
        if p_C: C.XERR([ctrl], p_C); C.DEP1([ctrl], p_C)
        for q in data:
            C.CSUM([(ctrl, q)], 1)
            if p_C: C.DEP2([(ctrl, q)], p_C)
        C.FD([ctrl])
        if p_C: C.DEP1([ctrl], p_C); C.XERR([ctrl], p_C)
        m = C.M([ctrl]); C.detector([(m[0], 1)], f"check{j}")
        zops[(int(j), 2)] = len(C.ops)
        C.ZERR(data, p_T)
        for c, t, w in reversed(net):
            C.CSUM([(col2q[c], col2q[t])], (-w) % d)
            if p_C: C.DEP2([(col2q[c], col2q[t])], p_C)
        C.FD(sq)
        if p_C: C.DEP1(sq, p_C); C.XERR(sq + aq, p_C)
        ms = C.M(sq)
        for i in ms: C.detector([(i, 1)], f"synX{j}")
        ma = C.M(aq)
        for i in ma: C.detector([(i, 1)], f"synZ{j}")
    C.zops = zops
    return C, inputs


def forced_verdicts(P: Protocol, configs, encoder=None):
    """configs: list of fault lists [(location, value)] in meas_sim location order (inputs, then per scheduled
    check 2n T locations).  Returns a bool array: accepted with a Z error on some output (circuit, noiseless
    otherwise)."""
    C, inputs = build(P, 0.0, 0.0, 0.0, encoder)
    n, N, d = P.n, P.n_in, P.d
    S = len(configs)
    forced = {}
    for shot, cfg in enumerate(configs):
        for loc, v in cfg:
            if loc < N:
                k, q = C.zops["inputs"], loc
            else:
                jj, l = divmod(loc - N, 2 * n)
                k, q = C.zops[(int(P.order[jj]), 1 + l // n)], l % n
            forced.setdefault((k, q), []).append((shot, v))
    fd = {}
    for (k, q), lst in forced.items():
        mask = np.zeros(S, dtype=np.int64); val = np.zeros((S, 1, 2), dtype=np.int64)
        for shot, v in lst:
            mask[shot] = 1; val[shot, 0, 1] = (val[shot, 0, 1] + v) % d
        op = C.ops[k]
        fd.setdefault(k, []).append((mask, [op[1][q]], val))
    sim = FrameSim(C, S, np.random.default_rng(0), gauge=False).run(forced=fd)
    ok = ~np.any(sim.detector_values(), axis=1)
    return ok & np.any(sim.z[:, inputs] % d != 0, axis=1)


def run(P: Protocol, p_in, p_T, p_C, shots, seed=0, chunk=20000, encoder=None):
    C, inputs = build(P, p_in, p_T, p_C, encoder)
    rng = np.random.default_rng(seed)
    acc = failZ = failXZ = badZ = badXZ = 0
    done = 0
    while done < shots:
        s = min(chunk, shots - done)
        sim = FrameSim(C, s, rng, gauge=False).run()
        ok = ~np.any(sim.detector_values(), axis=1)
        zb = (sim.z[:, inputs] % P.d) != 0
        xb = (sim.x[:, inputs] % P.d) != 0
        acc += int(ok.sum())
        failZ += int((ok & zb.any(1)).sum())
        failXZ += int((ok & (zb | xb).any(1)).sum())
        badZ += int((zb[ok]).sum())
        badXZ += int(((zb | xb)[ok]).sum())
        done += s
    return {"shots": shots, "accepted": acc, "fail_shots_Z": failZ, "fail_shots_XZ": failXZ, "bad_outputs_Z": badZ, "bad_outputs_XZ": badXZ,
            "noise_locations": C.noise_locations()}
