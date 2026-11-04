"""Circuit-level simulation of one distillation round with sdim (a Clifford-only qudit simulator).

The round.  W = n + k qudits hold the code state sum_{f in C} |f>, where C is spanned by [G0, 0] and
[H1, I] (columns: n inputs, then k outputs).  A noisy U^mu acts on each input, the inputs are measured
in the X basis, and the outputs keep the magic state up to a known Pauli frame.

Why a Clifford simulator is enough.  U is diagonal, so it commutes with every Z fault and with the Z errors of the
noisy inputs, and it is not needed to find out what those errors do.  X faults need one rule, below.  Every X^g
(g in C) and Z^h (h in C^perp) stabilises the code state and they all commute, so ``round_circuit`` skips U, reads
all of them out at the end with noiseless ancillas, and gets the X and Z parts of each fault jointly:

  Z part z (from X^g readout):  G0 . z_N != 0: detector.  H1[i] . z_N + z_{O_i} != 0: output i has a Z error.
  X part x (from Z^h readout, ``x_checks``):
      x in C: trivial, the state is invariant under shifts by C.
      x = c + y with c in C and y supported on outputs: U on the inputs sees only c, so this is a pure X error on the
          outputs, accepted, and the output has fidelity 1/p: an observable fires, not a detector.
      otherwise: x reaches the inputs, U turns it into a random Z-type error there, and a detector fires.  The model
          rejects the block.  The true escape rate is about max(p^-wt(x on inputs), p^-checks), since the random Z can be the
          identity: 1/p for one input.  ``x_escape_estimate`` sizes what this leaves out when the fault also hits an output.
      A codeword shift cannot move an output X onto the inputs: U acts between the two, so they are not equivalent.

The ancillas are bookkeeping, not hardware: the real round has no such readout, only M_X on the inputs.  Measurement
noise therefore acts on the inputs alone.  Faults after the inputs are measured (idling on the outputs) are not in
the circuit: add them by hand.  tests/test_sdim_round.py checks each single fault against exact state vectors with U.

An encoder is a plain list of operations, so a code-specific circuit is a few lines:
    ("H", q)          prepare |+> on q (every other qudit starts in |0>)
    ("ADD", t, c, f)  |x_t> -> |x_t + f x_c>, one two-qudit gate; f != 0
    ("SWAP", a, b)    wire relabelling, noiseless and not counted as a gate
    ("MUL", q, a)     |x_q> -> |a x_q>, a relabelling, noiseless and not counted as a gate
``check_encoder`` verifies a list against C without sdim.  ``sdim`` is imported only where it is used.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Noise:
    eps: float = 0.0         # input magic states: Z^k, k uniform nonzero, total probability eps
    p_gate: float = 0.0      # two-qudit depolarizing after every ADD (N2)
    p_prep: float = 0.0      # one-qudit depolarizing after every H (N1 'd')
    p_meas: float = 0.0      # outcome error: uniform nonzero shift of each measured value


def _ints(M, p):
    """A plain int64 matrix mod p: catalog and library routines may hand back galois field arrays."""
    return np.atleast_2d(np.array(np.asarray(M).view(np.ndarray), dtype=np.int64)) % p


def _rref(M, p):
    """(row-reduced copy without zero rows, pivot columns)."""
    A = _ints(M, p)
    piv, r = [], 0
    for c in range(A.shape[1]):
        if r == A.shape[0]:
            break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        A[[r, r + nz[0]]] = A[[r + nz[0], r]]
        A[r] = A[r] * pow(int(A[r, c]), -1, p) % p
        for i in range(A.shape[0]):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        piv.append(c)
        r += 1
    return A[:r], piv


def code_space(G0, H1, p):
    """Generator of C on n + k qudits (inputs first, outputs last)."""
    G0 = _ints(G0, p)
    H1 = _ints(H1, p)
    k = H1.shape[0]
    return np.vstack([np.hstack([G0, np.zeros((G0.shape[0], k), dtype=np.int64)]),
                      np.hstack([H1, np.eye(k, dtype=np.int64)])]) % p


def dual_basis(G, p):
    """Basis of the dual of the row space of G."""
    R, piv = _rref(G, p)
    W = R.shape[1]
    free = [c for c in range(W) if c not in piv]
    B = np.zeros((len(free), W), dtype=np.int64)
    for j, c in enumerate(free):
        B[j, c] = 1
        for i, q in enumerate(piv):
            B[j, q] = (-R[i, c]) % p
    assert not ((R @ B.T) % p).any()
    return B


def systematic_encoder(G, p):
    """Works for any code: |+> on a pivot per row of G, then one ADD per other nonzero entry of its row."""
    R, piv = _rref(G, p)
    ops = [("H", int(q)) for q in piv]
    for i, q in enumerate(piv):
        for j in np.nonzero(R[i])[0]:
            if j != q:
                ops.append(("ADD", int(j), int(q), int(R[i, j])))
    return ops


def check_encoder(ops, G, p):
    """Raise unless ``ops`` prepares the uniform superposition over the row space of G.  Propagates the span
    of the prepared registers through the linear gates: no sdim, works for thousands of qudits."""
    W = np.asarray(G).shape[1]
    hs = [op[1] for op in ops if op[0] == "H"]
    assert len(hs) == len(set(hs)), "H applied twice to a qudit"
    rows = {q: i for i, q in enumerate(hs)}
    state = np.zeros((W, len(hs)), dtype=np.int64)      # value of each qudit as a linear form in the registers
    prepared = set()
    for op in ops:
        if op[0] == "H":
            assert not state[op[1]].any(), f"H on qudit {op[1]} after a gate touched it"
            state[op[1], rows[op[1]]] = 1
            prepared.add(op[1])
        elif op[0] == "ADD":
            _, t, c, f = op
            assert t != c and f % p, f"bad ADD {op}"
            state[t] = (state[t] + f * state[c]) % p
        elif op[0] == "SWAP":
            state[[op[1], op[2]]] = state[[op[2], op[1]]]
        elif op[0] == "MUL":
            state[op[1]] = (state[op[1]] * op[2]) % p
        else:
            raise ValueError(op)
    A, _ = _rref(state.T, p)
    B, _ = _rref(G, p)
    assert A.shape == B.shape and (A == B).all(), "encoder prepares a different space"
    return {"two_qudit_gates": sum(op[0] == "ADD" for op in ops), "registers": len(hs)}


def _add(c, t, q, f, p):
    """|x_t> -> |x_t + f x_q>.  MUL around a CNOT is a relabelling, so this is one two-qudit gate."""
    f %= p
    if f != 1:
        c.add_gate("MUL", q, a=f)
    c.add_gate("CNOT", q, t)
    if f != 1:
        c.add_gate("MUL", q, a=pow(f, -1, p))


def encoder_circuit(ops, W, p, noise: Noise):
    """sdim Circuit for an encoder on W qudits (the caller's ancillas included)."""
    from sdim import Circuit

    c = Circuit(W, p)
    for op in ops:
        if op[0] == "H":
            c.add_gate("H", op[1])
            if noise.p_prep:
                c.add_gate("N1", [op[1]], noise_channel="d", prob=noise.p_prep)
        elif op[0] == "SWAP":
            c.add_gate("SWAP", op[1], op[2])
        elif op[0] == "MUL":
            c.add_gate("MUL", op[1], a=int(op[2]))
        else:
            _, t, q, f = op
            _add(c, t, q, f, p)
            if noise.p_gate:
                c.add_gate("N2", q, t, prob=noise.p_gate)
    return c


def x_checks(G0, H1, p):
    """Rows on n + k qudits for reading the X part: (detectors, observables).
    Detectors are the words of C^perp that vanish on the outputs.  The k observables have the identity on the outputs."""
    G0, H1 = _ints(G0, p), _ints(H1, p)
    k, n = H1.shape
    B = dual_basis(code_space(G0, H1, p), p)
    order = list(range(n, n + k)) + list(range(n))      # eliminate on the output columns first
    R, piv = _rref(B[:, order], p)
    assert piv[:k] == list(range(k)), "a word of C is supported on the outputs alone"
    return R[:, np.argsort(order)][k:], R[:, np.argsort(order)][:k]


def _read_x_type(c, row, anc, p):
    """Noiseless readout of X^row with ancilla anc."""
    c.add_gate("H", anc)
    for q in np.nonzero(row)[0]:
        _add(c, int(q), anc, int(row[q]), p)
    c.add_gate("M_X", anc)


def _read_z_type(c, row, anc, p):
    """Noiseless readout of Z^row with ancilla anc."""
    for q in np.nonzero(row)[0]:
        _add(c, anc, int(q), int(row[q]), p)
    c.add_gate("M", anc)


def round_circuit(ops, G0, H1, p, noise: Noise):
    """One round as an sdim circuit.  Measurement records, in order: X-type then Z-type readouts.  Detectors: all must be 0
    to accept.  Observables: k for Z errors on the outputs, then k for X errors."""
    G0, H1 = _ints(G0, p), _ints(H1, p)
    k, n = H1.shape
    W = n + k
    det, obs = x_checks(G0, H1, p)
    C = code_space(G0, H1, p)
    rows_x, rows_z = list(C), list(det) + list(obs)
    total = len(rows_x) + len(rows_z)
    c = encoder_circuit(ops, W + total, p, noise)
    if noise.eps:
        c.add_gate("N1", list(range(n)), noise_channel="p", prob=noise.eps)
    if noise.p_meas:
        c.add_gate("N1", list(range(n)), noise_channel="p", prob=noise.p_meas)
    for j, row in enumerate(rows_x):
        _read_x_type(c, row, W + j, p)
    for j, row in enumerate(rows_z):
        _read_z_type(c, row, W + len(rows_x) + j, p)
    rec = lambda j: f"rec[{j - total}]"
    for j in range(G0.shape[0]):
        c.add_gate("DETECTOR", expr=rec(j))
    for j in range(len(det)):
        c.add_gate("DETECTOR", expr=rec(len(rows_x) + j))
    for i in range(k):
        c.add_gate("LOGICAL_OBSERVABLE", expr=rec(G0.shape[0] + i))
    for i in range(k):
        c.add_gate("LOGICAL_OBSERVABLE", expr=rec(len(rows_x) + len(det) + i))
    return c


def sample_round(circuit, shots, p, k, seed=0, batch=100_000):
    """Counts from the detector error model of ``round_circuit``.  An output with an X error alone is wrong with
    probability 1 - 1/p (overlap 1/p); ``eps_out`` includes that factor."""
    from sdim.dem import DetectorErrorModel

    sampler = DetectorErrorModel.from_circuit(circuit).compile_sampler(seed=seed)
    acc = fz = fx = blocks = done = 0
    while done < shots:
        b = min(batch, shots - done)
        D, O = sampler.sample(b)
        ok = ~np.any(D, axis=1)
        z, x = O[ok, :k] != 0, O[ok, k:] != 0
        acc += int(ok.sum())
        fz += int(z.sum())
        fx += int((x & ~z).sum())
        blocks += int((z | x).any(axis=1).sum())
        done += b
    eps = (fz + (1 - 1 / p) * fx) / (k * acc) if acc else float("nan")
    return {"shots": shots, "accepted": acc, "P_acc": acc / shots, "fail_Z": fz, "fail_X_only": fx, "failed_blocks": blocks, "eps_out": eps}


def propagate(ops, loc, x, z, p, W):
    """Final (x, z) images on W qudits of a Pauli inserted after operation index ``loc``; x, z dicts q -> power."""
    X = np.zeros(W, dtype=np.int64)
    Z = np.zeros(W, dtype=np.int64)
    for q, a in x.items():
        X[q] = a % p
    for q, b in z.items():
        Z[q] = b % p
    for op in ops[loc + 1:]:
        if op[0] == "ADD":
            _, t, c, f = op
            X[t] = (X[t] + f * X[c]) % p
            Z[c] = (Z[c] - f * Z[t]) % p
        elif op[0] == "SWAP":
            X[[op[1], op[2]]] = X[[op[2], op[1]]]
            Z[[op[1], op[2]]] = Z[[op[2], op[1]]]
        elif op[0] == "MUL":
            X[op[1]] = (X[op[1]] * op[2]) % p
            Z[op[1]] = (Z[op[1]] * pow(int(op[2]), -1, p)) % p
    return X, Z


def classify(X, Z, G0, H1, p):
    """Two-sector prediction for one fault: ('rejected', None) or ('accepted', Z failures, X failures)."""
    G0 = _ints(G0, p)
    H1 = _ints(H1, p)
    n = H1.shape[1]
    det, obs = x_checks(G0, H1, p)
    if ((det @ X) % p).any() or ((G0 @ Z[:n]) % p).any():
        return "rejected", None, None
    zf = ((H1 @ Z[:n]) + Z[n:]) % p != 0
    xf = (obs @ X) % p != 0
    return "accepted", int(zf.sum()), int(xf.sum())


def x_escape_estimate(X, G0, H1, p):
    """First-order size of what the model drops for an X image it rejects: P(accept and wrong output) if the image also
    has an output part, taking the random Z on the inputs as independent uniform values."""
    G0 = _ints(G0, p)
    H1 = _ints(H1, p)
    n = H1.shape[1]
    det, obs = x_checks(G0, H1, p)
    if not ((det @ X) % p).any() or not X[n:].any():
        return 0.0
    w = int((X[:n] != 0).sum())
    return max(float(p) ** -w, float(p) ** -G0.shape[0]) * (1 - 1 / p)
