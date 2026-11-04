"""The two-sector Clifford model of a one-shot round (qmsd/sdim_round.py) against exact state vectors with U included."""
import itertools

import numpy as np
import pytest

from qmsd.sdim_round import check_encoder, classify, code_space, propagate, systematic_encoder

P = 7
W = np.exp(2j * np.pi / P)


def rs_code():
    """RM(1,1) on F_7, output at x = 0: [[6,1,2]]_7.  Inputs are x = 1..6, the output is last."""
    xs = np.arange(1, 7)
    return np.array([xs]) % P, np.array([np.ones(6, dtype=np.int64)]) % P


def add_gate(psi, t, c, f):
    """|x_t> -> |x_t + f x_c> on a tensor with one axis per qudit."""
    out = np.empty_like(psi)
    idx = np.indices(psi.shape, sparse=False)
    src = list(idx)
    src[t] = (idx[t] - f * idx[c]) % P
    out = psi[tuple(src)]
    return out


def pauli(psi, q, a, b):
    """Apply Z^b then X^a on qudit q."""
    phase = W ** (b * np.arange(P) % P)
    shape = [1] * psi.ndim
    shape[q] = P
    return np.roll(psi * phase.reshape(shape), a, axis=q)


def run_state(ops, Wn, fault=None):
    """Prepare the code state with ``ops``; fault = (loc, {q: (a, b)}) is applied after ops[loc]."""
    psi = np.zeros((P,) * Wn, dtype=complex)
    psi[(0,) * Wn] = 1
    for i, op in enumerate(ops):
        if op[0] == "H":
            psi = np.fft.ifft(psi, axis=op[1], norm="ortho")
        else:
            psi = add_gate(psi, op[1], op[2], op[3])
        if fault and fault[0] == i:
            for q, (a, b) in fault[1].items():
                psi = pauli(psi, q, a, b)
    return psi


def measure_inputs(psi, n, G0):
    """Apply U to the inputs and go to the X basis; returns the accepted outcomes' output vectors."""
    cube = np.array([pow(x, 3, P) for x in range(P)])
    for q in range(n):
        shape = [1] * psi.ndim
        shape[q] = P
        psi = psi * (W ** cube).reshape(shape)
        psi = np.fft.fft(psi, axis=q, norm="ortho")
    out = {}
    for b in itertools.product(range(P), repeat=n):
        if ((G0 @ np.array(b)) % P).any():
            continue
        out[b] = psi[b]
    return out


def exact_fault(ops, G0, n, reference, fault):
    """(P_accept, P_accept_and_wrong) of the output against the fault-free output for the same outcome."""
    got = measure_inputs(run_state(ops, n + 1, fault), n, G0)
    acc = wrong = 0.0
    for b, v in got.items():
        pr = float((abs(v) ** 2).sum())
        if pr < 1e-14:
            continue
        ref = reference[b]
        fid = abs(np.vdot(ref / np.linalg.norm(ref), v)) ** 2 / pr
        acc += pr
        wrong += pr * (1 - fid)
    return acc, wrong


@pytest.fixture(scope="module")
def setup():
    G0, H1 = rs_code()
    C = code_space(G0, H1, P)
    ops = systematic_encoder(C, P)
    check_encoder(ops, C, P)
    n = G0.shape[1]
    ref = measure_inputs(run_state(ops, n + 1), n, G0)
    clean = sum(float((abs(v) ** 2).sum()) for v in ref.values())
    return G0, H1, ops, n, ref, clean


def test_clean_round_is_deterministic(setup):
    G0, H1, ops, n, ref, clean = setup
    assert abs(clean - 1) < 1e-10           # every outcome passes the check
    cube = np.array([W ** pow(x, 3, P) for x in range(P)])
    for v in ref.values():                  # and leaves a magic state up to a Z power and inverse
        v = v / np.linalg.norm(v)
        best = max(abs(np.vdot(np.conj(cube) * W ** (k * np.arange(P) % P) / np.sqrt(P), v)) ** 2 for k in range(P))
        assert best > 1 - 1e-10


def expected(kind, zf, xf):
    """Output infidelity of an accepted fault: Z powers are orthogonal, an X power leaves overlap 1/p."""
    return 1.0 if zf else (1 - 1 / P if xf else 0.0)


def faults(ops):
    for loc, op in enumerate(ops):
        if op[0] == "SWAP":
            continue
        yield loc, ([op[1]] if op[0] == "H" else [op[1], op[2]])


def test_z_faults_match_the_model_exactly(setup):
    G0, H1, ops, n, ref, _ = setup
    for loc, qs in faults(ops):
        for q in qs:
            for b in (1, 3):
                X, Z = propagate(ops, loc, {}, {q: b}, P, n + 1)
                kind, zf, xf = classify(X, Z, G0, H1, P)
                acc, wrong = exact_fault(ops, G0, n, ref, (loc, {q: (0, b)}))
                if kind == "rejected":
                    assert acc < 1e-10
                else:
                    assert abs(acc - 1) < 1e-10 and abs(wrong - expected(kind, zf, xf)) < 1e-10


def test_x_faults_match_the_model(setup):
    """Accepted-by-model faults are exact; rejected ones escape only at the chance rate p^-checks (1/7 here)."""
    G0, H1, ops, n, ref, _ = setup
    kinds = {"trivial": 0, "output": 0, "rejected": 0}
    for loc, qs in faults(ops):
        for q in qs:
            X, Z = propagate(ops, loc, {q: 1}, {}, P, n + 1)
            kind, zf, xf = classify(X, Z, G0, H1, P)
            acc, wrong = exact_fault(ops, G0, n, ref, (loc, {q: (1, 0)}))
            if kind == "rejected":
                kinds["rejected"] += 1
                assert abs(acc - 1 / P) < 1e-10
            else:
                kinds["output" if xf else "trivial"] += 1
                assert abs(acc - 1) < 1e-10 and abs(wrong - expected(kind, zf, xf)) < 1e-10
    assert all(kinds.values()), kinds


def test_mixed_faults_match_the_model(setup):
    """Both parts on one qudit, and on the two qudits of a gate."""
    G0, H1, ops, n, ref, _ = setup
    for loc, qs in faults(ops):
        for a, b in ((1, 1), (2, 5)):
            X, Z = propagate(ops, loc, {qs[0]: a}, {qs[-1]: b}, P, n + 1)
            kind, zf, xf = classify(X, Z, G0, H1, P)
            fault = {qs[0]: (a, 0)}
            fault[qs[-1]] = (a if len(qs) == 1 else 0, b)
            acc, wrong = exact_fault(ops, G0, n, ref, (loc, fault))
            if kind == "accepted":
                assert abs(acc - 1) < 1e-10 and abs(wrong - expected(kind, zf, xf)) < 1e-10
            else:
                assert acc <= 1 / P + 1e-10


def test_sdim_noiseless_round_is_deterministic(setup):
    pytest.importorskip("sdim")
    from qmsd.sdim_round import Noise, round_circuit, sample_round

    G0, H1, ops, *_ = setup
    r = sample_round(round_circuit(ops, G0, H1, P, Noise()), 1000, P, 1)
    assert r["accepted"] == 1000 and r["fail_Z"] == 0 and r["fail_X_only"] == 0


def test_sdim_input_noise_matches_the_exact_count(setup):
    """[[6,1,2]]_7: acceptance and output error from direct enumeration of the input errors."""
    pytest.importorskip("sdim")
    from qmsd.sdim_round import Noise, round_circuit, sample_round

    G0, H1, ops, n, *_ = setup
    eps = 0.05
    weights = np.array([1 - eps] + [eps / (P - 1)] * (P - 1))
    pa = pf = 0.0
    for e in itertools.product(range(P), repeat=n):
        e = np.array(e)
        if ((G0 @ e) % P).any():
            continue
        q = float(np.prod(weights[e]))
        pa += q
        pf += q * float(((H1 @ e) % P).any())
    r = sample_round(round_circuit(ops, G0, H1, P, Noise(eps=eps)), 400_000, P, 1, seed=3)
    assert abs(r["P_acc"] - pa) < 5 * np.sqrt(pa * (1 - pa) / 400_000)
    assert abs(r["eps_out"] - pf / pa) < 5 * np.sqrt(pf / pa / r["accepted"])
