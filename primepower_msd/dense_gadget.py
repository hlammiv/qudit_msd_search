"""Exact statevector simulation of the [[7,1,3]] HHPW check gadget (scratch-grade script, see MEAS_DISTILL 5e) (single check, 8 qudits) with real T gates:
first-order Clifford-fault coefficient = sum_f w_f P(accept | f) (1 - F_out | f)."""
import sys, itertools as it, numpy as np
from joblib import Parallel, delayed
from primepower_msd.meas_distill import signed_rm
from primepower_msd.meas_sim import Protocol
from primepower_msd.meas_circuit import build
from primepower_msd.cultivation import escaping_carrier, modulus

d = int(sys.argv[1]) if len(sys.argv) > 1 else 4
jobs = int(sys.argv[2]) if len(sys.argv) > 2 else 8
S, tau = signed_rm(1, 3); bc = np.array([[0]]); P = Protocol(S, tau, d, bc)
C, inputs = build(P, 0.0, 0.0, 1e-3)
nq, inp = C.nq, inputs[0]
Mz = modulus(d); q_tab = escaping_carrier(d)
zeta = np.exp(2j * np.pi / Mz); om = np.exp(2j * np.pi / d)
Tm = np.diag(zeta ** q_tab); Tinv = Tm.conj()
Fm = np.array([[om ** (x * y) for y in range(d)] for x in range(d)]) / np.sqrt(d); FDm = Fm.conj().T
def Xp(a): return np.roll(np.eye(d), a, axis=0)            # X^a |y> = |y+a>
def Zp(b): return np.diag(om ** (b * np.arange(d)))
t_layers = {C.zops[(0, 1)]: -1, C.zops[(0, 2)]: +1}          # op index of ZERR(data) -> T^{sign*tau}

def apply1(psi, U, q):
    psi = np.tensordot(U, psi, axes=([1], [q])); return np.moveaxis(psi, 0, q)
def csum(psi, c, t, w):
    out = np.empty_like(psi)
    for x in range(d):
        idx = [slice(None)] * nq; idx[c] = x
        out[tuple(idx)] = np.roll(psi[tuple(idx)], (w * x) % d, axis=t if t < c else t - 1)
    return out
def measure(psi, q, expect):
    """project qudit q onto |expect>; return (prob, normalized state) -- prob 0 -> (0, None)."""
    idx = [slice(None)] * nq; idx[q] = expect
    sub = psi[tuple(idx)]; p = float(np.vdot(sub, sub).real)
    if p < 1e-15: return 0.0, None
    new = np.zeros_like(psi); new[tuple(idx)] = sub / np.sqrt(p); return p, new

def run(input_state, fault=None, expected=None):
    """fault: (op index, [(qudit, a, b)]) applied after that op.  expected: dict meas_index -> outcome (None: record
    the most likely outcome -- used for the noiseless reference).  Returns (accept prob, final state, outcomes)."""
    psi = np.zeros((d,) * nq, dtype=complex); psi[(0,) * nq] = 1.0
    acc, outcomes = 1.0, {}
    for k, op in enumerate(C.ops):
        name = op[0]
        if name == "R":
            for qq in op[1]:                                   # fresh qudits only (s = 1): assert |0>
                idx = [slice(None)] * nq; idx[qq] = 0
                assert abs(np.vdot(psi[tuple(idx)], psi[tuple(idx)]).real - 1) < 1e-9, "reset of a used qudit"
        elif name == "F":
            for qq in op[1]: psi = apply1(psi, Fm, qq)
        elif name == "FD":
            for qq in op[1]: psi = apply1(psi, FDm, qq)
        elif name == "CSUM":
            for (c, t), w in zip(op[1], op[2]): psi = csum(psi, c, t, w)
        elif name == "M":
            for qq, mi in zip(op[1], op[2]):
                if expected is None:
                    probs = [float(np.vdot(psi[tuple([slice(None)] * qq + [o])], psi[tuple([slice(None)] * qq + [o])]).real) for o in range(d)]
                    o = int(np.argmax(probs))
                else:
                    o = expected[mi]
                p, psi = measure(psi, qq, o); acc *= p; outcomes[mi] = o
                if psi is None: return 0.0, None, outcomes
        elif name == "ZERR":
            if k == C.zops["inputs"]:
                psi = apply1(psi, input_state, inp)            # input |+> -> given single-qudit unitary applied
            elif k in t_layers:
                sgn = t_layers[k]
                for c, qq in enumerate(op[1]):
                    psi = apply1(psi, Tm if sgn * tau[c] > 0 else Tinv, qq)
        if fault is not None and fault[0] == k:
            for qq, a, b in fault[1]: psi = apply1(psi, Xp(a) @ Zp(b), qq)
    return acc, psi, outcomes

def reduced(psi, q):
    psi_q = np.moveaxis(psi, q, 0).reshape(d, -1); return psi_q @ psi_q.conj().T

# reference: from |+> input, postselect, read off the magic state psi_T on the output qudit
acc0, psi0, out0 = run(np.eye(d), expected=None)
rho = reduced(psi0, inp); w, v = np.linalg.eigh(rho); psiT = v[:, -1]
print(f"d={d}: reference acceptance from |+> {acc0:.4f}, output purity {w[-1]:.6f}, outcomes {sorted(out0.values())}", flush=True)
# unitary preparing psiT from |+> (so the input ZERR site applies it)
plus = np.ones(d) / np.sqrt(d); Uprep = np.eye(d, dtype=complex) - 2 * np.outer(plus, plus) if False else None
# Householder: map |+> to psiT
vv = psiT - plus * np.vdot(plus, psiT) / abs(np.vdot(plus, psiT)) if abs(np.vdot(plus, psiT)) > 1e-12 else psiT - plus
Uprep = np.eye(d) - 2 * np.outer(vv, vv.conj()) / np.vdot(vv, vv)
Uprep = Uprep * (np.vdot(plus, psiT) / abs(np.vdot(plus, psiT))) if abs(np.vdot(plus, psiT)) > 1e-12 else Uprep
assert abs(abs(np.vdot(psiT, Uprep @ plus)) - 1) < 1e-9
acc1, psi1, out1 = run(Uprep, expected=out0)
F1 = float(np.vdot(psiT, reduced(psi1, inp) @ psiT).real)
print(f"noiseless with psi_T input: acceptance {acc1:.6f}, output fidelity {F1:.6f}", flush=True)
assert abs(acc1 - 1) < 1e-6 and abs(F1 - 1) < 1e-6

faults = []
for k, op in enumerate(C.ops):
    if op[0] in ("DEP1", "XERR"):
        for qq in op[1]:
            vals = [(a, b) for a in range(d) for b in range(d) if (a, b) != (0, 0)] if op[0] == "DEP1" else [(a, 0) for a in range(1, d)]
            wgt = 1.0 / len(vals)
            for a, b in vals: faults.append((k, [(qq, a, b)], wgt, op[0]))
    elif op[0] == "DEP2":
        for (c, t) in op[1]:
            vals = [v for v in it.product(range(d), repeat=4) if any(v)]
            wgt = 1.0 / len(vals)
            for v in vals: faults.append((k, [(c, v[0], v[1]), (t, v[2], v[3])], wgt, "DEP2"))
print("faults", len(faults), flush=True)

def one(f):
    k, fl, wgt, nm = f
    acc, psi, _ = run(Uprep, fault=(k, fl), expected=out0)
    if psi is None: return (k, nm, wgt, 0.0, 1.0)
    F = float(np.vdot(psiT, reduced(psi, inp) @ psiT).real)
    return (k, nm, wgt, acc, F)
res = Parallel(n_jobs=jobs, batch_size=64)(delayed(one)(f) for f in faults)
coef = sum(w * acc * (1 - F) for k, nm, w, acc, F in res)
rej = sum(w * (1 - acc) for k, nm, w, acc, F in res)
print(f"d={d}: exact first-order per-output error coefficient a = {coef:.4f}; rejection coefficient (locations-equivalent) {rej:.2f}")
# breakdown by region
from collections import defaultdict
by = defaultdict(float)
marks = sorted((v, kk) for kk, v in C.zops.items())
def region(k):
    prev = None
    for v, key in marks:
        if v > k: break
        prev = key
    return "prep" if prev is None else "enc" if prev == "inputs" else ("ctrl" if prev[1] == 1 else "dec")
for k, nm, w, acc, F in res: by[(region(k), nm)] += w * acc * (1 - F)
for key in sorted(by):
    if by[key] > 1e-6: print(f"  {key}: {by[key]:.4f}")
