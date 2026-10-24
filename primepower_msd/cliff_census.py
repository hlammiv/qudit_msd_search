"""Census of single Clifford faults: which locations give an undetected output error, by circuit position."""
import sys, itertools as it, json, numpy as np
from collections import Counter
from primepower_msd.meas_distill import signed_rm
from primepower_msd.meas_sim import Protocol
from primepower_msd.meas_circuit import build
from primepower_msd.zdsim import FrameSim
d = int(sys.argv[1]); which = sys.argv[2]
if which == "7":
    S, tau = signed_rm(1, 3); bc = np.array([[0, 1]])
elif which == "7s1":
    S, tau = signed_rm(1, 3); bc = np.array([[0]])
else:
    S, tau = signed_rm(1, 4); bc = np.array(list(it.combinations(range(8), 2)))
P = Protocol(S, tau, d, bc)
C, inputs = build(P, 0.0, 0.0, 1e-3)
inset = set(inputs)
# position markers
marks = sorted((v, k) for k, v in C.zops.items())
def region(k):
    """(check position in schedule, 'enc'|'ctrl'|'dec') for op index k."""
    prev = None
    for v, key in marks:
        if v > k: break
        prev = key
    if prev is None: return (0, "prep")
    if prev == "inputs": return (0, "enc")
    j, l = prev
    jj = list(P.order).index(j)
    if l == 1: return (jj, "ctrl")
    # after zops[(j,2)]: decoder of check jj until the next check's zops[(j',1)]
    nxt = [v for v, key in marks if v > k]
    return (jj, "dec")
# fault list: (op index, qudit(s), pauli vector)
faults = []
for k, op in enumerate(C.ops):
    name = op[0]
    if name == "DEP1" or name == "XERR":
        for q in op[1]:
            vals = [(a, b) for a in range(d) for b in range(d) if (a, b) != (0, 0)] if name == "DEP1" else [(a, 0) for a in range(1, d)]
            for v in vals: faults.append((k, (q,), (v,)))
    elif name == "DEP2":
        for (c, t) in op[1]:
            for v in it.product(range(d), repeat=4):
                if any(v): faults.append((k, (c, t), ((v[0], v[1]), (v[2], v[3]))))
S_ = len(faults); print("faults", S_, "ops", len(C.ops), flush=True)
CH = 2000
det = np.zeros(S_, bool); zout = np.zeros(S_, bool); xout = np.zeros(S_, bool); nbadZ = np.zeros(S_, int); nbad = np.zeros(S_, int)
for lo in range(0, S_, CH):
    hi = min(S_, lo + CH); forced = {}
    for s in range(lo, hi):
        k, qs, pv = faults[s]
        mask = np.zeros(hi - lo, np.int64); mask[s - lo] = 1
        arr = np.zeros((hi - lo, len(qs), 2), np.int64)
        for i, (a, b) in enumerate(pv): arr[s - lo, i] = (a, b)
        forced.setdefault(k, []).append((mask, list(qs), arr))
    sim = FrameSim(C, hi - lo, np.random.default_rng(0), gauge=False).run(scale=0.0, forced=forced)
    det[lo:hi] = np.any(sim.detector_values() != 0, axis=1)
    zout[lo:hi] = np.any(sim.z[:, inputs] % d != 0, axis=1); xout[lo:hi] = np.any(sim.x[:, inputs] % d != 0, axis=1)
    nbadZ[lo:hi] = np.count_nonzero(sim.z[:, inputs] % d, axis=1); nbad[lo:hi] = np.count_nonzero((sim.z[:, inputs] % d) | (sim.x[:, inputs] % d), axis=1)
bad = (~det) & (zout | xout)
cnt = Counter(); cntZ = Counter(); tot = Counter()
M = len(P.order)
for s, (k, qs, pv) in enumerate(faults):
    jj, reg = region(k); pos = "last" if jj == M - 1 else ("first" if jj == 0 else "mid")
    on_in = any(q in inset for q in qs)
    key = (pos, reg, C.ops[k][0], "on_input" if on_in else "block")
    tot[key] += 1
    if bad[s]: cnt[key] += 1
    if bad[s] and zout[s]: cntZ[key] += 1
print("undetected output-error single faults (count / total locations*paulis), Z-only in brackets:")
for key in sorted(tot):
    if cnt[key]: print(f"  {key}: {cnt[key]} / {tot[key]}  [Z {cntZ[key]}]")
print("total bad", int(bad.sum()), "of", S_, "Z-type", int((bad & zout).sum()))
# per-location weights: DEP1 has d^2-1 paulis each prob p/(d^2-1); DEP2 d^4-1; XERR d-1.  First-order coefficient:
coef = 0.0; coefZ = 0.0
for s, (k, qs, pv) in enumerate(faults):
    name = C.ops[k][0]; w = 1.0 / ((d * d - 1) if name == "DEP1" else (d ** 4 - 1) if name == "DEP2" else (d - 1))
    if bad[s]: coef += w * nbad[s]
    if bad[s] and zout[s]: coefZ += w * nbadZ[s]
print(f"first-order expected bad outputs per run: any {coef:.3f}, Z {coefZ:.3f}; per output: any {coef/P.n_in:.3f}, Z {coefZ/P.n_in:.3f}")
