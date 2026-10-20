"""Exact dense validation of the Clifford proxy (T -> S^sigma) for the D=3 colour-code T-check over Z_d.

Data-level model: 7 data qudits, the check  C = A^dag X^{(x)n} A,  A = prod_i G^{-tau_i}  (G = T exact / S proxy),
measured ideally (ancilla/cat errors are Clifford-only and already covered by the frame simulator).
Faults: single-qudit Paulis X^a Z^b at 4 slots per data qudit (before A, after A, before A^dag, after A^dag).
After the check: post-select on the trivial stabilizer syndrome (codespace projector), as cultivation does.
Undetected-error weight u(E) = ||P_code Pi_C E psi||^2 - |<psi|Pi_C E psi>|^2   (accepted but wrong).
Proxy (Clifford) gives u in {0, 1}; exact T gives fractional values.  Compare sums over fault sets.
"""
import itertools as it, sys, json, time
import numpy as np
from primepower_msd.cultivation import color_code_plaquettes, honeycomb_signs, escaping_carrier
from primepower_msd.phase_tables import modulus

d = int(sys.argv[1]) if len(sys.argv) > 1 else 4
H, xy = color_code_plaquettes(3); tau = honeycomb_signs(H, xy); n = H.shape[1]
MZ = (H * tau) % d
M = modulus(d); zeta = np.exp(2j * np.pi / M)
X = np.array(list(it.product(range(d), repeat=n)), dtype=np.int64)       # basis index -> x
W = d ** np.arange(n - 1, -1, -1)
idx = lambda xs: (xs % d) @ W
q = escaping_carrier(d)
tabs = {"T": np.asarray(q) % M, "S": (2 * np.arange(d) ** 2) % M}      # S = zeta_{2d}^{x^2}
code_mask = ~((X @ MZ.T) % d).any(1)                                       # Z-stab +1
shiftH = [idx(X + h) for h in H]                                           # X-stab shifts
shift1 = idx(X + 1)                                                        # X^{(x)n}

def diag_A(G, sign):                                                       # prod_i G^{sign*(-tau_i)}
    ph = np.zeros(len(X), dtype=np.int64)
    for i in range(n):
        ph += sign * (-int(tau[i])) * tabs[G][X[:, i]]
    return zeta ** (ph % M)

def P_code(v):
    v = v * code_mask
    for s in shiftH:                                                       # (1/d) sum_k X^{k h}
        acc = v.copy(); w = v
        for _ in range(d - 1):
            w = w[s]; acc = acc + w                                        # w(x) = v(x+h)  (X^{-h})
        v = acc / d
    return v

plus_L = code_mask.astype(complex); plus_L /= np.linalg.norm(plus_L)

def pauli(i, a, b):
    """returns (perm, phase) for X^a Z^b on qudit i: (Pv)(x) = omega^{b(x_i - a)} v(x - a e_i)."""
    e = np.zeros(n, dtype=np.int64); e[i] = a
    perm = idx(X - e)
    ph = np.exp(2j * np.pi * b * ((X[:, i] - a) % d) / d)
    return perm, ph

def apply(v, E):
    for (perm, ph) in E:
        v = ph * v[perm]
    return v

def run(G, faults_by_slot):
    Af, Ab = diag_A(G, +1), diag_A(G, -1)
    psi = TARGET[G]
    v = apply(psi, faults_by_slot[0]); v = Af * v; v = apply(v, faults_by_slot[1])
    acc = np.zeros_like(v); w = v
    for c in range(d):
        acc += apply(Ab * apply(w, faults_by_slot[2]), faults_by_slot[3]); w = w[idx(X - 1)]
    phi = acc / d
    pc = P_code(phi)
    acc_p = np.vdot(pc, pc).real; good = abs(np.vdot(psi, phi)) ** 2
    return acc_p, acc_p - good

def check_op(G, v):
    Af, Ab = diag_A(G, +1), diag_A(G, -1)
    return Ab * (Af * v)[idx(X - 1)]

TARGET = {}
rng = np.random.default_rng(1)
for G in "TS":
    r = P_code(rng.normal(size=len(X)) + 1j * rng.normal(size=len(X)))
    acc = r.copy(); w = r
    for _ in range(d - 1):
        w = check_op(G, w); acc = acc + w
    t = P_code(acc / d); t /= np.linalg.norm(t)                         # codespace C=1 eigenvector
    Ct = check_op(G, t)
    assert np.linalg.norm(Ct - t) < 1e-9 and np.linalg.norm(P_code(t) - t) < 1e-9, G
    TARGET[G] = t
print("targets built; hosted check verified on codespace", flush=True)

# sanity: no fault -> accept 1, undetected 0
for G in "TS":
    a, u = run(G, [[], [], [], []]); assert abs(a - 1) < 1e-9 and abs(u) < 1e-9, (G, a, u)
paulis = [(a, b) for a in range(d) for b in range(d) if (a, b) != (0, 0)]
singles = [(s, i, a, b) for s in range(4) for i in range(n) for (a, b) in paulis]
P = {f: pauli(*f[1:]) for f in singles}

def fault_set(fs):
    slots = [[], [], [], []]
    for f in fs:
        slots[f[0]].append(P[f])
    return slots

out = {"d": d, "n": n}
t0 = time.time()
for name, fsets in [("single", [(f,) for f in singles]),
                    ("pair", list(it.combinations(singles, 2)) if d <= 2 or len(sys.argv) > 3 else [])]:
    res = {}
    for G in "TS":
        U = np.array([run(G, fault_set(fs))[1] for fs in fsets]); A = np.array([run(G, fault_set(fs))[0] for fs in fsets]) if name == "single" else None
        res[G] = {"sum_undetected": float(U.sum()), "n_undetected>1e-9": int((U > 1e-9).sum()), "max_u": float(U.max()) if len(U) else None,
                  "mean_accept": float(A.mean()) if A is not None else None}
    # disagreement classes
    out[name] = {"count": len(fsets), **res}
    print(name, len(fsets), res, round(time.time() - t0), "s", flush=True)
# weight-3 faults on logical supports, any slots (the leading logical-failure class): exact vs proxy
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 0
sups = set()
for sup in it.combinations(range(n), 3):
    e = np.zeros(n, dtype=np.int64); e[list(sup)] = 1
    for Hc in (MZ, H):
        for v in it.product(range(1, d), repeat=3):
            e[list(sup)] = v
            if not ((Hc @ e) % d).any():
                sups.add(sup)
combos = [(sup, sl) for sup in sorted(sups) for sl in it.product(range(4), repeat=3)]
rng = np.random.default_rng(0)
if NS and NS < len(combos):                      # sample (support, slots); always ALL Paulis within a combo
    combos = [combos[j] for j in rng.choice(len(combos), NS, replace=False)]
fs_all = [tuple((s, i, a, b) for s, i, (a, b) in zip(sl, sup, pv)) for sup, sl in combos for pv in it.product(paulis, repeat=3)]
rows = []
for fs in fs_all:
    aT, uT = run("T", fault_set(fs)); aS, uS = run("S", fault_set(fs))
    rows.append((aT, uT, aS, uS))
R = np.array(rows)
worse = R[:, 1] > R[:, 3] + 1e-6
out["weight3_logical_supports"] = {"supports": len(sups), "count": len(R), "sampled": bool(NS),
    "combos": len(combos), "per_combo_sum_equal": bool(np.allclose(R[:, 1].reshape(len(combos), -1).sum(1), R[:, 3].reshape(len(combos), -1).sum(1), atol=1e-6)),
    "sum_u_T": float(R[:, 1].sum()), "sum_u_S": float(R[:, 3].sum()),
    "n_u_T>0": int((R[:, 1] > 1e-6).sum()), "n_u_S>0": int((R[:, 3] > 1e-6).sum()),
    "n_exact_worse": int(worse.sum()), "max_excess": float((R[:, 1] - R[:, 3]).max()),
    "sum_accept_T": float(R[:, 0].sum()), "sum_accept_S": float(R[:, 2].sum()),
    "examples_exact_worse": [ {"faults": fs_all[j], "u_T": R[j, 1], "u_S": R[j, 3]} for j in np.nonzero(worse)[0][:8]]}
print("weight3", {k: v for k, v in out["weight3_logical_supports"].items() if k != "examples_exact_worse"}, round(time.time() - t0), "s", flush=True)
json.dump(out, open(f"proxy_validate_d{d}.json", "w"), indent=1, default=str)
