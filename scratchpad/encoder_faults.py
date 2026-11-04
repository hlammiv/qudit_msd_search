"""Single-fault analysis of the axis-by-axis one-shot encoder (CODE_METRICS.md section 6), with the
circuit run on PHYSICAL qudits.  The circuit is linear (add-multiple gates from the Gauss-Jordan
decomposition of each line map), so a fault propagates exactly: an X^a after gate g becomes the shift
a L_g e_q on the final code state, a Z^b becomes b (L_g^{-1})^T e_q, with L_g = product of the later gates.
Classification of a final Z image z (grid vector), C = RM(r,m) on the grid, O = outputs, N = the rest:
  trivial    G z = 0               (z in C^perp: acts trivially on the code state)
  detected   G0 z_N != 0           (changes the syndrome of the X measurements)
  HARMFUL    G0 z_N = 0, G z != 0  (undetected and changes the outputs: first-order failure)
For X images: a shift by a codeword is trivial, so x is equivalent to x' = x - c with c in C, c|_O = x|_O,
i.e. a shift supported on N; after the T layer it acts on those qudits like Z-type errors of random values,
which are undetected only if they form a kernel word of weight >= d.  We report the support size of x'.
usage: encoder_faults.py p m r ARTIFACT_ID OUT.json"""
import sys, json, itertools, time, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.results import load_result_catalog
from scratchpad.fast_encoder import line_matrix
p, m, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); aid = sys.argv[4]; t0 = time.time()
rec = next(x for x in load_result_catalog() if x.artifact_id == aid)
O = sorted(c - 1 for c in rec.puncture_columns); W = p ** m
P = np.asarray(points(m, p), dtype=np.int64); index = {tuple(int(t) for t in x): i for i, x in enumerate(P)}
def gate_seq(M):
    """Add-multiple gates (target, control, f) implementing x -> M x, in application order."""
    A = M.copy() % p; n = A.shape[0]; ops = []; perm = list(range(n))
    for c in range(n):
        piv = next(i for i in range(c, n) if A[i, c])
        if piv != c: A[[c, piv]] = A[[piv, c]]; ops.append(("swap", c, piv))
        if A[c, c] != 1: s = pow(int(A[c, c]), -1, p); A[c] = (A[c] * s) % p; ops.append(("scale", c, s))
        for i in range(n):
            if i != c and A[i, c]: f = int(A[i, c]); A[i] = (A[i] - f * A[c]) % p; ops.append(("add", i, c, (-f) % p))
    # ops reduce M to I: R_last ... R_1 M = I  =>  M = R_1^-1 ... R_last^-1: apply inverses in reverse order
    out = []
    for op in reversed(ops):
        if op[0] == "add": out.append(("add", op[1], op[2], (-op[3]) % p))
        elif op[0] == "scale": out.append(("scale", op[1], pow(op[2], -1, p)))
        else: out.append(op)
    return out
seqs = {s: gate_seq(line_matrix(p, s)) for s in range(1, p + 1)}
# build the full gate list on grid qudits; qudit index = grid point whose i-th coordinate is the line position
gates = []
for i in range(m):
    for other in itertools.product(range(p), repeat=m - 1):
        later = other[i:]                           # exponents still to come on axes i+1..m-1
        rem = r - sum(later); s = 0 if rem < 0 else min(p, rem + 1)
        if s == 0: continue
        wires = [index[tuple(other[:i]) + (t,) + tuple(other[i:])] for t in range(p)]
        for op in seqs[s]:
            if op[0] == "add": gates.append((wires[op[1]], wires[op[2]], op[3]))
            elif op[0] == "swap": gates.append(("swap", wires[op[1]], wires[op[2]]))
print(f"p={p} m={m} r={r}: {sum(1 for g in gates if g[0] != 'swap')} add gates ({time.time()-t0:.0f}s)", flush=True)
G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
Nmask = np.ones(W, dtype=bool); Nmask[O] = False
from qmsd.triorthogonal import shorten_matrix
G0 = np.asarray(shorten_matrix(G, O, p), dtype=np.int64) % p
# interpolating codewords for X images: H (k x W) with H|_O = I
from scratchpad.fast_matrices import _rref
R_, piv = _rref(np.ascontiguousarray(np.hstack([G[:, O], G])), p); k = len(O)
rowsO = [i for i in range(R_.shape[0]) if piv[i] < k]; A = R_[rowsO][:, k:]; E = A[:, O]
Ei = _rref(np.ascontiguousarray(np.hstack([E, np.eye(k, dtype=np.int64)])), p)[0][:, k:]; H = (Ei @ A) % p
# backward sweep: Lx = later-gate map (acts on X shifts), Lz = (L^-1)^T (acts on Z)
Lx = np.eye(W, dtype=np.int64); Lz = np.eye(W, dtype=np.int64)
counts = {"Z_trivial": 0, "Z_detected": 0, "Z_harmful": 0, "Z_harmful_direct_output": 0}
xw = []; zsupp_harm = []
outset = set(O)
def classify_z(z, q):
    if not (G @ z % p).any(): counts["Z_trivial"] += 1; return
    if (G0 @ z[Nmask] % p).any(): counts["Z_detected"] += 1; return
    counts["Z_harmful"] += 1; zsupp_harm.append(int((z != 0).sum()))
    if (z[Nmask] == 0).all() and int((z != 0).sum()) == 1: counts["Z_harmful_direct_output"] += 1
# X images (corrected 2026-11-04): x = c + y, c in C, y on outputs is accepted and wrong; test x_N against the dual of C|_N.
# Random combinations of that dual catch any other image except with probability p^-24, then each hit is confirmed.
from qmsd.sdim_round import dual_basis
DN = dual_basis(G[:, Nmask], p); rng_ = np.random.default_rng(0)
DR = (rng_.integers(0, p, (min(24, DN.shape[0]), DN.shape[0])) @ DN) % p
counts.update({"X_in_C": 0, "X_harmful": 0, "X_harmful_direct_output": 0, "X_escape_sum": 0.0})
def classify_x(x):
    xp = (x - x[O] @ H) % p; xw.append(int((xp[Nmask] != 0).sum()))
    if (DR @ x[Nmask] % p).any():
        if x[O].any(): counts["X_escape_sum"] += max(float(p) ** -int((x[Nmask] != 0).sum()), float(p) ** -G0.shape[0]) * (1 - 1 / p)
        return
    assert not (DN @ x[Nmask] % p).any()
    if not (x[O] - 0).any() and not x[Nmask].any(): counts["X_in_C"] += 1; return
    # y = residual on outputs after removing the codeword that matches x on the inputs
    Ey = np.vstack([G, x[None, :]]) % p
    if _rref(np.ascontiguousarray(Ey), p)[0].shape[0] == _rref(np.ascontiguousarray(G), p)[0].shape[0]: counts["X_in_C"] += 1; return
    counts["X_harmful"] += 1
    if (x != 0).sum() == 1 and x[O].any(): counts["X_harmful_direct_output"] += 1
nf = 0
for g in reversed(gates):
    if g[0] == "swap":
        _, a, b = g; Lx[:, [a, b]] = Lx[:, [b, a]]; Lz[:, [a, b]] = Lz[:, [b, a]]; continue
    t, c, f = g
    # faults right AFTER this gate on its two wires (images through the later gates = current Lx/Lz)
    for q in (t, c):
        classify_z(Lz[:, q].copy(), q); classify_x(Lx[:, q].copy()); nf += 1
    # move to before this gate: L <- L E_g ; E_g: x_t += f x_c  -> column c += f column t
    Lx[:, c] = (Lx[:, c] + f * Lx[:, t]) % p
    Lz[:, t] = (Lz[:, t] - f * Lz[:, c]) % p
for q in range(W):      # preparation faults (before the first gate)
    classify_z(Lz[:, q].copy(), q); classify_x(Lx[:, q].copy()); nf += 1
xw = np.array(xw); d = rec.d
res = {"code": rec.label, "p": p, "m": m, "r": r, "fault_locations": nf, "k": k, **counts,
       "harmful_per_output": (counts["Z_harmful"] + counts["X_harmful"]) / k,
       "Z_harmful_fraction": counts["Z_harmful"] / nf, "Z_harmful_from_spreading": counts["Z_harmful"] - counts["Z_harmful_direct_output"],
       "Z_harmful_support_median": float(np.median(zsupp_harm)) if zsupp_harm else None,
       "X_image_N_support_median": float(np.median(xw)), "X_images_with_N_support_below_d": int((xw < d).sum()), "d": d, "seconds": round(time.time() - t0)}
json.dump(res, open(sys.argv[5], "w"), indent=1); print(json.dumps(res), flush=True)
