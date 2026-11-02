"""Clifford (two-qudit gate) cost of one distillation round, for catalog codes.

Gate model (standard triorthogonal round, upper-bound counting, no circuit optimisation): after the
transversal T layer, measure each X-type check (row of G0) with one two-qudit gate per nonzero entry
and decode each output with one two-qudit gate per nonzero of its logical row.  Matrices are put in
systematic form first: G0 -> reduced row echelon (the cheapest basis we can find quickly; the true
minimum-weight basis is NP-hard), H1 -> each logical row reduced against the G0 pivots.
Reported per round: checks (rows), max / mean check weight, total gates, gates per output state
(gates / k) and per *accepted* output at the eps where the code is used (gates / (k P_acc)).
For codes too large to row-reduce here (G0 > 1.5e8 entries) the systematic weight is extrapolated
from the dense-fraction law nnz ~ r (n - r)(1 - 1/p), checked on the codes that do fit.
usage: clifford_cost.py OUT.json"""
import sys, json, math, time, numpy as np, numba as nb
from qmsd.results import load_result_catalog
from qmsd.performance import record_matrices
from qmsd.cyclic import construction_B

@nb.njit(cache=True)
def rref_mod(A, p):
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

def reduce_rows(H, G, piv, p):
    H = H.copy() % p
    for i, c in enumerate(piv):
        f = H[:, c].copy()
        H = (H - np.outer(f, G[i])) % p
    return H

IDS = ["imported-p7-n1745-k271-03f7dbe035d6", "imported-p7-n1756-k260-9d263f2240e4", "imported-p7-n1746-k270-271d5366cc68", "new-p7-m4-n2093-k308",
       "imported-p11-n9-k2-85962d023b9c", "imported-p11-n100-k21-20afd92fd535", "imported-p11-n13332-k1309-c7aae3253e15",
       "imported-p17-n13-k4-d1c98fb7704a", "new-p17-d6", "imported-p17-n78455-k5066-2d7880c18fa5",
       "imported-p19-n15-k4-f9a71e9685fd", "imported-p19-n296-k65-211e1fbb4b1a", "imported-p19-n123284-k7037-526d98fa2549",
       "imported-p29-n22-k7-cb6f781d0a7b", "imported-p29-n666-k175-4f8ad426d19b", "paper-table3-p5-m4-n519-k106", "paper-table3-p3-m5-n206-k37"]
recs = {r.artifact_id: r for r in load_result_catalog()}
rows = []
def account(name, p, n, k, d, G0, H1, eps_use):
    t0 = time.time(); G0 = np.ascontiguousarray(np.asarray(G0, dtype=np.int64) % p); H1 = np.asarray(H1, dtype=np.int64) % p
    G0 = G0[[t for t in range(G0.shape[0]) if G0[t].any()]]
    if G0.size <= 1.5e8:
        Gs, piv = rref_mod(G0, p); Hs = reduce_rows(H1, Gs, piv, p)
        w = (Gs != 0).sum(axis=1); wl = (Hs != 0).sum(axis=1); r = Gs.shape[0]; how = "rref"
        nnz_g, nnz_h = int(w.sum()), int(wl.sum()); wmax, wmean = int(w.max()), float(w.mean())
        dense_pred = r * (n - r) * (1 - 1 / p) + r
    else:
        r = G0.shape[0]; how = "extrapolated (dense fraction)"
        nnz_g = int(r * (n - r) * (1 - 1 / p) + r); nnz_h = int(k * (n - r) * (1 - 1 / p)); wmax = int((n - r) * (1 - 1 / p) + 1); wmean = float(wmax); dense_pred = nnz_g
    gates = nnz_g + nnz_h; pacc = (1 - eps_use) ** n
    row = {"code": name, "p": p, "n": n, "k": k, "d": d, "checks": int(r), "max_check_weight": wmax, "mean_check_weight": round(wmean, 1),
           "check_gates": nnz_g, "decode_gates": nnz_h, "gates_per_round": gates, "gates_per_output": round(gates / k, 1),
           "eps_use": eps_use, "P_acc": pacc, "gates_per_accepted_output": round(gates / (k * pacc), 1), "dense_fraction_vs_law": round(nnz_g / dense_pred, 3), "method": how, "seconds": round(time.time() - t0)}
    rows.append(row); print(json.dumps({kk: row[kk] for kk in ("code", "checks", "max_check_weight", "mean_check_weight", "gates_per_round", "gates_per_output", "gates_per_accepted_output", "dense_fraction_vs_law", "method")}), flush=True)
# eps at which each code is used: 1e-3 for small codes, 1e-6 for the large ones (their place in the cost-optimal protocols)
G7, H7 = construction_B(7, 8, {1, 4, 7}); account("[[7,1,3]]_7 cyc-B", 7, 7, 1, 3, G7, H7, 1e-2)
for aid in IDS:
    r = recs[aid]
    big = r.n > 10 ** 4
    if big:
        # avoid forming the 7304 x 78455 / 10k x 123k matrices: extrapolate
        if r.n > 5 * 10 ** 4:
            class Dummy: pass
            rr = r.n; k = r.k; p = r.p
            from qmsd.reedmuller import r_max  # dim A from the record geometry when available
            dimA = (r.geometry or {}).get("residual_checks") or {"[[123284,7037,>=57]]_19": 17270 - 7037}.get(r.label)
            account(r.label, p, rr, k, r.d, np.zeros((1, 1), dtype=np.int64) + 1, np.zeros((k, 1), dtype=np.int64), 1e-6) if False else None
            n = rr; rdim = dimA; nnz_g = int(rdim * (n - rdim) * (1 - 1 / p) + rdim); nnz_h = int(k * (n - rdim) * (1 - 1 / p)); g = nnz_g + nnz_h; pacc = (1 - 1e-6) ** n
            row = {"code": r.label, "p": p, "n": n, "k": k, "d": r.d, "checks": int(rdim), "max_check_weight": int((n - rdim) * (1 - 1 / p) + 1), "mean_check_weight": round((n - rdim) * (1 - 1 / p) + 1, 1),
                   "check_gates": nnz_g, "decode_gates": nnz_h, "gates_per_round": g, "gates_per_output": round(g / k, 1), "eps_use": 1e-6, "P_acc": pacc,
                   "gates_per_accepted_output": round(g / (k * pacc), 1), "dense_fraction_vs_law": 1.0, "method": "extrapolated (dense fraction)", "seconds": 0}
            rows.append(row); print(json.dumps({kk: row[kk] for kk in ("code", "checks", "max_check_weight", "gates_per_round", "gates_per_output", "gates_per_accepted_output", "method")}), flush=True)
            continue
    print(f"... {r.label}", flush=True); G0, H1 = (__import__("scratchpad.fast_matrices", fromlist=["x"]).fast_matrices(r) if r.n > 5000 and (r.construction or {}).get("family") == "punctured_rm" else record_matrices(r))
    account(r.label, r.p, r.n, r.k, r.d, G0, H1, 1e-6 if r.n > 1000 else (1e-4 if r.n > 150 else 1e-3))
json.dump(rows, open(sys.argv[1], "w"), indent=1)
