"""Tests of the assumption that U, U^2, U^3 raw states at p = 7 are interchangeable at equal eps0.
  A. resource: mana of |U^a>, Wigner-positivity noise level, mana bound on conversion
  B. noise: eps and the Z-power distribution of |U^a> after a channel, Pauli and non-Pauli
  C. code: exact state-vector round of a small mixed-species weighted code vs the G0, H1 prediction
Writes runs/weighted/species_verify.json.
"""
import itertools, json, math
import numpy as np

P = 7
W = np.exp(2j * np.pi / P)
CUBE = np.array([pow(x, 3, P) for x in range(P)])
X = np.roll(np.eye(P), 1, axis=0).astype(complex)
Z = np.diag([W ** x for x in range(P)])
INV2 = pow(2, -1, P)


def mpow(M, k):
    return np.linalg.matrix_power(M, k % P)


def magic(a):
    return np.array([W ** (a * int(CUBE[x]) % P) for x in range(P)]) / math.sqrt(P)


def phase_points():
    """Gross's odd-d phase-point operators A_(q,r) = T_(q,r) A_00 T_(q,r)^dag, A_00 = parity."""
    par = np.zeros((P, P), dtype=complex)
    for x in range(P):
        par[(-x) % P, x] = 1
    out = {}
    for q in range(P):
        for r in range(P):
            T = W ** (-INV2 * q * r % P) * mpow(Z, r) @ mpow(X, q)
            out[q, r] = T @ par @ T.conj().T
    return out


PP = phase_points()


def wigner(rho):
    return np.array([[np.trace(rho @ PP[q, r]).real / P for r in range(P)] for q in range(P)])


def sum_neg(rho):
    w = wigner(rho)
    assert abs(w.sum() - 1) < 1e-10
    return float(-w[w < 0].sum())


def noisy(a, e):
    """(1 - e)|U^a><U^a| + e/(p-1) sum_{k != 0} Z^k |U^a><U^a| Z^-k"""
    psi = magic(a)
    rho = (1 - e) * np.outer(psi, psi.conj())
    for k in range(1, P):
        v = mpow(Z, k) @ psi
        rho = rho + e / (P - 1) * np.outer(v, v.conj())
    return rho


def part_a():
    rows = {}
    for a in range(1, P):
        rho = noisy(a, 0.0)
        sn = sum_neg(rho)
        lo, hi = 0.0, (P - 1) / P
        for _ in range(60):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if sum_neg(noisy(a, mid)) > 1e-13 else (lo, mid)
        rows[a] = {"sum_negativity": sn, "mana": math.log(2 * sn + 1), "eps_wigner_positive": hi}
    m1 = rows[1]["mana"]
    conv = {a: rows[a]["mana"] / m1 for a in rows}
    return {"by_a": rows, "mana_ratio_to_U": conv}


def distribution(a, channel):
    """P(k) in the basis Z^k|U^a> after channel (a function rho -> rho), ie. after twirling by X_a powers."""
    psi = magic(a)
    rho = channel(np.outer(psi, psi.conj()))
    out = []
    for k in range(P):
        v = mpow(Z, k) @ psi
        out.append(float(np.vdot(v, rho @ v).real))
    return np.array(out)


def kraus(Ks):
    return lambda r: sum(K @ r @ K.conj().T for K in Ks)


def amp_damp(g):
    """Bosonic-style loss on the ladder, truncated at d = 7: first-order Kraus with sqrt(n) rates, exact generator."""
    from scipy.linalg import expm
    a = np.diag(np.sqrt(np.arange(1, P)), 1).astype(complex)
    I = np.eye(P)
    L = np.kron(a.conj(), a) - 0.5 * np.kron(I, a.conj().T @ a) - 0.5 * np.kron((a.conj().T @ a).T, I)
    E = expm(g * L)
    return lambda r: (E @ r.flatten(order="F")).reshape(P, P, order="F")


def dephase_n(g):
    """Lindblad L = n (photon-number dephasing): diagonal, so commutes with U^a."""
    from scipy.linalg import expm
    n = np.diag(np.arange(P)).astype(complex)
    I = np.eye(P)
    L = np.kron(n.conj(), n) - 0.5 * np.kron(I, n @ n) - 0.5 * np.kron((n @ n).T, I)
    E = expm(g * L)
    return lambda r: (E @ r.flatten(order="F")).reshape(P, P, order="F")


def part_b():
    q = 1e-3
    pauli = [(i, j) for i in range(P) for j in range(P) if (i, j) != (0, 0)]
    dep = kraus([math.sqrt(1 - q) * np.eye(P)] + [math.sqrt(q / len(pauli)) * mpow(X, i) @ mpow(Z, j) for i, j in pauli])
    onlyX = kraus([math.sqrt(1 - q) * np.eye(P), math.sqrt(q) * X])
    onlyXZ = kraus([math.sqrt(1 - q) * np.eye(P), math.sqrt(q) * X @ mpow(Z, 2)])
    onlyZ = kraus([math.sqrt(1 - q) * np.eye(P), math.sqrt(q) * Z])
    chans = {"depolarizing": dep, "X only": onlyX, "X Z^2": onlyXZ, "Z only": onlyZ, "loss g=1e-3": amp_damp(1e-3), "n dephasing g=1e-3": dephase_n(1e-3)}
    # coherent over-rotation of a diagonal gate: U^a with each level phase off by a Gaussian-free fixed pattern
    theta = 1e-2
    n = np.arange(P)
    chans["coherent n^2 phase 1e-2"] = kraus([np.diag(np.exp(1j * theta * n * n / P))])
    from scipy.linalg import expm
    Vx = expm(1j * theta * (X + X.conj().T))
    chans["coherent X rot 1e-2"] = kraus([Vx])
    out = {}
    for name, ch in chans.items():
        row = {}
        for a in (1, 2, 3, 4, 5, 6):
            d = distribution(a, ch)
            row[a] = {"eps": float(1 - d[0]), "orbit_124": float(d[1] + d[2] + d[4]), "orbit_356": float(d[3] + d[5] + d[6])}
        e = [row[a]["eps"] for a in row]
        out[name] = {"by_a": row, "eps_spread_relative": (max(e) - min(e)) / max(max(e), 1e-300)}
    return out


def exact_round(G0, H1, mu, lam, eps):
    """Exact round: state sum over the row space, gate U^mu, Z-error sum, X-syndrome projection, output
    compared with |U^lam>.  Returns (P_accept, P_output_wrong_and_accepted) from the state itself."""
    n = G0.shape[1]
    r = G0.shape[0]
    words, ys = [], []
    for s in itertools.product(range(P), repeat=r):
        base = (np.array(s) @ G0) % P
        for y in range(P):
            words.append((base + y * H1[0]) % P)
            ys.append(y)
    words = np.array(words)
    ys = np.array(ys)
    phase = (CUBE[words] * mu).sum(axis=1) % P
    assert np.all(phase == (lam * CUBE[ys]) % P), "not triorthogonal with these charges"
    amp0 = W ** phase / math.sqrt(len(words))
    tgt = magic(lam)
    pa = pf = 0.0
    weights = np.array([1 - eps] + [eps / (P - 1)] * (P - 1))
    for e in itertools.product(range(P), repeat=n):
        e = np.array(e)
        pe = float(np.prod(weights[e]))
        if pe == 0:
            continue
        amp = amp0 * W ** ((words @ e) % P)
        # projecting onto trivial X-syndrome = keep the y marginal after averaging over the stabilizer sum
        out = np.array([amp[ys == y].sum() for y in range(P)]) / math.sqrt(P ** r)
        nrm = float((abs(out) ** 2).sum())
        if nrm < 1e-12:
            continue
        pa += pe * nrm
        fid = abs(np.vdot(tgt, out)) ** 2 / nrm
        pf += pe * nrm * (1 - fid)
    return pa, pf


def predicted(G0, H1, eps):
    n = G0.shape[1]
    weights = np.array([1 - eps] + [eps / (P - 1)] * (P - 1))
    pa = pf = 0.0
    for e in itertools.product(range(P), repeat=n):
        e = np.array(e)
        if np.any((G0 @ e) % P):
            continue
        pe = float(np.prod(weights[e]))
        pa += pe
        if (H1 @ e) % P:
            pf += pe
    return pa, pf


def small_weighted_code():
    """Search m = 1 charges on F_7 with deg g <= 2 (RM(1,1)) for one with all three species among inputs."""
    best = None
    for c2, c1, c0 in itertools.product(range(P), repeat=3):
        g = [(c2 * x * x + c1 * x + c0) % P for x in range(P)]
        if c2 == 0 and c1 == 0:
            continue
        S = [x for x in range(P) if g[x]]
        for out in S:
            lam = g[out]
            if lam not in (1, 6):
                continue
            keep = [x for x in S if x != out]
            kinds = {min(g[x], P - g[x]) for x in keep}
            if kinds == {1, 2, 3} and (best is None or len(keep) < len(best[1])):
                best = (g, keep, out, lam)
    g, keep, out, lam = best
    # f(x) = (x - x0) span gives G0 row (vanishes on the output); H1 = delta on output via constant + slope
    G0 = np.array([[(x - out) % P for x in keep]], dtype=np.int64)
    H1 = np.array([[1] * len(keep)], dtype=np.int64)
    mu = np.array([g[x] for x in keep], dtype=np.int64)
    return {"g": g, "keep": keep, "out": out, "lam": lam, "G0": G0, "H1": H1, "mu": mu}


def part_c():
    c = small_weighted_code()
    G0, H1, mu, lam = c["G0"], c["H1"], c["mu"], c["lam"]
    # sign of the logical: phase = sum mu (s(x-out) + y)^3 must be -lam y^3 from the point sum rule
    lam_eff = (-lam) % P
    rows = {}
    for eps in (1e-3, 1e-2, 5e-2):
        pr = predicted(G0, H1, eps)
        row = {"predicted": pr}
        for scale in (1, 2, 3):
            row[f"scale {scale}"] = exact_round(G0, H1, (scale * mu) % P, (scale * lam_eff) % P, eps)
        rows[eps] = row
    return {"n": int(G0.shape[1]), "g": c["g"], "mu": [int(x) for x in mu], "lambda": int(lam_eff), "rows": rows}


def main():
    r = {"A": part_a()}
    print("A", json.dumps(r["A"], indent=0)[:1500], flush=True)
    r["B"] = part_b()
    for k, v in r["B"].items():
        print("B", k, "spread", f"{v['eps_spread_relative']:.3g}", {a: (round(x["eps"], 8), round(x["orbit_124"], 8), round(x["orbit_356"], 8)) for a, x in v["by_a"].items() if a in (1, 2, 3)}, flush=True)
    r["C"] = part_c()
    print("C", r["C"], flush=True)
    json.dump(r, open("runs/weighted/species_verify.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
