"""Does eps depend on the species at the gate level?  p = 7.
 1. loss before vs after an ideal diagonal gate.
 2. pulse-level SNAP on Fock 0..6 of a dispersively coupled cavity: selective pi pulse, then pi pulse with
    tone phases set by the target, Gaussian envelopes, phases calibrated per target by iterating on the
    coherent block unitaries; then Lindblad with cavity loss, transmon T1, Tphi.
 Units: microseconds, rad/us.  Writes runs/weighted/species_snap.json.
"""
import json, math
import numpy as np
from scipy.linalg import expm

P = 7
W = np.exp(2j * np.pi / P)


def target_phases(a):
    return np.array([2 * np.pi * (a * pow(n, 3, P) % P) / P for n in range(P)])


def plus():
    return np.ones(P, dtype=complex) / math.sqrt(P)


def eps_of(rho_cav, a):
    psi = plus() * np.exp(1j * target_phases(a))
    return float(1 - np.vdot(psi, rho_cav @ psi).real)


def loss_channel(g):
    A = np.diag(np.sqrt(np.arange(1, P)), 1).astype(complex)
    I = np.eye(P)
    n = A.conj().T @ A
    L = np.kron(A.conj(), A) - 0.5 * np.kron(I, n) - 0.5 * np.kron(n.T, I)
    E = expm(g * L)
    return lambda r: (E @ r.flatten(order="F")).reshape(P, P, order="F")


def order_test(g=1e-3):
    ch = loss_channel(g)
    out = {}
    for a in range(1, P):
        D = np.diag(np.exp(1j * target_phases(a)))
        r0 = np.outer(plus(), plus().conj())
        before = D @ ch(r0) @ D.conj().T
        r1 = D @ r0 @ D.conj().T
        after = ch(r1)
        half = loss_channel(g / 2)
        split = half(D @ half(r0) @ D.conj().T)
        out[a] = {"loss_before_gate": eps_of(before, a), "loss_after_gate": eps_of(after, a), "half_and_half": eps_of(split, a)}
    return out


def drive(t, A, sigma, t1, t2, phi1, phi2, chi):
    """Complex amplitude multiplying |e><g| for block n (returned per tone, sum over tones)."""
    env1 = A * math.exp(-((t - t1) ** 2) / (2 * sigma ** 2))
    env2 = A * math.exp(-((t - t2) ** 2) / (2 * sigma ** 2))
    tones = np.exp(1j * chi * np.arange(P) * t)
    return env1 * np.sum(np.exp(1j * phi1) * tones) + env2 * np.sum(np.exp(1j * phi2) * tones)


def coherent_blocks(phi2, sigma, chi, dt):
    """<g|V|g> for each Fock block under the two-pulse sequence."""
    T = 12 * sigma
    t1, t2 = 3 * sigma, 9 * sigma
    A = math.pi / (sigma * math.sqrt(2 * math.pi))
    phi1 = np.zeros(P)
    steps = int(T / dt)
    psi = np.zeros((P, 2), dtype=complex)
    psi[:, 0] = 1
    ns = np.arange(P)
    h = T / steps

    def deriv(t, s):
        e = drive(t, A, sigma, t1, t2, phi1, phi2, chi) / 2
        out = np.empty_like(s)
        out[:, 0] = -1j * np.conj(e) * s[:, 1]
        out[:, 1] = -1j * (e * s[:, 0] - chi * ns * s[:, 1])
        return out

    t = 0.0
    for _ in range(steps):
        k1 = deriv(t, psi); k2 = deriv(t + h / 2, psi + h / 2 * k1)
        k3 = deriv(t + h / 2, psi + h / 2 * k2); k4 = deriv(t + h, psi + h * k3)
        psi = psi + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
    return psi[:, 0]


def calibrate(a, sigma, chi, dt, iters=4):
    theta = target_phases(a)
    # second pulse with relative phase pi + x gives phase -/+ x on g; find the sign once
    phi = np.pi + theta
    g = coherent_blocks(phi, sigma, chi, dt)
    err_p = np.abs(np.angle(g * np.exp(-1j * theta))).sum()
    phi_m = np.pi - theta
    g_m = coherent_blocks(phi_m, sigma, chi, dt)
    err_m = np.abs(np.angle(g_m * np.exp(-1j * theta))).sum()
    s = 1 if err_p < err_m else -1
    phi = np.pi + s * theta
    for _ in range(iters):
        g = coherent_blocks(phi, sigma, chi, dt)
        phi = phi - s * np.angle(g * np.exp(-1j * theta))
    return phi, coherent_blocks(phi, sigma, chi, dt)


def coherent_eps(g, a):
    psi = plus() * np.exp(1j * target_phases(a))
    return float(1 - abs(np.vdot(psi, plus() * g)) ** 2)


def lindblad(a, phi2, sigma, chi, dt, T1c, T1q, Tphi):
    T = 12 * sigma
    t1, t2 = 3 * sigma, 9 * sigma
    A = math.pi / (sigma * math.sqrt(2 * math.pi))
    phi1 = np.zeros(P)
    steps = int(T / dt)
    h = T / steps
    d = 2 * P
    ac = np.kron(np.diag(np.sqrt(np.arange(1, P)), 1), np.eye(2)).astype(complex)
    sm = np.kron(np.eye(P), np.array([[0, 1], [0, 0]])).astype(complex)         # index 0 = g, 1 = e
    sz = np.kron(np.eye(P), np.diag([0.0, 1.0])).astype(complex)
    n = ac.conj().T @ ac
    H0 = -chi * n @ sz
    eg = sm.conj().T
    Ls = [math.sqrt(1 / T1c) * ac, math.sqrt(1 / T1q) * sm, math.sqrt(2 / Tphi) * sz]
    LdL = [L.conj().T @ L for L in Ls]
    psi0 = np.kron(plus(), np.array([1, 0], dtype=complex))
    rho = np.outer(psi0, psi0.conj())

    def rhs(t, r):
        e = drive(t, A, sigma, t1, t2, phi1, phi2, chi) / 2
        H = H0 + e * eg + np.conj(e) * sm
        out = -1j * (H @ r - r @ H)
        for L, M in zip(Ls, LdL):
            out = out + L @ r @ L.conj().T - 0.5 * (M @ r + r @ M)
        return out

    t = 0.0
    for _ in range(steps):
        k1 = rhs(t, rho); k2 = rhs(t + h / 2, rho + h / 2 * k1)
        k3 = rhs(t + h / 2, rho + h / 2 * k2); k4 = rhs(t + h, rho + h * k3)
        rho = rho + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
    tgt = np.kron(plus() * np.exp(1j * target_phases(a)), np.array([1, 0], dtype=complex))
    fid = float(np.vdot(tgt, rho @ tgt).real)
    return 1 - fid, float(np.trace(rho).real)


def main():
    r = {"order": order_test()}
    for a, v in r["order"].items():
        print(a, {k: round(x, 7) for k, x in v.items()}, flush=True)
    chi = 2 * math.pi * 1.0
    dt = 0.004
    noise = {"T1c": 1000.0, "T1q": 100.0, "Tphi": 100.0}
    r["snap"] = {}
    for x in (6.0, 10.0, 16.0):
        sigma = x / chi
        rows = {}
        for a in range(1, P):
            phi, g = calibrate(a, sigma, chi, dt)
            ec = coherent_eps(g, a)
            eps, tr = lindblad(a, phi, sigma, chi, dt, **noise)
            pop = float(1 - np.min(np.abs(g)) ** 2)
            res = float(np.max(np.abs(np.angle(g * np.exp(-1j * target_phases(a))))))
            rows[a] = {"eps_coherent": ec, "eps_with_decoherence": eps, "trace": tr, "max_block_population_error": pop, "max_phase_residual": res}
            print(f"sigma chi = {x} a = {a}: coherent {ec:.3e} (pop err {pop:.2e}, phase residual {res:.2e}), with decoherence {eps:.5f}", flush=True)
        ee = [v["eps_with_decoherence"] for v in rows.values()]
        cc = [v["eps_coherent"] for v in rows.values()]
        r["snap"][x] = {"gate_time_us": 12 * sigma, "by_a": rows, "spread_with_decoherence": (max(ee) - min(ee)) / max(ee), "spread_coherent": (max(cc) - min(cc)) / max(cc)}
        print(f"  T = {12 * sigma:.2f} us: spread with decoherence {r['snap'][x]['spread_with_decoherence']:.3g}, coherent {r['snap'][x]['spread_coherent']:.3g}", flush=True)
    r["noise_us"] = noise
    json.dump(r, open("runs/weighted/species_snap.json", "w"), indent=1)


if __name__ == "__main__":
    main()
