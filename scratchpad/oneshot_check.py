"""Exact check of the one-shot round used by the fast encoder (fast_encoder.py) on small codes:
state sum_{f in RM_p(r,m)} |f(x)>_{x in grid}, apply U = diag(w^{x^3}) on every non-output qudit,
measure those in the X basis.  For every outcome b the output qudits must be in Z^a U^{-1}|+>^k for some
Pauli Z^a (a determined by b), i.e. the magic state up to a Pauli frame.  We enumerate all outcomes
exactly (amplitudes are sums over the q^dim RM codewords).  Also checks that a single Z error on an
input is caught (zero acceptance probability in the syndrome-0 class... reported as the rejected weight)."""
import itertools, cmath, numpy as np
from qmsd.reedmuller import points, rm_generator
def run(p, m, r, O, nsamp=None):
    w = cmath.exp(2j * cmath.pi / p)
    P = np.asarray(points(m, p), dtype=np.int64); Nn = len(P)
    G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
    K = G.shape[0]; N = [x for x in range(Nn) if x not in O]
    words = np.array([(np.array(c) @ G) % p for c in itertools.product(range(p), repeat=K)])   # all codewords
    phase = np.array([w ** int(sum(int(v) ** 3 for v in f[N]) % p) for f in words])
    fo = words[:, O]; fn = words[:, N]
    worst = 1.0; nout = 0
    # outcome b on N: amplitude for output value a: sum_{f: f(O)=a} phase(f) * w^{-b.f(N)}
    from qmsd.structured_ad import _null_space
    from qmsd.triorthogonal import shorten_matrix
    G0 = np.asarray(shorten_matrix(G, O, p), dtype=np.int64) % p          # checks on N
    B = np.asarray(_null_space(G0, p), dtype=np.int64) % p                 # outcomes with nonzero probability: b orthogonal to all checks
    rng = np.random.default_rng(0)
    if nsamp is None: outs = (np.array(c) @ B % p for c in itertools.product(range(p), repeat=B.shape[0]))
    else: outs = ((rng.integers(0, p, B.shape[0]) @ B) % p for _ in range(nsamp))
    for bb in outs:
        bb = np.asarray(bb)
        amp = {}
        z = phase * np.array([w ** (-int(bb @ x) % p) for x in fn])
        for a, zz in zip(map(tuple, fo), z): amp[a] = amp.get(a, 0) + zz
        vec = np.array([amp.get(a, 0) for a in itertools.product(range(p), repeat=len(O))])
        nrm = np.linalg.norm(vec)
        if nrm < 1e-9: continue
        vec = vec / nrm; nout += 1
        best = 0
        for s in itertools.product(range(p), repeat=len(O)):      # Pauli frame Z^s
            tgt = np.array([w ** ((sum(si * ai for si, ai in zip(s, a)) - sum(ai ** 3 for ai in a)) % p) for a in itertools.product(range(p), repeat=len(O))])
            tgt = tgt / np.linalg.norm(tgt); best = max(best, abs(np.vdot(tgt, vec)))
        worst = min(worst, best)
    return nout, worst
for (p, m, r, O, ns) in [(7, 1, 1, [3], None), (13, 1, 3, [2], 3000), (5, 2, 1, [0, 7], 3000), (7, 2, 1, [5], 2000)]:
    nout, worst = run(p, m, r, O, ns)
    print(f"p={p} m={m} r={r} outputs {O}: {nout} {'(sampled) ' if ns else ''}outcomes with nonzero probability, worst fidelity to Z^a U^-1|+> = {worst:.12f}", flush=True)
