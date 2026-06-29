"""Low-overhead search: k-output single-ququart MSD codes, minimizing gamma = log(n/k)/log(d).

gamma < 1 is IMPOSSIBLE at k=1 (distance <= n forces gamma >= 1), so we must encode k >= 2 logical
ququarts.  For the resource to stay SINGLE-QUDIT (not an entangling logical CCZ — the field route the
project rejects), the transversal gate must induce a SEPARABLE logical gate diag = T_L^{⊗k}, i.e. the
induced phase factorizes  lambda(c_1..c_k) = sum_j lambda_j(c_j), with each lambda_j a strict level-3
single-ququart gate.

Construction (scalable, no d^n): M_X self-orthogonal under phi; greedily add logical-X generators g_j
that (a) have order d mod the current module, (b) induce a strict level-3 single-ququart gate, and
(c) are logically orthogonal (no cross phase) to the previous generators.  Distance via the bounded
ring-distance engine; gamma from (n, k, distance).
"""

from __future__ import annotations

import numpy as np

from .ringlinalg import enumerate_module, right_kernel_howell, in_span, howell_form, module_size
from .ring_distance import min_logical_weight, code_gamma
from .clifford_ring import is_level3, is_clifford


def _Phi(state, phi, N):
    return int(sum(int(phi[int(x)]) for x in state) % N)


def _induced_single(g, mx, phi, N, d, n):
    """The single-ququart gate induced by logical generator g (its phases over l=0..d-1), or None."""
    lam0 = _Phi(tuple((0 * g) % d), phi, N)
    lam = [(_Phi(tuple((ell * g) % d), phi, N) - lam0) % N for ell in range(d)]
    UL = np.diag([np.exp(2j * np.pi * l / N) for l in lam])
    return lam, UL


def _separable(g_new, g_list, mx, phi, N, d, n):
    """Logical orthogonality: Phi(c*g_new + c'*g_j) == Phi(c*g_new) + Phi(c'*g_j) for all c,c', j."""
    for gj in g_list:
        for c in range(1, d):
            cg = (c * g_new) % d
            base_c = _Phi(tuple(cg), phi, N)
            for cp in range(1, d):
                cpg = (cp * gj) % d
                lhs = _Phi(tuple((cg + cpg) % d), phi, N)
                if lhs != (base_c + _Phi(tuple(cpg), phi, N)) % N:
                    return False
    return True


def _in_V(r, mx, phi, N, d):
    """r is in the transversality region V: Phi constant on the M_X-coset of r."""
    base = _Phi(r, phi, N)
    return all(_Phi(tuple((np.array(r) + np.array(a)) % d), phi, N) == base for a in mx)


def _coset_reps(gens, n, d):
    from itertools import product as _product
    reps = []
    for co in _product(range(d), repeat=len(gens)):
        r = np.zeros(n, dtype=int)
        for cc, gg in zip(co, gens):
            r = (r + cc * np.array(gg)) % d
        reps.append(tuple(int(x) for x in r))
    return reps


def construct_code(rng, n, d, phi, N, k_target, mx_rank=(1, 3), tries_g=60, dist_cap=6):
    """Build a k-output single-ququart MSD code at block length n (incremental transversality). Returns dict|None.

    Conditions enforced: (self) Phi vanishes on M_X; each logical generator g induces a strict level-3
    single-ququart gate, has order d in the quotient, is logically orthogonal (separable) to the prior
    generators, AND keeps the whole submodule L = M_X + <g_1..g_j> inside the transversality region V
    (so the transversal gate D^{⊗n} induces the separable logical gate T_L^{⊗k}).  Distance via the
    ring-distance engine; gamma from (n, k, distance).
    """
    rx = int(rng.integers(mx_rank[0], mx_rank[1] + 1))
    A = rng.integers(0, d, size=(rx, n))
    mx = enumerate_module(A, d)
    mx_set = set(mx)
    if any(_Phi(a, phi, N) != 0 for a in mx):        # (self) Phi vanishes on M_X
        return None
    gens = []
    cur_set = set(mx)
    for _ in range(tries_g):
        if len(gens) >= k_target:
            break
        wt = int(rng.integers(2, n + 1))
        pos = rng.choice(n, size=wt, replace=False)
        g = np.zeros(n, dtype=int)
        g[pos] = rng.integers(1, d, size=wt)
        if tuple(int(x) for x in g) in cur_set:
            continue
        if next((t for t in range(1, d + 1) if tuple((t * g) % d) in cur_set), None) != d:
            continue
        lam, UL = _induced_single(g, mx, phi, N, d, n)
        if not (is_level3(UL, d) and not is_clifford(UL, d)):       # induces strict level-3
            continue
        if not _separable(g, gens, mx, phi, N, d, n):               # logically orthogonal to prior gens
            continue
        if not all(_in_V(r, mx, phi, N, d) for r in _coset_reps(gens + [g], n, d)):  # incremental transversality
            continue
        gens.append(g)
        cur_set = set(enumerate_module(np.vstack([A] + [np.array(x) for x in gens]), d))
    k = len(gens)
    if k < 2:
        return None
    Lgen = np.vstack([A] + [np.array(x) for x in gens])
    if module_size(Lgen, d) != (d ** k) * len(mx):                  # free rank-k logical quotient
        return None
    mz_gens = right_kernel_howell(Lgen, d)
    dist = min_logical_weight(howell_form(A, d), mz_gens, n, d, cap=dist_cap)
    if dist < 2:
        return None
    return {"n": n, "k": k, "d": dist, "gamma": code_gamma(n, k, dist),
            "A": A.tolist(), "gens": [g.tolist() for g in gens]}


def lg_worker(seed: int, trials: int):
    """Parallel worker: sample codes, keep the Pareto-best by (gamma, then distance)."""
    import numpy as _np
    from primepower_msd.ring_transversal import extract_phase
    from primepower_msd.hierarchy_search import search as _level3_family
    d = 4
    rng = _np.random.default_rng(seed)
    fam = _level3_family(d)
    gate_phis = [(gi, *extract_phase(U, d)) for gi, (a, b, s, U) in enumerate(fam)]
    best = None
    found = 0
    for n in (9, 11, 13, 15):
        for _ in range(trials):
            gi, phi, N = gate_phis[int(rng.integers(0, len(gate_phis)))]
            res = construct_code(rng, n, d, list(phi), N, k_target=int(rng.integers(2, 5)))
            if res is None:
                continue
            found += 1
            res["gi"] = gi
            if best is None or res["gamma"] < best["gamma"]:
                best = res
    return {"best": best, "found": found}
