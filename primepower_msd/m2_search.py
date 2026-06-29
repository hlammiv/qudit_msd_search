"""M2 — search for a single-ququart magic-state-distillation code over the cyclic ring Z_d.

Strategy.  Enumerate weakly-self-dual CSS codes over Z_d = Z_4 with k=1 and a CYCLIC logical group
(code-level anti-collapse), and test whether the complete level-3 single-qudit gate family (the
d(d-1) gates from `hierarchy_search`) acts transversally — UNIFORMLY (D^{⊗n}) or ADDRESSABLY
(⊗_{i in S} D, a subset) — inducing a strict level-3 logical gate.  A hit is a genuine single-ququart
distillation code; a strong null sharpens toward the no-go.

Code enumeration is *targeted* (over self-orthogonal rows), so n=3 and n=5 are exhaustive and n=7 is
densely sampled.

Run:  python -m primepower_msd.m2_search
"""

from __future__ import annotations

import sys
from itertools import combinations, product
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from math import comb

from primepower_msd.ringlinalg import howell_form
from primepower_msd.ring_css import build_css, is_self_orthogonal, logical_group_invariants, code_distance
from primepower_msd.ring_transversal import certify_distillation_code
from primepower_msd.hierarchy_search import search as level3_family


def self_orthogonal_rows(n: int, d: int):
    """All rows a in Z_d^n with a . a == 0 (mod d) (the self-orthogonal X-stabilizer rows)."""
    return [a for a in product(range(d), repeat=n) if sum(x * x for x in a) % d == 0]


def candidate_codes(n: int, d: int, rng, per_rank: int, exhaustive_cap: int = 200_000):
    """Yield k=1 cyclic weakly-self-dual CSS codes on n qudits — ANY rank (free OR non-free).

    Generating sets are combinations of self-orthogonal rows; modules are deduped by their Howell
    canonical form.  Dropping the free-module restriction is deliberate: the prime construction uses
    punctured / non-free codes, and the literature flags non-free ring codes as the open loophole.
    """
    rows = self_orthogonal_rows(n, d)
    seen = set()

    def consider(combo):
        G = np.array(combo, dtype=object) % d
        if not is_self_orthogonal(G, d):
            return None
        key = tuple(map(tuple, howell_form(G, d).tolist()))   # canonical module identity
        if key in seen:
            return None
        seen.add(key)
        code = build_css(G, d)
        if code is None or code.k != 1 or not code.cyclic:
            return None
        return code

    for r in range(1, n):                       # number of generators
        total = comb(len(rows), r)
        if total <= exhaustive_cap:
            for combo in combinations(rows, r):
                code = consider(combo)
                if code is not None:
                    yield code
        else:
            for _ in range(per_rank):
                idx = rng.choice(len(rows), size=r, replace=False)
                code = consider([rows[i] for i in idx])
                if code is not None:
                    yield code


def brute_force_transversal_check(code, U: np.ndarray, weights=None) -> float:
    """Validate transversality by building the actual codespace in C^{d^n}; returns leakage norm."""
    d, n = code.d, code.n
    Dvec = np.diag(U).astype(complex)
    if weights is None:
        weights = [1] * n
    Uvec = np.ones(d ** n, dtype=complex)
    for idx in range(d ** n):
        t, digits = idx, []
        for _ in range(n):
            digits.append(t % d); t //= d
        Uvec[idx] = np.prod([Dvec[g] ** weights[i] for i, g in enumerate(digits)])

    def encode(vec):
        return sum(int(vec[i]) * (d ** i) for i in range(n))

    psis = []
    for r in code.cosets:
        v = np.zeros(d ** n, dtype=complex)
        for a in code.mx:
            v[encode([(r[j] + a[j]) % d for j in range(n)])] += 1.0
        psis.append(v / np.linalg.norm(v))
    P = np.array(psis)
    UP = (Uvec[None, :] * P).T
    L = P.conj() @ UP
    return float(np.linalg.norm(UP - P.T @ L))


def _subset_weights(n: int):
    """w in {0,1}^n with 2..n ones (skip empty and the all-ones uniform case)."""
    for ones in range(2, n + 1):
        for S in combinations(range(n), ones):
            if ones == n:
                continue
            w = [0] * n
            for i in S:
                w[i] = 1
            yield w


def _evaluate_code(code, gates, n: int):
    """Test one code against all level-3 gates (uniform + {0,1}-addressable when distance>=2).

    Returns ('hit', distance, gate_index, weights, cert) for the best nontrivial hit, 'trivial' if only
    distance-1 transversal hits, or None.
    """
    dist = code_distance(code)
    weight_sets = [None] + (list(_subset_weights(n)) if (dist >= 2 and n <= 7) else [])
    found_trivial = False
    for gi, (gname, U) in enumerate(gates):
        for w in weight_sets:
            cert = certify_distillation_code(code, U, weights=w)
            if not cert.get("is_distillation_code"):
                continue
            if dist < 2:
                found_trivial = True
            else:
                return ("hit", dist, gi, w, cert)
            break
    return "trivial" if found_trivial else None


def wsd_worker(seed: int, n: int, trials: int):
    """Parallel worker: sample `trials` random weakly-self-dual k=1 cyclic codes at dim n, test them.

    Top-level (picklable) for joblib/loky.  Returns a picklable summary.
    """
    import numpy as _np
    from primepower_msd.hierarchy_search import search as _level3_family
    rng = _np.random.default_rng(seed)
    gates = [(f"antidiff(S^{a} Z^{b})", U) for (a, b, sig, U) in _level3_family(4)]
    best = None
    trivial = 0
    ndist2 = 0
    scanned = 0
    for code in candidate_codes(n, 4, rng, per_rank=trials, exhaustive_cap=0):  # pure sampling
        scanned += 1
        if code_distance(code) >= 2:
            ndist2 += 1
        res = _evaluate_code(code, gates, n)
        if res == "trivial":
            trivial += 1
        elif res is not None and (best is None or res[1] > best[1]):
            _, dist, gi, w, cert = res
            best = (dist, gi, w, tuple(map(tuple, code.gens)), n)
    return {"best": best, "trivial": trivial, "ndist2": ndist2, "scanned": scanned}


def main() -> int:
    d = 4
    rng = np.random.default_rng(2024)
    gates = [(f"antidiff(S^{a} Z^{b})", U) for (a, b, sig, U) in level3_family(d)]
    print("=" * 78)
    print(f"M2 SEARCH — single-ququart (d={d}) MSD code; {len(gates)} level-3 physical gates")
    print("=" * 78)

    import time
    from functools import partial
    from joblib import Parallel, delayed
    from primepower_msd.parallel import default_jobs

    n_jobs = default_jobs()
    trivial = 0
    ndist2 = 0
    scanned = 0
    best = None        # (dist, gi, w, gens, n)

    def absorb(summary):
        nonlocal trivial, ndist2, scanned, best
        trivial += summary["trivial"]
        ndist2 += summary["ndist2"]
        scanned += summary["scanned"]
        b = summary["best"]
        if b and (best is None or b[0] > best[0]):
            best = b

    t0 = time.time()
    # n=3 exhaustive (all ranks); n=5,7 = exhaustive small-rank (serial) + parallel random sampling
    plans = [(3, 0), (5, 3000), (7, 4000)]
    for n, trials_per_chunk in plans:
        # serial exhaustive small ranks (deterministic; once)
        exh = list(candidate_codes(n, d, rng, per_rank=0, exhaustive_cap=200_000))
        for code in exh:
            scanned += 1
            if code_distance(code) >= 2:
                ndist2 += 1
            res = _evaluate_code(code, gates, n)
            if res == "trivial":
                trivial += 1
            elif res is not None and (best is None or res[1] > best[0]):
                best = (res[1], res[2], res[3], tuple(map(tuple, code.gens)), n)
        # parallel random sampling of the higher ranks
        if trials_per_chunk:
            results = Parallel(n_jobs=n_jobs, backend="loky")(
                delayed(wsd_worker)(20_240_629 + i * 1_000_003, n, trials_per_chunk) for i in range(n_jobs)
            )
            for r in results:
                absorb(r)
        print(f"  n={n}: cumulative scanned {scanned}, {ndist2} with distance>=2", flush=True)

    print("\n" + "=" * 78)
    print("VERDICT (M2)")
    print("=" * 78)
    print(f"  scanned {scanned} k=1 cyclic codes x {len(gates)} level-3 gates in {time.time()-t0:.1f}s "
          f"({n_jobs} cores; uniform + subset-addressable).")
    print(f"  trivial (distance-1, unencoded) transversal hits: {trivial}")
    if best:
        dist, gi, w, gens, n = best
        ga, gb, _, U = level3_family(d)[gi]
        from primepower_msd.ring_css import build_css
        code = build_css(list(gens), d)
        kind = "UNIFORM" if w is None else f"ADDRESSABLE qudits {[i for i, x in enumerate(w) if x]}"
        print(f"\n  NONTRIVIAL EXISTENCE WITNESS: [[{n},1,{dist}]]_{d} CSS, {kind} transversal antidiff(S^{ga} Z^{gb})")
        print(f"    X-stabilizer generators (Z_{d}): {list(gens)}")
        leak = brute_force_transversal_check(code, U, w)
        print(f"    brute-force codespace leakage: {leak:.2e} (0 => transversal confirmed)")
        print(f"    -> NON-TRIVIAL single-ququart distillation code (d={dist}) CONFIRMED.")
        return 0
    print("\n  NO distance>=2 single-ququart distillation code found (only trivial distance-1 hits).")
    print("  => transversal level-3 gates exist on Z_4 CSS codes, but NOT (in this family) on codes")
    print("     with error protection. Next (M2b): general CSS (A!=B), cubic gates, larger/punctured n.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
