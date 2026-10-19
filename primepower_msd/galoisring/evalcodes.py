"""Step 1 of the diverse-gate plan: FULL-ALPHABET ring codes from Galois-ring evaluation (Reed-Solomon /
Reed-Muller type) constructions over GR(4, m), as the ring analog of Krishna-Tillich's prime-field RS
triorthogonal codes ("does a larger alphabet help, as for primes?").

Construction.  R = GR(4, m), points = Teichmueller set T (2^m elements, incl. 0) or a subset; the
evaluation vector ev(f) = (f(t))_{t in points} of a polynomial f in R[x] is a vector of n ring elements,
each a tuple of m Z_4 components, i.e. a vector in Z_4^{m n} on N = m n ququarts.  CSS code:
    M_X  = span{ ev(x^j) * w : j in stabilizer band }   (component embedding; optional unit weights w)
    L    = ev(x^j) for j in the logical band            (0/1-free, full alphabet)
    M_Z  = (M_X + L)^perp over Z_4
plus an optional 2-torsion X-layer 2 * ev(x^j), j in a torsion band.  The transversal single-qudit layer
is q (the carrier) or a per-qudit choice among {q, q o neg} (coordinate sign flips are Clifford
relabelings), optimized by a randomized local search on the S1-S3 violation count.

Everything is evaluated with the generic machinery: rgt.hosts (exact for CSS), induced logical escape,
z4_distance MITM, gamma.  Results -> runs/diverse/catalog_greval.json.
"""

from __future__ import annotations

import itertools as it
import json
import math
import random
import time

import numpy as np

from .. import phase_tables as PT
from .. import rgt
from ..ringlinalg import right_kernel_howell
from ..symplectic import Z4Code
from ..z4_distance import code_distance
from .engine import GR, gr


# --------------------------------------------------------------------------- evaluation vectors
def eval_rows(R: GR, points, degrees, weights=None, torsion: bool = False) -> np.ndarray:
    """Component vectors (length m n) of ev(x^j) (times unit weights w_t), j in degrees."""
    rows = []
    for j in degrees:
        comps = []
        for idx, t in enumerate(points):
            v = R.pow(t, j) if j > 0 else R.one
            if weights is not None:
                v = R.mul(v, weights[idx])
            comps.append(R.xpart(v))
        row = np.concatenate(comps) % 4
        rows.append((2 * row) % 4 if torsion else row)
    return np.array(rows, dtype=np.int64).reshape(len(degrees), -1)


def build_eval_css(R: GR, points, stab_degs, log_degs, torsion_degs=(), weights=None, label: str = "") -> Z4Code | None:
    MX = eval_rows(R, points, stab_degs, weights) if len(stab_degs) else np.zeros((0, R.m * len(points)), dtype=np.int64)
    if len(torsion_degs):
        MX = np.vstack([MX, eval_rows(R, points, torsion_degs, weights, torsion=True)])
    L = eval_rows(R, points, log_degs, weights)
    N = R.m * len(points)
    span = np.vstack([MX, L]) if MX.size else L
    MZ = right_kernel_howell(span % 4, 4)
    MZ = np.array([z for z in MZ if np.any(np.array(z) % 4)], dtype=np.int64).reshape(-1, N)
    try:
        code = Z4Code.from_css(MX, MZ, N, 4, label or f"GR(4,{R.m}) eval n={len(points)} S{list(stab_degs)} L{list(log_degs)}")
    except AssertionError:
        return None
    code.eval_info = {"m": R.m, "n_points": len(points), "stab_degs": list(stab_degs), "log_degs": list(log_degs),
                      "torsion_degs": list(torsion_degs)}
    return code


# --------------------------------------------------------------------------- per-qudit sign choice
def sign_layer(n: int, q, signs) -> rgt.TransversalLayer:
    q = np.asarray(q, dtype=np.int64)
    qneg = q[(-np.arange(4)) % 4]
    return rgt.TransversalLayer.per_qudit([q if s == 0 else qneg for s in signs], 4, "signed carrier")


def violations(code: Z4Code, q, signs) -> int:
    """S1 + S2 + S3 violation count for the signed layer (x = 0 coset conditions)."""
    layer = sign_layer(code.n, q, signs)
    xg = [np.array(g) for g in code.x_gens]
    zero = np.zeros(code.n, dtype=np.int64)
    bad = 0
    for a in xg:
        bad += int(layer.phi(a) != layer.phi(zero))
    for i in range(len(xg)):
        for j in range(i, len(xg)):
            bad += int((layer.phi((xg[i] + xg[j]) % 4) - layer.phi(xg[i]) - layer.phi(xg[j]) + layer.phi(zero)) % 16 != 0)
    return bad


def optimize_signs(code: Z4Code, q, rng: random.Random, iters: int = 300) -> tuple[list[int], int]:
    n = code.n
    signs = [0] * n
    best = violations(code, q, signs)
    if best == 0:
        return signs, 0
    for _ in range(iters):
        i = rng.randrange(n)
        signs[i] ^= 1
        v = violations(code, q, signs)
        if v <= best:
            best = v
            if best == 0:
                break
        else:
            signs[i] ^= 1
    return signs, best


# --------------------------------------------------------------------------- evaluation
def evaluate_code(code: Z4Code, rng: random.Random, d_max: int = 4, tables=None) -> dict:
    t0 = time.time()
    rep = {"label": code.label, "n": code.n, "k": code.k, "k_genuine": code.k_genuine, "soft": code.soft,
           "free": code.free_quotient_ok, "info": getattr(code, "eval_info", {})}
    tables = tables if tables is not None else [np.array(PT.CARRIER_Q)]
    hits = []
    for q in tables:
        screens_uniform = rgt.uniform_screens(code, q)
        signs, viol = optimize_signs(code, q, rng) if not screens_uniform["all"] else ([0] * code.n, 0)
        if viol != 0:
            continue
        layer = sign_layer(code.n, q, signs)
        h = rgt.hosts(code, layer, "support")
        if not h.hosted:
            continue
        entry = {"q": [int(v) for v in q], "signs": int(sum(signs)), "hosted": True}
        if 1 <= code.k_genuine <= 6:
            cl = PT.classify(rgt.induced_logical(code, layer))
            entry["logical_escapes"] = cl["escapes_A2"]; entry["logical_level"] = cl["level"]
        hits.append(entry)
    rep["hosting_hits"] = hits
    rep["any_escaping_logical"] = any(e.get("logical_escapes") for e in hits)
    if hits:
        try:
            dmin, cert, det = code_distance(code, d_max=d_max)
        except MemoryError:
            dmin, cert, det = None, False, {}
        rep.update({"distance": dmin, "distance_certified": cert, "distance_detail": det})
        rep["gamma"] = (math.log(code.n / code.k) / math.log(dmin)) if (dmin and dmin > 1 and code.k > 0) else None
    rep["seconds"] = round(time.time() - t0, 1)
    return rep


def sweep(m_values=(2, 3), seed: int = 0, max_points: int = 16, d_max: int = 4, out: str | None = None,
          verbose: bool = True) -> list[dict]:
    rng = random.Random(seed)
    tables = [np.array(q) for q in PT.strict_level3_single(4) if __import__("primepower_msd.a2", fromlist=["escapes_A2"]).escapes_A2(np.array(q))]
    results = []
    for m in m_values:
        R = gr(m)
        T = R.teichmueller()
        point_sets = [("T", T), ("T*", T[1:])]
        for name, pts in point_sets:
            if len(pts) > max_points:
                continue
            n = len(pts)
            for r0 in range(1, min(n, 5)):
                for r1 in range(r0 + 1, min(n, r0 + 4)):
                    for tors in ([], list(range(r1, min(n, r1 + 2)))):
                        code = build_eval_css(R, pts, list(range(r0)), list(range(r0, r1)), tors,
                                              label=f"GR(4,{m}) {name} n={n} S<{r0} L[{r0},{r1}) T{tors}")
                        if code is None or code.k_genuine < 1:
                            continue
                        rep = evaluate_code(code, rng, d_max=d_max, tables=tables[:8])
                        results.append(rep)
                        if verbose:
                            print({k: v for k, v in rep.items() if k not in ("info",)}, flush=True)
    if out:
        with open(out, "w") as f:
            json.dump(results, f, indent=1, default=str)
    return results


if __name__ == "__main__":
    import sys
    ms = tuple(int(x) for x in sys.argv[1].split(",")) if len(sys.argv) > 1 else (2, 3)
    sweep(ms, out="runs/diverse/catalog_greval.json")
