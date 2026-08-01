"""Structure-aware puncture-set sampling: cap sets (no 3 collinear points).

Reed-Muller minimum-weight codewords are supported on affine lines, so puncturing three
*collinear* points creates a short dual codeword that collapses the code distance. A **cap
set** -- a point set in F_p^m with no 3 points collinear -- structurally avoids this, and
the paper's own search codes are caps (the [[72,9,3]]_3 puncture points form a 9-cap in
F_3^4). Sampling biased toward cap sets reliably reaches the rare high-distance puncture
sets where uniform-random sampling stalls at d=2. See ../SAMPLING_INVESTIGATION.md.

Three points a, b, c are collinear in F_p^m iff (c-a) is parallel to (b-a); for p=3 this is
the familiar a+b+c == 0 (mod 3) line condition.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np

from .puncture import column_to_point, point_to_column


def collinear(a, b, c, p) -> bool:
    """True if distinct points a, b, c in F_p^m are collinear ((c-a) parallel to (b-a))."""
    m = len(a)
    u = [(b[i] - a[i]) % p for i in range(m)]
    v = [(c[i] - a[i]) % p for i in range(m)]
    s = None
    for i in range(m):
        if u[i] != 0:
            s = (v[i] * pow(u[i], p - 2, p)) % p  # scalar with v = s*u on this coord
            break
    if s is None:
        return True  # u == 0 means a == b (degenerate); treat as collinear
    return all((s * u[i] - v[i]) % p == 0 for i in range(m))


def cap_extends(chosen, x, p) -> bool:
    """Does adding point x to the cap ``chosen`` keep it a cap (no 3 collinear)?

    ``chosen`` must already be a cap; only new triples {a, b, x} can violate it.
    """
    n = len(chosen)
    for i in range(n):
        a = chosen[i]
        for j in range(i + 1, n):
            if collinear(a, chosen[j], x, p):
                return False
    return True


def is_cap(points, p) -> bool:
    """True if no 3 of ``points`` are collinear in F_p^m (full O(k^3) check)."""
    n = len(points)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if collinear(points[i], points[j], points[k], p):
                    return False
    return True


def all_points(m, p) -> list:
    """All p^m points of F_p^m as int tuples, ordered so entry c-1 <-> 1-indexed column c."""
    return [tuple(int(v) for v in column_to_point(c, m, p)) for c in range(1, p ** m + 1)]


def random_cap(m, p, k, rng, allpts=None):
    """A greedy random cap of ``k`` points (list of point tuples), or None if it stalls.

    Shuffles the points and greedily admits each one that keeps the cap property, until k
    are chosen. Returns None if a single greedy pass cannot reach size k (the caller retries
    with a fresh shuffle). Note caps have a bounded maximum size, so k must be within it.
    """
    if allpts is None:
        allpts = all_points(m, p)
    order = list(allpts)
    rng.shuffle(order)
    chosen: list = []
    for x in order:
        if cap_extends(chosen, x, p):
            chosen.append(x)
            if len(chosen) == k:
                return chosen
    return None


def points_to_columns(points, p) -> tuple:
    """Sorted 1-indexed columns for a list of F_p^m points."""
    return tuple(sorted(point_to_column(x, p) for x in points))


# --- plane-spread sampling: caps that additionally keep <=3 points per 2-flat -------------
# Empirical basis: a cap (no 3 collinear) kills the LINE-supported low-weight codewords, but at
# higher m the binding low-weight dual codewords sit on 2-flats. The paper's [[230,13,6]]_3 (d=6)
# puncture set has max 3 points on any 2-flat, whereas every random d=5 cap has 4. Enforcing
# "no 4 coplanar" reproduces the d=6 code deterministically (see tests).
def _fp_rank(rows, p) -> int:
    """Rank over F_p of a small integer matrix given as a list of rows."""
    M = np.array(rows, dtype=np.int64) % p
    nr, nc = M.shape
    rr = 0
    for c in range(nc):
        piv = next((i for i in range(rr, nr) if M[i, c] % p), None)
        if piv is None:
            continue
        M[[rr, piv]] = M[[piv, rr]]
        M[rr] = (M[rr] * pow(int(M[rr, c]), p - 2, p)) % p
        for i in range(nr):
            if i != rr and M[i, c] % p:
                M[i] = (M[i] - M[i, c] * M[rr]) % p
        rr += 1
        if rr == nr:
            break
    return rr


def four_coplanar(a, b, c, x, p) -> bool:
    """True if a, b, c, x lie on one affine 2-flat in F_p^m (difference vectors have rank <= 2)."""
    m = len(a)
    dirs = [[(b[i] - a[i]) % p for i in range(m)],
            [(c[i] - a[i]) % p for i in range(m)],
            [(x[i] - a[i]) % p for i in range(m)]]
    return _fp_rank(dirs, p) <= 2


def plane_spread_extends(chosen, x, p) -> bool:
    """Does adding x keep ``chosen`` a cap (no 3 collinear) AND no 4 points coplanar (<=3/2-flat)?"""
    n = len(chosen)
    for i in range(n):                       # cap: no 3 collinear
        for j in range(i + 1, n):
            if collinear(chosen[i], chosen[j], x, p):
                return False
    for i in range(n):                       # no 4 coplanar
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if four_coplanar(chosen[i], chosen[j], chosen[k], x, p):
                    return False
    return True


def new_collinear_count(chosen, x, p) -> int:
    """Number of collinear triples {a, b, x} that adding x to a cap-ish ``chosen`` would create."""
    n = len(chosen)
    cnt = 0
    for i in range(n):
        a = chosen[i]
        for j in range(i + 1, n):
            if collinear(a, chosen[j], x, p):
                cnt += 1
    return cnt


def random_near_cap(m, p, k, max_triples, rng, allpts=None):
    """A greedy random NEAR-cap of ``k`` points with at most ``max_triples`` collinear triples
    total, or None if it stalls. Near-caps build at k near the max-cap size where strict caps
    stall (e.g. k=43 in AG(5,3), whose max cap is 45) and reach codes like [[200,43,3]]_3."""
    if allpts is None:
        allpts = all_points(m, p)
    order = list(allpts)
    rng.shuffle(order)
    chosen: list = []
    used = 0
    for x in order:
        nc = new_collinear_count(chosen, x, p)
        if used + nc <= max_triples:
            chosen.append(x)
            used += nc
            if len(chosen) == k:
                return chosen
    return None


def random_plane_spread(m, p, k, rng, allpts=None):
    """A greedy random plane-spread cap of ``k`` points (no 3 collinear, no 4 coplanar), or None if
    it stalls. These sets carry the higher distance the plain cap misses (e.g. d=6 vs d=5 at m=5)."""
    if allpts is None:
        allpts = all_points(m, p)
    order = list(allpts)
    rng.shuffle(order)
    chosen: list = []
    for x in order:
        if plane_spread_extends(chosen, x, p):
            chosen.append(x)
            if len(chosen) == k:
                return chosen
    return None


# --- unified flat-spread / arc-order sampling ---------------------------------------------
# The cap / plane_spread / near_cap samplers are one idea at different strictness: bound the
# occupancy of every j-flat. An "order-t arc" keeps <=(j+1) points on every j-flat for
# j=1..t: order 1 == cap (no 3 collinear); order 2 == plane_spread (no 4 coplanar);
# order 3 == no 5 on a 3-flat (-> d=7). Higher order = stronger = higher distance but only
# builds at smaller k. random_flat_spread auto-picks the max feasible order for the given k,
# falling back to a near-cap when even a cap stalls (k near the max-cap size).
def _flat_occupied(chosen, x, p, j) -> bool:
    """True if adding x would put a (j+2)-th point on some j-flat that already holds j+1 chosen
    points (an affinely-independent (j+1)-subset spanning the flat, with x on it)."""
    m = len(x)
    for combo in combinations(chosen, j + 1):
        base = combo[0]
        dirs = [[(c[t] - base[t]) % p for t in range(m)] for c in combo[1:]]
        if _fp_rank(dirs, p) < j:                      # combo degenerate: not a spanning j-flat
            continue
        xd = [(x[t] - base[t]) % p for t in range(m)]
        if _fp_rank(dirs + [xd], p) <= j:              # x lies on the j-flat the combo spans
            return True
    return False


def arc_extends(chosen, x, p, order) -> bool:
    """Does adding x keep <=(j+1) points on every j-flat for j=1..order? (an order-`order` arc)."""
    n = len(chosen)
    for a in range(n):                                 # j=1 (lines): no 3 collinear
        for b in range(a + 1, n):
            if collinear(chosen[a], chosen[b], x, p):
                return False
    for j in range(2, order + 1):                      # j>=2: no (j+2) on a j-flat
        if _flat_occupied(chosen, x, p, j):
            return False
    return True


def random_arc(m, p, k, order, rng, allpts=None):
    """A greedy random order-`order` arc of ``k`` points, or None if it stalls.
    order 1 == cap, 2 == plane_spread, 3 == space-spread (d=7). Higher order builds only at
    smaller k but reaches higher distance."""
    if allpts is None:
        allpts = all_points(m, p)
    order_pts = list(allpts)
    rng.shuffle(order_pts)
    chosen: list = []
    for x in order_pts:
        if arc_extends(chosen, x, p, order):
            chosen.append(x)
            if len(chosen) == k:
                return chosen
    return None


_ORDER_CACHE: dict = {}


def max_feasible_order(m, p, k, allpts=None, max_order=3, attempts=4) -> int:
    """Highest arc order (<=max_order) at which a size-k arc builds; 0 if even a cap (order 1)
    stalls (k near/above the max-cap size). Cached per (m,p,k,max_order); deterministic."""
    key = (m, p, k, max_order)
    if key not in _ORDER_CACHE:
        if allpts is None:
            allpts = all_points(m, p)
        import random as _random
        rng = _random.Random(0x5EED ^ (k * 131 + p))
        found = 0
        for order in range(max_order, 0, -1):
            if any(random_arc(m, p, k, order, rng, allpts) is not None for _ in range(attempts)):
                found = order
                break
        _ORDER_CACHE[key] = found
    return _ORDER_CACHE[key]


def random_flat_spread(m, p, k, rng, allpts=None, max_order=3):
    """UNIFIED adaptive sampler: build the highest-order arc that fits size ``k`` (auto-selected,
    cached), or fall back to a near-cap with an auto-grown triple budget when even a cap stalls.
    Subsumes capset (order 1), plane_spread (order 2), and near_cap (fallback) as operating
    points; adds order 3 (d=7) at small k."""
    if allpts is None:
        allpts = all_points(m, p)
    order = max_feasible_order(m, p, k, allpts, max_order)
    if order >= 1:
        return random_arc(m, p, k, order, rng, allpts)
    budget = k                                          # cap stalls -> near-cap, grow budget
    for _ in range(6):
        pts = random_near_cap(m, p, k, budget, rng, allpts)
        if pts is not None:
            return pts
        budget *= 2
    return None
