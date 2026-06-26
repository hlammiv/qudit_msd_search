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
