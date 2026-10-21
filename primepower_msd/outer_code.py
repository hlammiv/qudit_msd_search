"""Outer codes for HHPW-style distillation over Z_d (sensitivity is alphabet-independent via unique neighbours).

s = 2 (target order 5): every input lies in two checks, so the outer code is a w-regular graph G on the checks
(inputs = edges).  For an error set V of edges, |Mv| >= #vertices of V-degree 1, and with adversarial Z_d values
non-unique checks can cancel.  The exact condition 2|Mv| + |v| >= 5 for |v| <= 4 holds for every d iff G is
simple with no triangles and no 4-cycles (girth >= 5):  double edge (|v|=2, U=0), triangle (|v|=3, U=0) and
4-cycle (|v|=4, U=0) violate it; every other configuration of <= 4 edges has 2U + |v| >= 5 (checked
exhaustively in tests).  The number of inputs is n_outer = m w / 2 for m checks; a w-regular girth-5 graph
needs m >= w^2 + 1 (Moore bound).
"""

from __future__ import annotations

import itertools as it

import numpy as np


def girth5_regular(m: int, w: int, seed: int = 0, tries: int = 20, max_stall: int = 200000):
    """Randomised greedy + edge switching for a w-regular graph on m vertices with girth >= 5.
    Returns a list of adjacency sets or None."""
    rng = np.random.default_rng(seed)
    for _ in range(tries):
        adj = [set() for _ in range(m)]

        def ok(u, v):
            if u == v or v in adj[u]:
                return False
            Nu, Nv = adj[u], adj[v]
            if Nu & Nv:                       # triangle / 4-cycle through a common neighbour
                return False
            for a in Nu:                      # 4-cycle u-a-b-v
                if adj[a] & Nv:
                    return False
            return True

        for _ in range(max_stall):
            open_ = [u for u in range(m) if len(adj[u]) < w]
            if not open_:
                return adj
            u = open_[rng.integers(len(open_))]
            cand = [v for v in open_ if ok(u, v)]
            if cand:
                v = cand[rng.integers(len(cand))]
                adj[u].add(v); adj[v].add(u)
                continue
            # stuck at u: remove a random edge at u's 2-neighbourhood blocker and retry
            a = list(adj[u])[rng.integers(len(adj[u]))] if adj[u] else int(rng.integers(m))
            if adj[a]:
                b = list(adj[a])[rng.integers(len(adj[a]))]
                adj[a].discard(b); adj[b].discard(a)
        # failed; retry
    return None


def check_girth5(adj) -> bool:
    for u, Nu in enumerate(adj):
        for v in Nu:
            if adj[v] & Nu:
                return False
        for a, b in it.combinations(Nu, 2):
            if len((adj[a] & adj[b]) - {u}) > 0:
                return False
    return True


def sensitivity_graph(adj, D: int = 5) -> bool:
    """Exhaustive check of 2 U + |V| >= D over all connected edge sets V with |V| <= D-1 (small graphs only)."""
    edges = sorted({(min(u, v), max(u, v)) for u in range(len(adj)) for v in adj[u]})
    eidx = {e: i for i, e in enumerate(edges)}
    inc = [[] for _ in adj]
    for i, (u, v) in enumerate(edges):
        inc[u].append(i); inc[v].append(i)

    def U(V):
        deg = {}
        for i in V:
            for x in edges[i]:
                deg[x] = deg.get(x, 0) + 1
        return sum(1 for c in deg.values() if c == 1)

    seen = set()

    def grow(V):
        key = frozenset(V)
        if key in seen:
            return True
        seen.add(key)
        if 2 * U(V) + len(V) < D:
            return False
        if len(V) == D - 1:
            return True
        verts = {x for i in V for x in edges[i]}
        for x in verts:
            for j in inc[x]:
                if j not in V and not grow(V | {j}):
                    return False
        return True

    return all(grow(frozenset([i])) for i in range(len(edges)))


def biaffine(q: int, xs=None, ms=None):
    """Biaffine-plane incidence graph over prime q: points (x,y) with x in xs, lines [m,b] with m in ms,
    (x,y) ~ [m,b] iff y = m x + b.  |xs| = |ms| = w gives a w-regular bipartite graph of girth >= 6 on 2wq vertices."""
    xs = list(range(q)) if xs is None else list(xs)
    ms = list(range(q)) if ms is None else list(ms)
    P = {(x, y): i for i, (x, y) in enumerate(it.product(xs, range(q)))}
    L = {(m, b): len(P) + i for i, (m, b) in enumerate(it.product(ms, range(q)))}
    adj = [set() for _ in range(len(P) + len(L))]
    for (x, y), i in P.items():
        for m in ms:
            j = L[(m, (y - m * x) % q)]
            adj[i].add(j); adj[j].add(i)
    return adj


def augment_girth5(adj, w: int, seed: int = 0, iters: int = 2_000_000):
    """Add edges to reach w-regular while keeping girth >= 5 (random greedy with local switches)."""
    rng = np.random.default_rng(seed)
    adj = [set(a) for a in adj]
    base = {(min(u, v), max(u, v)) for u in range(len(adj)) for v in adj[u]}
    m = len(adj)

    def ok(u, v):
        if u == v or v in adj[u] or adj[u] & adj[v]:
            return False
        for a in adj[u]:
            if adj[a] & adj[v]:
                return False
        return True

    for _ in range(iters):
        open_ = [u for u in range(m) if len(adj[u]) < w]
        if not open_:
            return adj
        u = open_[rng.integers(len(open_))]
        cand = [v for v in open_ if ok(u, v)]
        if cand:
            v = cand[rng.integers(len(cand))]
            adj[u].add(v); adj[v].add(u)
            continue
        added = [b for b in adj[u] if (min(u, b), max(u, b)) not in base]
        pool = [a for a in range(m) if any((min(a, b), max(a, b)) not in base for b in adj[a])]
        a = added[rng.integers(len(added))] if added and rng.random() < 0.5 else (pool[rng.integers(len(pool))] if pool else None)
        if a is None:
            return None
        extra = [b for b in adj[a] if (min(a, b), max(a, b)) not in base]
        if extra:
            b = extra[rng.integers(len(extra))]
            adj[a].discard(b); adj[b].discard(a)
    return None
