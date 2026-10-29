"""Signed Z_d rotated planar surface code: memory experiment, detector-error model and a 2-adic matching decoder.

Code (d = 2^k ring qudits).  Data qudits on a D x D grid; checks on the faces (checkerboard: X-type when
row+col of the NW corner is even), weight-2 X checks on the left/right boundary, weight-2 Z checks on the
top/bottom boundary.  Signs: an X check acts as X^{(-1)^row}, a Z check as Z^{(-1)^col}; two overlapping checks
share one row-pair or one column-pair, so the signs cancel and the group is abelian over any Z_d.  Logical Z =
Z on one column (constant sign), logical X = X on one row; they share one qudit, so the logical qudit is a genuine
Z_d qudit.

Syndrome extraction: one ancilla per check; Z checks by CSUM(data -> ancilla, w); X checks by F, CSUM(ancilla ->
data, w), F^dag.  Schedule X: NW NE SW SE, Z: NW SW NE SE (the usual N / Z patterns): overlapping checks touch the
shared pair in the same order, and hook errors run perpendicular to the logical they could flip.

Noise (every location is one "p"): DEP1 on data while the ancillas reset, XERR after reset, DEP1 after F / F^dag,
DEP2 after every CSUM, XERR before every measurement.  The same convention (uniform non-identity Pauli with total
probability p per location) is used by the distillation counts, so the thresholds here combine with them directly.

Decoder: the detector-error model is built from the circuit itself by propagating one basis fault (X^1 or Z^1 on
one qudit) per noise location through the frame simulator; every non-identity Pauli at the location is then a Z_d
combination of those.  Faults with <= 2 detectors in a sector become edges; the Z_d coefficient pattern along the
edge (c_a, c_b, o) says how a unit flow moves charge and shifts the observable.  Decoding peels 2-adically: level j
(j = 0..k-1) matches the detectors whose charge has 2-adic valuation j with binary MWPM (PyMatching; edge weights
from the probability of a fault of that valuation on the edge), walks each matched pair along a shortest path
carrying the Z_d charge, and leaves the residue (now divisible by 2^{j+1}) to the next level.  Hyperedges are left
out of the graph (their probability share is reported).
"""
from __future__ import annotations

import itertools as it
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

from .zdsim import Circuit, FrameSim, check_determinism

STEPS_X = {"NW": 0, "SW": 1, "NE": 2, "SE": 3}
STEPS_Z = {"NW": 0, "NE": 1, "SW": 2, "SE": 3}


@dataclass
class Check:
    typ: str                      # 'X' or 'Z'
    legs: list                    # [(qudit, weight, step)]


def rotated_code(D: int, d: int):
    """Checks of the signed Z_d rotated surface code with distance D."""
    checks = []
    for i in range(-1, D):
        for j in range(-1, D):
            corners = {"NW": (i, j), "NE": (i, j + 1), "SW": (i + 1, j), "SE": (i + 1, j + 1)}
            inside = {k: (0 <= r < D and 0 <= c < D) for k, (r, c) in corners.items()}
            n_in = sum(inside.values())
            typ = "X" if (i + j) % 2 == 0 else "Z"
            if n_in == 2:
                horizontal = i in (-1, D - 1)          # top/bottom boundary -> pair in one row
                if (horizontal and typ == "X") or (not horizontal and typ == "Z"):
                    continue
            elif n_in != 4:
                continue
            steps = STEPS_X if typ == "X" else STEPS_Z
            legs = []
            for k, (r, c) in corners.items():
                if not inside[k]:
                    continue
                w = (-1) ** (r if typ == "X" else c)
                legs.append((r * D + c, w % d, steps[k]))
            checks.append(Check(typ, legs))
    return checks


def verify_code(D: int, d: int) -> dict:
    """Commutation of all checks and of the logicals with the checks (mod d)."""
    checks = rotated_code(D, d)
    n = D * D
    X = [c for c in checks if c.typ == "X"]
    Z = [c for c in checks if c.typ == "Z"]
    bad = 0
    for cx in X:
        wx = dict((q, w) for q, w, _ in cx.legs)
        for cz in Z:
            s = sum(wx.get(q, 0) * w for q, w, _ in cz.legs) % d
            bad += s != 0
    logZ = {i * D for i in range(D)}                   # column 0
    logX = {j for j in range(D)}                       # row 0
    badL = sum(sum(w for q, w, _ in cx.legs if q in logZ) % d != 0 for cx in X)
    badL += sum(sum(w for q, w, _ in cz.legs if q in logX) % d != 0 for cz in Z)
    return dict(n=n, nX=len(X), nZ=len(Z), k=n - len(checks), noncommuting=bad, logical_violations=badL)


@dataclass
class Memory:
    circ: Circuit
    d: int
    D: int
    rounds: int
    basis: str
    det_sector: np.ndarray          # 'X' (detects X errors, i.e. Z checks) or 'Z' per detector
    obs: list                       # [(meas_index, coeff)] -- the logical measured at the end
    noise_locs: list = field(default_factory=list)  # (op_index, kind, qudits)


def memory_circuit(D: int, d: int, rounds: int, p: float, basis: str = "Z", idle: bool = True,
                   code_capacity: bool = False) -> Memory:
    """code_capacity: one DEP1(p) layer on the data before a single noiseless syndrome round (perfect measurement)."""
    checks = rotated_code(D, d)
    c = Circuit(d)
    data = c.q(D * D)
    anc = c.q(len(checks))
    xanc = [anc[i] for i, ch in enumerate(checks) if ch.typ == "X"]
    sector = {"X": "Z", "Z": "X"}       # a Z-type check detects X errors (sector 'X')
    det_sector = []
    locs = []

    def noise(kind, qs, pr, always=False):
        if pr > 0 and (always or not code_capacity):
            getattr(c, kind)(qs, pr)
            locs.append((len(c.ops) - 1, kind, list(qs)))

    # data preparation
    c.R(data); noise("XERR", data, p)
    if basis == "X":
        c.F(data); noise("DEP1", data, p)
    if code_capacity:
        noise("DEP1", data, p, always=True)
    prev = None
    for r in range(rounds):
        if idle and r > 0:
            noise("DEP1", data, p)
        c.R(anc); noise("XERR", anc, p)
        c.F(xanc); noise("DEP1", xanc, p)
        for step in range(4):
            pairs, ws = [], []
            for i, ch in enumerate(checks):
                for q, w, s in ch.legs:
                    if s != step:
                        continue
                    if ch.typ == "Z":
                        pairs.append((data[q], anc[i])); ws.append(w)
                    else:
                        pairs.append((anc[i], data[q])); ws.append(w)
            c.CSUM(pairs, ws); noise("DEP2", pairs, p)
        c.FD(xanc); noise("DEP1", xanc, p)
        noise("XERR", anc, p)
        m = c.M(anc)
        for i, ch in enumerate(checks):
            if prev is None:
                if ch.typ == basis:          # deterministic in the first round
                    c.detector([(m[i], 1)], f"{ch.typ}{i} r0"); det_sector.append(sector[ch.typ])
            else:
                c.detector([(m[i], 1), (prev[i], -1)], f"{ch.typ}{i} r{r}"); det_sector.append(sector[ch.typ])
        prev = m
    # final data measurement in the memory basis
    if basis == "X":
        c.FD(data); noise("DEP1", data, p)
    noise("XERR", data, p)
    md = c.M(data)
    for i, ch in enumerate(checks):
        if ch.typ != basis:
            continue
        sign = 1 if basis == "Z" else -1         # Z ancilla records +sum w x; X ancilla records -sum w z
        c.detector([(md[q], sign * w) for q, w, _ in ch.legs] + [(prev[i], -1)], f"{ch.typ}{i} final")
        det_sector.append(sector[ch.typ])
    if basis == "Z":
        obs = [(md[i * D], 1) for i in range(D)]          # logical Z = column 0
    else:
        obs = [(md[j], 1) for j in range(D)]              # logical X = row 0 (after F^dag, x records z)
    return Memory(c, d, D, rounds, basis, np.array(det_sector), obs, locs)


# ----------------------------------------------------------------------------------------------------------------
# detector-error model from basis faults


def _basis_faults(mem: Memory):
    """One (op_index, qudit, axis) basis fault per Pauli generator at each noise location; returns also, for each
    location, the list of its basis-fault indices in the order (x_q1, z_q1, x_q2, z_q2) / (x_q) for XERR."""
    faults, per_loc = [], []
    for k, kind, qs in mem.noise_locs:
        if kind == "DEP2":
            for a, b in qs:
                ids = []
                for q in (a, b):
                    for ax in ("x", "z"):
                        ids.append(len(faults)); faults.append((k, q, ax))
                per_loc.append((kind, ids))
        else:
            for q in qs:
                ids = []
                axes = ("x", "z") if kind == "DEP1" else (("x",) if kind == "XERR" else ("z",))
                for ax in axes:
                    ids.append(len(faults)); faults.append((k, q, ax))
                per_loc.append((kind, ids))
    return faults, per_loc


def propagate_basis(mem: Memory, chunk: int = 2000):
    """Detector vector (mod d) and observable coefficient of every basis fault, by frame-sim runs with one fault per
    shot (noise off, gauge off), in chunks of `chunk` faults to bound the dense forced-fault arrays."""
    faults, per_loc = _basis_faults(mem)
    c, d = mem.circ, mem.d
    V = np.zeros((len(faults), len(c.detectors)), dtype=np.int16)      # detector values mod d; int16 keeps D = 9 two-patch models in a few hundred MB
    O = np.zeros(len(faults), dtype=np.int64)
    for lo in range(0, len(faults), chunk):
        part = faults[lo:lo + chunk]
        forced = {}
        for s, (k, q, ax) in enumerate(part):
            forced.setdefault(k, []).append((s, q, ax))
        forced_arr = {}
        for k, lst in forced.items():
            qs = sorted(set(q for _, q, _ in lst))
            col = {q: i for i, q in enumerate(qs)}
            P = np.zeros((len(part), len(qs), 2), dtype=np.int64)
            mask = np.zeros(len(part), dtype=np.int64)
            for s, q, ax in lst:
                P[s, col[q], 0 if ax == "x" else 1] = 1; mask[s] = 1
            forced_arr[k] = [(mask, qs, P)]
        sim = FrameSim(c, len(part), np.random.default_rng(0), gauge=False).run(scale=0.0, forced=forced_arr)
        V[lo:lo + len(part)] = sim.detector_values()
        for i, cf in mem.obs:
            O[lo:lo + len(part)] += cf * sim.meas[:, i]
    return faults, per_loc, V, O % d


def _val2(x: int, d: int) -> int:
    """2-adic valuation of x mod d (= log2 d for x == 0)."""
    x %= d
    if x == 0:
        return int(round(math.log2(d)))
    v = 0
    while x % 2 == 0:
        x //= 2; v += 1
    return v


@dataclass
class DEM:
    d: int
    n_det: int
    edges: dict                       # (a, b) with a < b or b == -1 -> {'pat': (c_a, c_b, o), 'P': [per level]}
    hyper_prob: float
    total_prob: float
    undetectable_prob: float
    conflicts: int


def _same_class(a, b, d):
    return any(all((lam * x - y) % d == 0 for x, y in zip(a, b)) for lam in range(1, d))


def build_dem(mem: Memory, sector: str, exclude=None) -> DEM:
    """Project every fault onto the detectors of one sector; group graphlike faults into edges with per-level
    probabilities at unit location probability (scale by p when decoding; uniform over the non-identity Paulis of
    the location).  A fault with > 2 detectors is decomposed into its single-generator components (X^a or Z^b on one qudit); each
    graphlike component's edge receives the fault's probability (stim's decompose_errors).  When several fault
    classes share an edge with patterns that are not scalar multiples (a boundary degeneracy), the most probable
    class fixes the pattern.  exclude: detector indices that are postselected; every fault firing one of them is left out
    of the model (those shots are discarded), and its probability is reported as DEM.discarded_prob."""
    faults, per_loc, V, O = propagate_basis(mem)
    Vx = V[:, exclude] if exclude is not None and len(exclude) else None
    discarded = 0.0
    d = mem.d
    k = int(round(math.log2(d)))
    sel = np.nonzero(mem.det_sector == sector)[0]
    Vs = V[:, sel]
    n_det = len(sel)
    edges, hyper, total, undet, conflicts = {}, 0.0, 0.0, 0.0, 0

    def add(v, o, pr):
        """v: detector vector (sector), o: observable coefficient.  Returns False if not graphlike."""
        nonlocal undet
        dets = np.nonzero(v)[0]
        if len(dets) == 0:
            if o % d:
                undet += pr
            return True
        if len(dets) > 2:
            return False
        if len(dets) == 1:
            a, b = int(dets[0]), -1
            ca, cb = int(v[a]), 0
        else:
            a, b = int(dets[0]), int(dets[1])
            ca, cb = int(v[a]), int(v[b])
        lev = _val2(ca, d)
        if b != -1 and _val2(cb, d) != lev:
            return False
        e = edges.setdefault((a, b), {"P": [0.0] * k, "classes": []})
        e["P"][lev] += pr
        pat = (ca, cb, int(o) % d)
        for cl in e["classes"]:
            if _same_class(cl[0], pat, d):
                cl[1] += pr
                if _val2(pat[0], d) < _val2(cl[0][0], d):
                    cl[0] = pat
                break
        else:
            e["classes"].append([pat, pr])
        return True

    leftover = []                      # hyperedge faults (v, o, pr) for the second pass
    singles = {}                       # every kept single-fault pattern ((det, value), ...) -> {observable: probability}
    for kind, ids in per_loc:
        m = len(ids)
        combos = np.array(list(it.product(range(d), repeat=m)), dtype=np.int64)[1:]     # non-identity
        pr = 1.0 / len(combos)
        total += 1.0
        Vl = combos @ Vs[ids] % d
        Ol = combos @ O[ids] % d
        cnt = (Vl != 0).sum(axis=1)
        drop = np.any(combos @ Vx[ids] % d != 0, axis=1) if Vx is not None else np.zeros(len(combos), dtype=bool)
        for row in range(len(combos)):
            if drop[row]:
                discarded += pr
                continue
            if cnt[row]:
                key = tuple((int(j), int(Vl[row, j])) for j in np.nonzero(Vl[row])[0])
                sd = singles.setdefault(key, {}); sd[int(Ol[row])] = sd.get(int(Ol[row]), 0.0) + pr
            if cnt[row] <= 2 and add(Vl[row], Ol[row], pr):
                continue
            # decompose into single-generator components
            comps = [(combos[row, i] * Vs[ids[i]] % d, combos[row, i] * O[ids[i]] % d)
                     for i in range(m) if combos[row, i] % d]
            if all(len(np.nonzero(v)[0]) <= 2 for v, _ in comps):
                for v, o in comps:
                    add(v, o, pr)
            else:
                for v, o in comps:
                    if len(np.nonzero(v)[0]) <= 2:
                        add(v, o, pr)
                    else:
                        leftover.append((v, int(o), pr))
    for e in edges.values():
        e["classes"].sort(key=lambda cl: -cl[1])
        e["pat"] = e["classes"][0][0]
    # second pass (stim's decompose_errors): write each hyperedge as a sum of existing graphlike edges whose supports
    # partition its detector set, with consistent Z_d coefficients and observable; take the most probable decomposition
    def partitions(items):
        if not items:
            yield []
            return
        a, rest = items[0], items[1:]
        for p_ in partitions(rest):
            yield [(a, -1)] + p_
        for k, b in enumerate(rest):
            for p_ in partitions(rest[:k] + rest[k + 1:]):
                yield [(a, b)] + p_
    hyperedges = {}                    # undecomposed: ((det, value), ...) -> [observable, probability]
    def keep_hyper(v, o, pr, dets):
        key = tuple((j, int(v[j])) for j in dets)
        h = hyperedges.setdefault(key, [int(o), 0.0]); h[1] += pr
    for v, o, pr in leftover:
        dets = [int(x) for x in np.nonzero(v)[0]]
        if len(dets) > 4:
            hyper += pr; keep_hyper(v, o, pr, dets)
            continue
        best, best_lp = None, -math.inf
        for part in partitions(dets):
            if not all(key in edges for key in part):
                continue
            lams, lp, osum, ok = [], 0.0, 0, True
            for key in part:
                e = edges[key]; ca, cb, oe = e["pat"]
                lam = next((l for l in range(1, d) if (l * ca - v[key[0]]) % d == 0 and (key[1] == -1 or (l * cb - v[key[1]]) % d == 0)), None)
                if lam is None:
                    ok = False
                    break
                lams.append(lam); osum += lam * oe
                P = sum(e["P"])
                lp += math.log(P) if P > 0 else -50.0
            if ok and (osum - o) % d == 0 and lp > best_lp:
                best, best_lp = part, lp
        if best is None:
            hyper += pr; keep_hyper(v, o, pr, dets)
            continue
        lev = _val2(int(v[dets[0]]), d)
        for key in best:
            edges[key]["P"][lev] += pr
    for e in edges.values():
        e["classes"].sort(key=lambda cl: -cl[1])
        e["pat"] = e["classes"][0][0]
        conflicts += len(e["classes"]) - 1
    dem = DEM(d, n_det, edges, hyper, total, undet, conflicts)
    dem.discarded_prob = discarded
    dem.hyperedges = [(k, o, pr) for k, (o, pr) in hyperedges.items()]
    dem.singles = singles
    return dem


def misdecoded_patterns(dem: DEM, dec) -> dict:
    """First-order check of a decoder on `dem`: every kept single-fault pattern whose most probable observable the
    decoder does not reproduce (pattern -> probability per unit p).  Discarding shots that contain these patterns
    (`pattern_mask`) removes every first-order failure."""
    d, n = dem.d, dem.n_det
    bad = {}
    for key, od in getattr(dem, "singles", {}).items():
        o = max(od, key=od.get)
        sigma = np.zeros(n, dtype=np.int64)
        for j, val in key:
            sigma[j] = val
        corr, _ = dec.decode(sigma)
        if (corr + o) % d != 0:
            bad[key] = sum(od.values())
    return bad


def _pair_worker(dem, p, keys, singles, i_lo, i_hi, keys_b=None):
    dec = PeelDecoder(dem, p)
    d, n = dem.d, dem.n_det
    bad, undet, checked = {}, 0.0, 0
    obs = {k: max(od, key=od.get) for k, od in singles.items()}
    mass = {k: sum(od.values()) for k, od in singles.items()}
    a_index = {k: i for i, k in enumerate(keys)}
    for i in range(i_lo, i_hi):
        k1 = keys[i]
        # partner set: the rest of `keys` after i (triangle), plus every key of keys_b outside `keys`
        partners = keys[i + 1:] if keys_b is None else (keys[i + 1:] + [k for k in keys_b if k not in a_index])
        for k2 in partners:
            sigma = np.zeros(n, dtype=np.int64)
            for jj, val in k1:
                sigma[jj] += val
            for jj, val in k2:
                sigma[jj] += val
            sigma %= d
            o = (obs[k1] + obs[k2]) % d
            m = mass[k1] * mass[k2]
            checked += 1
            nz = np.nonzero(sigma)[0]
            if len(nz) == 0:
                if o:
                    undet += m
                continue
            corr, _ = dec.decode(sigma)
            if (corr + o) % d != 0:
                key = tuple((int(jj), int(sigma[jj])) for jj in nz)
                bad[key] = bad.get(key, 0.0) + m
    return bad, undet, checked


def misdecoded_pairs(dem: DEM, p: float, keys, n_jobs: int = 6, keys_b=None, exclude=None):
    """Second-order check of the peeling decoder: every unordered pair of the single-fault patterns `keys` (subset of
    dem.singles) — and, if keys_b is given, every pair (a in keys, b in keys_b outside keys) — is summed and decoded;
    returns (combined pattern -> probability mass per unit p^2 of the pairs that are mis-decoded, undetectable-pair mass,
    pairs checked).  Discarding shots that contain a bad combined pattern (`pattern_mask`) removes those second-order
    failures at a p^2 discard cost (plus first-order cost when the combined pattern is itself a single-fault pattern)."""
    from joblib import Parallel, delayed
    exclude = set(exclude or ())                  # e.g. the first-order mis-decoded singles (already postselected)
    keys = [k for k in keys if k not in exclude]
    keys_b = None if keys_b is None else [k for k in keys_b if k not in set(keys) and k not in exclude]
    singles = {k: dem.singles[k] for k in keys + (keys_b or [])}
    N = len(keys)
    extra = len(keys_b or [])
    # balance the workload: cut points with equal numbers of pairs
    cuts = [0]
    total = N * (N - 1) // 2 + N * extra
    nchunks = max(1, 4 * n_jobs)
    acc = 0
    for i in range(N):
        acc += N - 1 - i + extra
        if acc >= total * len(cuts) / nchunks and cuts[-1] < i + 1:
            cuts.append(i + 1)
    if cuts[-1] != N:
        cuts.append(N)
    res = Parallel(n_jobs=n_jobs)(delayed(_pair_worker)(dem, p, keys, singles, a, b, keys_b) for a, b in zip(cuts[:-1], cuts[1:]))
    bad, undet, checked = {}, 0.0, 0
    for b_, u_, c_ in res:
        for k, m in b_.items():
            bad[k] = bad.get(k, 0.0) + m
        undet += u_; checked += c_
    return bad, undet, checked


def pattern_mask(patterns, Dv: np.ndarray) -> np.ndarray:
    """Shots whose syndrome contains (exactly, on its own detectors) one of `patterns` (((det, value), ...))."""
    hit = np.zeros(len(Dv), dtype=bool)
    for key in patterns:
        idx = [j for j, _ in key]; vals = np.array([v for _, v in key])
        hit |= np.all(Dv[:, idx] == vals, axis=1)
    return hit


def hyperedge_mask(dem: DEM, Dv: np.ndarray) -> np.ndarray:
    """Shots containing one of the undecomposed hyperedge patterns of `dem`."""
    return pattern_mask([k for k, _, _ in getattr(dem, "hyperedges", [])], Dv)


# ----------------------------------------------------------------------------------------------------------------
# 2-adic peeling decoder


def boundary_potential(dem: DEM, lev: int):
    """Split the boundary for the complementary-gap decode at 2-adic level `lev`.  A potential x in {0,1}^n on the
    detectors with o(a,b) + x_a + x_b = 0 (mod 2) on every interior edge (o = the edge pattern's observable, odd classes
    taken mod 2) absorbs the interior observable flips; a boundary edge (a, B) then goes to boundary node B_label with
    label = x_a + o (mod 2), and every correction's logical parity equals the parity of its matchings to B_1 (plus a
    syndrome-only term common to all corrections).  Returns (x, labels {edge key: 0/1}, conflicts) where conflicts counts
    edges present at this level (P[lev] > 0) whose constraint cannot be met (undetectable cycles)."""
    n = dem.n_det
    adj = [[] for _ in range(n)]
    for (a, b), e in dem.edges.items():
        if b != -1 and e["P"][lev] > 0:
            o = e["pat"][2] & 1
            adj[a].append((b, o)); adj[b].append((a, o))
    x = -np.ones(n, dtype=np.int64)
    conflicts = 0
    for s0 in range(n):
        if x[s0] >= 0:
            continue
        x[s0] = 0
        stack = [s0]
        while stack:
            a = stack.pop()
            for b, o in adj[a]:
                want = (x[a] + o) % 2
                if x[b] < 0:
                    x[b] = want; stack.append(b)
                elif x[b] != want:
                    conflicts += 1
    labels = {}
    for (a, b), e in dem.edges.items():
        if b == -1:
            labels[(a, b)] = int((x[a] + (e["pat"][2] & 1)) % 2)
    return x, labels, conflicts // 2


class PeelDecoder:
    def __init__(self, dem: DEM, p: float, gap: bool = False):
        import pymatching
        self.d, self.n = dem.d, dem.n_det
        self.k = int(round(math.log2(self.d)))
        self.edges = dem.edges
        self.gap = gap
        self.trace = None                             # when a list: (edge key, multiplier) of every traversed edge
        B = self.n                                    # boundary node for the shortest-path graph
        self.matchers, self.pred, self.gap_matchers, self.potential, self.conflicts = [], [], [], [], []
        wmax = 0.0
        for lev in range(self.k):
            Ps = [p * e["P"][lev] for e in dem.edges.values() if e["P"][lev] > 0]
            wmax = max([math.log((1 - P) / P) for P in Ps] + [1.0])
            M = pymatching.Matching()
            MG = pymatching.Matching() if gap else None
            if gap:
                x, labels, conf = boundary_potential(dem, lev)
                self.potential.append(x); self.conflicts.append(conf)
            rows, cols, ws = [], [], []
            for (a, b), e in dem.edges.items():
                P = p * e["P"][lev]
                w = math.log((1 - P) / P) if P > 0 else 2 * wmax + 1.0
                if b == -1:
                    M.add_boundary_edge(a, weight=w, merge_strategy="smallest-weight")
                    rows.append(a); cols.append(B); ws.append(w)
                    if gap:
                        MG.add_edge(a, self.n + labels[(a, b)], weight=w, merge_strategy="smallest-weight")
                else:
                    M.add_edge(a, b, weight=w, merge_strategy="smallest-weight")
                    rows.append(a); cols.append(b); ws.append(w)
                    if gap:
                        MG.add_edge(a, b, weight=w, merge_strategy="smallest-weight")
            if M.num_nodes < self.n:                   # detectors no edge touches (e.g. code-capacity final round)
                M.add_boundary_edge(self.n - 1, weight=2 * wmax + 1.0)
            if gap:
                for bn in (self.n, self.n + 1):          # make sure both boundary nodes exist
                    MG.add_edge(bn, self.n + 2, weight=2 * wmax + 1.0, merge_strategy="smallest-weight")
                MG.add_boundary_edge(self.n + 2, weight=0.0)
            G = coo_matrix((ws, (rows, cols)), shape=(B + 1, B + 1)).tocsr()
            _, pred = dijkstra(G, directed=False, return_predecessors=True)
            self.matchers.append(M); self.pred.append(pred); self.gap_matchers.append(MG)
        # w with w*c + q == 0 mod d (first solution)
        self.solve = {}
        for cc in range(self.d):
            for q in range(self.d):
                sol = [w for w in range(self.d) if (w * cc + q) % self.d == 0]
                self.solve[(cc, q)] = sol[0] if sol else None

    def _path(self, pred, u, v):
        path = [v]
        while path[-1] != u:
            path.append(pred[u, path[-1]])
            if path[-1] < 0:
                raise RuntimeError("no path")
        return path[::-1]

    def _step(self, a, b, q, obs, d):
        """Traverse edge a -> b (b may be a boundary node >= n) carrying charge q to cancel at a; returns (w, new q, obs)."""
        B = self.n
        if b >= B:
            key = (a, -1)
        elif a >= B:
            key = (b, -1)
        else:
            key = (a, b) if a < b else (b, a)
        ca, cb, o = self.edges[key]["pat"]
        if key[0] != a and a < B:
            ca, cb = cb, ca
        if a >= B:                                   # leaving a boundary: unit fault on the boundary edge
            if self.trace is not None:
                self.trace.append((key, 1))
            return 1, ca % d, obs + o
        w = self.solve[(ca % d, q % d)]
        if w is None:
            return None, q, obs
        if self.trace is not None:
            self.trace.append((key, w))
        return w, (w * cb) % d, obs + w * o

    def _walk_edges(self, lev, edge_arr, binsyn, sigma, obs):
        """Charge walk over an explicit edge set (pymatching decode_to_edges_array of the split-boundary graph).
        Trails start at defects and end at a defect or a boundary node; leftover boundary-to-boundary trails (the
        forced logical strings) are walked with a unit charge."""
        d, n = self.d, self.n
        adj = {}
        for a, b in edge_arr:
            a = int(a); b = int(b)
            if b == -1:
                b = n + 2
            adj.setdefault(a, []).append([b, False]); adj.setdefault(b, []).append([a, False])

        def take(cur):
            for e in adj.get(cur, []):
                if not e[1]:
                    e[1] = True
                    for e2 in adj[e[0]]:
                        if e2[0] == cur and not e2[1]:
                            e2[1] = True
                            break
                    return e[0]
            return None

        def trail(start, q):
            nonlocal obs
            cur = start
            while True:
                nxt = take(cur)
                if nxt is None:
                    return cur, q
                w, q, obs = self._step(cur, nxt, q, obs, d)
                if w is None:
                    return None, q
                cur = nxt
                if cur >= n:
                    return cur, q
        for u in np.nonzero(binsyn)[0]:
            u = int(u)
            if not any(not e[1] for e in adj.get(u, [])):
                continue
            end, q = trail(u, int(sigma[u]))
            if end is None:
                return obs, sigma, False
            sigma[u] = 0
            if end < n:
                sigma[end] = (sigma[end] + q) % d
        for bnode in (n, n + 1):
            while any(not e[1] for e in adj.get(bnode, [])):
                end, q = trail(bnode, 0)
                if end is None:
                    return obs, sigma, False
                if end < n:
                    sigma[end] = (sigma[end] + q) % d
        return obs, sigma, True

    def decode(self, sigma: np.ndarray):
        """sigma: detector charges in Z_d (one shot).  Returns (observable correction mod d, leftover flag)."""
        d, B = self.d, self.n
        sigma = sigma.copy()
        obs = 0
        gap_min = math.inf
        for lev in range(self.k):
            binsyn = ((sigma >> lev) & 1).astype(np.uint8)
            if not binsyn.any():
                continue
            if self.gap:
                # force the parity of matchings to boundary node B_1 (node n+1): the two logical classes
                MG = self.gap_matchers[lev]
                par = int(binsyn.sum() % 2)
                sols = []
                for s1 in (0, 1):
                    s0 = (par + s1) % 2
                    syn = np.concatenate([binsyn, [s0, s1, 0]]).astype(np.uint8)
                    _, w = MG.decode(syn, return_weight=True)
                    sols.append((w, syn))
                sols.sort(key=lambda t: t[0])
                gap_min = min(gap_min, sols[1][0] - sols[0][0])
                # walk the charges along the edge set the matcher chose (its strings may end on either boundary)
                obs, sigma, ok = self._walk_edges(lev, MG.decode_to_edges_array(sols[0][1]), binsyn, sigma, obs)
                if not ok:
                    return obs % d, 0.0
                continue
            pairs = self.matchers[lev].decode_to_matched_dets_array(binsyn)
            pred = self.pred[lev]
            for u, v in pairs:
                u = int(u); v = int(v)
                tgt = B if v == -1 else v
                q = int(sigma[u])
                path = self._path(pred, u, tgt)
                for a, b in zip(path[:-1], path[1:]):
                    key = (a, b) if a < b else (b, a)
                    if b == B:
                        key = (a, -1)
                    elif a == B:
                        key = (b, -1)
                    ca, cb, o = self.edges[key]["pat"]
                    if key[0] != a:                   # traversing from the 'b' end
                        ca, cb = cb, ca
                    w = self.solve[(ca % d, q % d)]
                    if w is None:
                        return obs % d, True
                    if self.trace is not None:
                        self.trace.append((key, w))
                    obs += w * o
                    q = (w * cb) % d
                sigma[u] = 0
                if v != -1:
                    sigma[v] = (sigma[v] + q) % d
        if self.gap:
            return obs % d, (gap_min if gap_min < math.inf else 0.0)
        return obs % d, bool(sigma.any())


# ----------------------------------------------------------------------------------------------------------------
# experiment


P_BUILD = 0.01     # the circuit is built at this p; FrameSim.run(scale) rescales every location


class Experiment:
    """Circuit + detector-error model for one (D, d, basis); run() simulates and decodes at any p."""

    def __init__(self, D: int, d: int, basis: str = "Z", rounds: int | None = None, idle: bool = True,
                 code_capacity: bool = False):
        self.D, self.d, self.basis = D, d, basis
        self.rounds = (1 if code_capacity else D) if rounds is None else rounds
        self.mem = memory_circuit(D, d, self.rounds, P_BUILD, basis, idle, code_capacity)
        bad = check_determinism(self.mem.circ, 64, 0)
        assert bad == [], f"non-deterministic detectors {bad}"
        self.sector = "X" if basis == "Z" else "Z"
        self.dem = build_dem(self.mem, self.sector)
        self.sel = np.nonzero(self.mem.det_sector == self.sector)[0]
        self.locations = sum(len(qs) for _, _, qs in self.mem.noise_locs)

    def run(self, p: float, shots: int = 10_000, seed: int = 0, chunk: int = 10_000) -> dict:
        dec = PeelDecoder(self.dem, p)
        rng = np.random.default_rng(seed)
        d = self.d
        fails = leftovers = done = 0
        while done < shots:
            n = min(chunk, shots - done)
            sim = FrameSim(self.mem.circ, n, rng).run(scale=p / P_BUILD)
            Dv = sim.detector_values()[:, self.sel]
            O = np.zeros(n, dtype=np.int64)
            for i, cf in self.mem.obs:
                O += cf * sim.meas[:, i]
            O %= d
            for s in range(n):
                corr, left = dec.decode(Dv[s])
                leftovers += left
                fails += ((O[s] + corr) % d) != 0
            done += n
        return dict(D=self.D, d=d, p=p, rounds=self.rounds, basis=self.basis, shots=shots, fails=int(fails),
                    p_fail=fails / shots, leftover=int(leftovers), n_det=int(len(self.sel)), edges=len(self.dem.edges),
                    hyper_share=self.dem.hyper_prob / self.dem.total_prob, conflicts=self.dem.conflicts,
                    undetectable_share=self.dem.undetectable_prob / self.dem.total_prob, locations=self.locations)


def logical_error_rate(D: int, d: int, p: float, rounds: int | None = None, shots: int = 10_000, basis: str = "Z",
                       seed: int = 0, idle: bool = True) -> dict:
    return Experiment(D, d, basis, rounds, idle).run(p, shots, seed)


if __name__ == "__main__":
    import sys
    D = int(sys.argv[1]); d = int(sys.argv[2]); p = float(sys.argv[3]); shots = int(sys.argv[4]) if len(sys.argv) > 4 else 2000
    basis = sys.argv[5] if len(sys.argv) > 5 else "Z"
    print(verify_code(D, d))
    print(logical_error_rate(D, d, p, shots=shots, basis=basis))
