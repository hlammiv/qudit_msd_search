"""Ring generalized triorthogonality (RGT): does a Z_d stabilizer code HOST a transversal diagonal
gate layer, and what logical gate does it induce?

Gate layer.  Phi(x) = sum_i q_i(x_i) + sum_{(i,j)} T_ij(x_i, x_j) + sum_{(i,j,k)} T_ijk(x_i, x_j, x_k)
(mod M = 4d), each component a level-(<=3) phase table (see phase_tables).  The physical gate is
G = diag(zeta_M^{Phi}).  Transversal single-qudit layers (q_i only), ring-CS / ring-CCZ layers (mod-d
tables i^{c x y}, i^{c x y z}), ARBITRARY 2-body tables (e.g. the entangling A2-escapers found in the
m=2 catalog) and across-block pairs/triples on m copies of a base code are all instances.

Hosting criterion (Phase 1 M1 -> Phase 5, generalized to tables).  G X^a G^dag = diag(zeta^{Delta_a Phi})
X^a, so G preserves the codespace of a code whose X-stabilizer parts are proj_X(S) iff
    Delta_a Phi(x) == 0 (mod M)   for all a in proj_X(S) and all x in the relevant point set:
  * mode 'support' (strict, sufficient): x in supp(C) = (pure-Z stabilizers)^perp;
  * mode 'cosets' (exact when the logical quotient is free and has no soft directions): x in the
    logical-representative cosets x_p + M_X, x_p = sum_j p_j a_j over the genuine logical-X axes.
Since Phi is level-3, Delta_a Phi is a CLIFFORD-type table (constant second differences), so its
vanishing on a coset base + <s_1..s_r> is decided at the points base + {0, s, 2s, s+s'} (generalizes
cubic_solve.supp_spanning from monomials to tables); and a generator set of proj_X(S) suffices because
Delta_{a+a'}Phi(x) = Delta_a Phi(x+a') + Delta_{a'} Phi(x) with x+a' in the same module.

Induced logical.  For a hosted layer, F(p) = Phi(x_p) - Phi(0) on Z_d^k is the logical diagonal gate's
table, classified by phase_tables.classify (T / CS / CCZ content, A2 escape = universality, odd third
difference = genuine order-d cubic content).

Fault tolerance bookkeeping.  `degree_report` gives, per qudit, the number of mixed (multi-qudit)
physical terms it participates in and whether every mixed term has its legs in pairwise distinct code
blocks; a layer counts as a PROTOCOL (transversal / across-block) iff every qudit is in at most one
mixed term and all mixed terms are cross-block.

Explicit Z_4 "triorthogonality" screens for a UNIFORM single-qudit table q on X-generators s (the
x = 0 coset conditions; necessary for hosting, cheap, and the sampler's filter):
    S1  sum_i q(s_i) == 0 (mod M)                       (q-weight; carrier: #1 + 4#2 + 13#3 == 0 mod 16)
    S2  sum_i B_q(s_i, s'_i) == 0,  B_q(a,b) = q(a+b) - q(a) - q(b) + q(0)
    S3  sum_i C_q(s_i, s'_i, s''_i) == 0,  C_q = mixed third difference of q at 0
        (for the carrier C_q(a,b,c) = 4abc  =>  sum_i s_i s'_i s''_i == 0 mod 4: ring triorthogonality)
S1-S3 over all generator (multi)sets  <=>  Phi constant on <x_gens> (the x = 0 coset).
"""

from __future__ import annotations

import itertools as it
from dataclasses import dataclass, field

import numpy as np

from . import phase_tables as PT
from .phase_tables import modulus
from .symplectic import Z4Code


# --------------------------------------------------------------------------- the gate layer
@dataclass
class TransversalLayer:
    n: int
    d: int = 4
    qudit_tables: list = field(default_factory=list)      # len n; None or array (d,)
    pair_tables: list = field(default_factory=list)       # (i, j, array (d,d))
    triple_tables: list = field(default_factory=list)     # (i, j, k, array (d,d,d))
    label: str = ""

    def __post_init__(self):
        if not self.qudit_tables:
            self.qudit_tables = [None] * self.n
        assert len(self.qudit_tables) == self.n
        M = modulus(self.d)
        self.qudit_tables = [None if q is None else (np.asarray(q, dtype=np.int64) % M) for q in self.qudit_tables]
        self.pair_tables = [(i, j, np.asarray(T, dtype=np.int64) % M) for (i, j, T) in self.pair_tables]
        self.triple_tables = [(i, j, k, np.asarray(T, dtype=np.int64) % M) for (i, j, k, T) in self.triple_tables]

    # ---- constructors
    @classmethod
    def uniform(cls, n: int, q, d: int = 4, label: str = "") -> "TransversalLayer":
        return cls(n, d, [np.asarray(q)] * n, [], [], label or "uniform")

    @classmethod
    def per_qudit(cls, tables, d: int = 4, label: str = "") -> "TransversalLayer":
        return cls(len(tables), d, list(tables), [], [], label)

    @classmethod
    def cs_layer(cls, n: int, pairs, d: int = 4, label: str = "") -> "TransversalLayer":
        """ring-CS terms i^{c x_i x_j}: pairs = [(i, j, c)]."""
        grid = np.indices((d, d))
        pt = [(i, j, (4 * ((c * grid[0] * grid[1]) % d)) % modulus(d)) for (i, j, c) in pairs]
        return cls(n, d, [None] * n, pt, [], label or "cs")

    @classmethod
    def ccz_layer(cls, n: int, triples, d: int = 4, label: str = "") -> "TransversalLayer":
        """ring-CCZ terms i^{c x_i x_j x_k}: triples = [(i, j, k, c)]."""
        grid = np.indices((d, d, d))
        tt = [(i, j, k, (4 * ((c * grid[0] * grid[1] * grid[2]) % d)) % modulus(d)) for (i, j, k, c) in triples]
        return cls(n, d, [None] * n, [], tt, label or "ccz")

    @classmethod
    def across_block_pairs(cls, n0: int, T2, d: int = 4, label: str = "") -> "TransversalLayer":
        """On 2 blocks of n0 qudits: the 2-body table T2 between qudit i of block 0 and qudit i of block 1."""
        T2 = np.asarray(T2)
        return cls(2 * n0, d, [None] * (2 * n0), [(i, n0 + i, T2) for i in range(n0)], [], label or "across2")

    @classmethod
    def across_block_ccz(cls, n0: int, c: int = 1, d: int = 4, label: str = "") -> "TransversalLayer":
        """On 3 blocks: ring-CCZ i^{c x y z} between corresponding qudits (Haah-Hastings geometry)."""
        return cls.ccz_layer(3 * n0, [(i, n0 + i, 2 * n0 + i, c) for i in range(n0)], d, label or "across_ccz")

    def plus(self, other: "TransversalLayer") -> "TransversalLayer":
        assert self.n == other.n and self.d == other.d
        M = modulus(self.d)
        qt = []
        for a, b in zip(self.qudit_tables, other.qudit_tables):
            if a is None:
                qt.append(b)
            elif b is None:
                qt.append(a)
            else:
                qt.append((a + b) % M)
        return TransversalLayer(self.n, self.d, qt, self.pair_tables + other.pair_tables,
                                self.triple_tables + other.triple_tables, f"{self.label}+{other.label}")

    # ---- evaluation
    def phi(self, x) -> int:
        return int(self.phi_batch(np.asarray(x, dtype=np.int64)[None, :])[0])

    def phi_batch(self, X) -> np.ndarray:
        """Phi at each row of X (shape (N, n)), mod M."""
        X = np.asarray(X, dtype=np.int64) % self.d
        M = modulus(self.d)
        out = np.zeros(X.shape[0], dtype=np.int64)
        for i, q in enumerate(self.qudit_tables):
            if q is not None:
                out += q[X[:, i]]
        for (i, j, T) in self.pair_tables:
            out += T[X[:, i], X[:, j]]
        for (i, j, k, T) in self.triple_tables:
            out += T[X[:, i], X[:, j], X[:, k]]
        return out % M

    def delta_batch(self, a, X) -> np.ndarray:
        X = np.asarray(X, dtype=np.int64)
        a = np.asarray(a, dtype=np.int64)
        return (self.phi_batch((X + a) % self.d) - self.phi_batch(X)) % modulus(self.d)

    # ---- structure
    def component_levels(self) -> int:
        """Max hierarchy level over the component tables (the layer's level is <= this)."""
        lv = 0
        for q in self.qudit_tables:
            if q is not None:
                lv = max(lv, PT.level(q, self.d))
        for (_, _, T) in self.pair_tables:
            lv = max(lv, PT.level(T, self.d))
        for (_, _, _, T) in self.triple_tables:
            lv = max(lv, PT.level(T, self.d))
        return lv

    def mixed_terms(self) -> list[tuple[int, ...]]:
        return [(i, j) for (i, j, _) in self.pair_tables] + [(i, j, k) for (i, j, k, _) in self.triple_tables]

    def degree_report(self, block_of=None) -> dict:
        """Per-qudit count of mixed terms, max degree, cross-block flag and the PROTOCOL flag
        (every qudit in <= 1 mixed term and every mixed term has legs in pairwise distinct blocks).
        block_of: list mapping qudit -> block id (default: one block)."""
        if block_of is None:
            block_of = [0] * self.n
        deg = np.zeros(self.n, dtype=int)
        cross = True
        for term in self.mixed_terms():
            for q in term:
                deg[q] += 1
            blocks = [block_of[q] for q in term]
            if len(set(blocks)) != len(blocks):
                cross = False
        maxdeg = int(deg.max()) if self.n else 0
        return {"max_mixed_degree": maxdeg, "all_mixed_cross_block": cross,
                "n_mixed_terms": len(self.mixed_terms()),
                "protocol": (maxdeg <= 1) and cross}

    def escapes_physically(self) -> bool:
        """Does the PHYSICAL layer escape all sectors (some component table does)?"""
        from .a2 import escapes_all_sectors
        for q in self.qudit_tables:
            if q is not None and escapes_all_sectors(q, self.d):
                return True
        for (_, _, T) in self.pair_tables:
            if escapes_all_sectors(T, self.d):
                return True
        for (_, _, _, T) in self.triple_tables:
            if escapes_all_sectors(T, self.d):
                return True
        return False


# --------------------------------------------------------------------------- spanning points
def spanning_points(base, gens, d: int = 4) -> np.ndarray:
    """base + {0, s, 2s, s+s'} over the generator rows: enough to decide whether a Clifford-type table
    (constant second differences) vanishes on the whole coset base + <gens>."""
    base = np.asarray(base, dtype=np.int64)
    G = [np.asarray(g, dtype=np.int64) for g in gens]
    pts = [base.copy()]
    for g in G:
        pts.append(base + g)
        pts.append(base + 2 * g)
    for i in range(len(G)):
        for j in range(i + 1, len(G)):
            pts.append(base + G[i] + G[j])
    return np.array(pts, dtype=np.int64) % d


def cubic_spanning_points(gens, d: int = 4) -> np.ndarray:
    """All sum c_i s_i with <= 3 nonzero c_i in {1..d-1}: decides constancy of a LEVEL-3 table on <gens>."""
    G = [np.asarray(g, dtype=np.int64) for g in gens]
    n = len(G[0]) if G else 0
    pts = [np.zeros(n, dtype=np.int64)]
    r = len(G)
    for size in (1, 2, 3):
        for idx in it.combinations(range(r), size):
            for cs in it.product(range(1, d), repeat=size):
                v = np.zeros(n, dtype=np.int64)
                for c, i in zip(cs, idx):
                    v += c * G[i]
                pts.append(v)
    return np.array(pts, dtype=np.int64) % d


# --------------------------------------------------------------------------- hosting
@dataclass
class HostingVerdict:
    hosted: bool
    mode: str
    exact: bool
    reason: str = ""
    witness: dict | None = None

    def __bool__(self):
        return self.hosted


def delta_vanishes_on_coset(layer: TransversalLayer, a, base, gens) -> tuple[bool, dict | None]:
    pts = spanning_points(base, gens, layer.d) if len(gens) else np.asarray(base, dtype=np.int64)[None, :] % layer.d
    D = layer.delta_batch(a, pts)
    bad = np.nonzero(D)[0]
    if len(bad):
        i = int(bad[0])
        return False, {"a": [int(t) for t in a], "x": [int(t) for t in pts[i]], "residue": int(D[i])}
    return True, None


def phi_constant_on_module(layer: TransversalLayer, gens) -> bool:
    if not len(gens):
        return True
    pts = cubic_spanning_points(gens, layer.d)
    v = layer.phi_batch(pts)
    return bool(np.all(v == v[0]))


def logical_rep_points(code: Z4Code, max_nonzero: int = 2) -> list[np.ndarray]:
    """x_p = sum_j p_j a_j for p with <= max_nonzero nonzero entries (values 1..d-1)."""
    axes = [np.asarray(a, dtype=np.int64) for a in code.logical_axes]
    k = len(axes)
    d = code.d
    pts = [np.zeros(code.n, dtype=np.int64)]
    for size in range(1, min(max_nonzero, k) + 1):
        for idx in it.combinations(range(k), size):
            for cs in it.product(range(1, d), repeat=size):
                v = np.zeros(code.n, dtype=np.int64)
                for c, j in zip(cs, idx):
                    v += c * axes[j]
                pts.append(v % d)
    return pts


def hosts(code: Z4Code, layer: TransversalLayer, mode: str = "support") -> HostingVerdict:
    """Does diag(zeta^Phi) preserve the codespace?  mode 'support' (strict/sufficient) or 'cosets'
    (exact iff code.free_quotient_ok and code.soft == 0)."""
    assert layer.n == code.n and layer.d == code.d
    xg = code.x_gens
    if mode == "support":
        base_sets = [(np.zeros(code.n, dtype=np.int64), code.support_gens)]
        exact = code.is_css and code.free_quotient_ok and code.soft == 0
    elif mode == "cosets":
        base_sets = [(x, xg) for x in logical_rep_points(code)]
        exact = code.free_quotient_ok and code.soft == 0
    else:
        raise ValueError(mode)
    for a in xg:
        for base, gens in base_sets:
            ok, wit = delta_vanishes_on_coset(layer, a, base, gens)
            if not ok:
                return HostingVerdict(False, mode, exact, "Delta_a Phi != 0", wit)
    return HostingVerdict(True, mode, exact, "all Delta_a Phi vanish on the test cosets")


# --------------------------------------------------------------------------- screens (uniform q)
def q_weight(q, s, d: int = 4) -> int:
    q = np.asarray(q, dtype=np.int64)
    return int(q[np.asarray(s, dtype=np.int64) % d].sum() % modulus(d))


def pair_defect(q, s, t, d: int = 4) -> int:
    q = np.asarray(q, dtype=np.int64); s = np.asarray(s) % d; t = np.asarray(t) % d
    return int((q[(s + t) % d] - q[s] - q[t] + q[0]).sum() % modulus(d))


def triple_defect(q, s, t, u, d: int = 4) -> int:
    q = np.asarray(q, dtype=np.int64); s = np.asarray(s) % d; t = np.asarray(t) % d; u = np.asarray(u) % d
    C = (q[(s + t + u) % d] - q[(s + t) % d] - q[(s + u) % d] - q[(t + u) % d] + q[s] + q[t] + q[u] - q[0])
    return int(C.sum() % modulus(d))


def uniform_screens(code: Z4Code, q) -> dict:
    """S1-S3 on the X-generators for a uniform single-qudit table q (necessary for hosting)."""
    xg = code.x_gens
    d = code.d
    s1 = all(q_weight(q, s, d) == 0 for s in xg)
    s2 = all(pair_defect(q, xg[i], xg[j], d) == 0 for i in range(len(xg)) for j in range(i, len(xg)))
    s3 = all(triple_defect(q, xg[i], xg[j], xg[k], d) == 0
             for i in range(len(xg)) for j in range(i, len(xg)) for k in range(j, len(xg)))
    return {"S1_qweight": s1, "S2_pair": s2, "S3_triple": s3, "all": s1 and s2 and s3}


# --------------------------------------------------------------------------- induced logical
def induced_logical(code: Z4Code, layer: TransversalLayer, k_max: int = 6) -> np.ndarray:
    """F(p) = Phi(x_p) - Phi(0) over Z_d^k, p over the genuine logical axes (k <= k_max)."""
    axes = np.array(code.logical_axes, dtype=np.int64)
    k = len(axes)
    if k > k_max:
        raise ValueError(f"k={k} too large for a dense logical table (k_max={k_max})")
    d = code.d
    if k == 0:
        return np.zeros((), dtype=np.int64)
    P = np.array(list(it.product(range(d), repeat=k)), dtype=np.int64)
    X = (P @ axes) % d
    F = layer.phi_batch(X)
    return PT.normalize(F.reshape((d,) * k), d)


def induced_logical_report(code: Z4Code, layer: TransversalLayer, k_max: int = 6) -> dict:
    F = induced_logical(code, layer, k_max)
    rep = PT.classify(F, code.d) if F.ndim else {"level": 0}
    rep["table"] = F.reshape(-1).tolist() if F.ndim else []
    return rep


# --------------------------------------------------------------------------- dense cross-check
def diagonal_leakage(code: Z4Code, layer: TransversalLayer) -> float:
    """||(I - Pi) G Pi|| with the dense signed-stabilizer projector (n <= 5 non-CSS; n <= 6 at most)."""
    from .noncss_transversal import codespace_projector
    n, d = code.n, code.d
    if d ** n > 4096:
        raise ValueError("dense leakage only for d^n <= 4096")
    Pi = codespace_projector([list(map(int, r)) for r in code.rows], d, n)
    if abs(np.trace(Pi).real) < 0.5:
        return float("nan")                       # empty codespace
    X = np.array(list(it.product(range(d), repeat=n)), dtype=np.int64)
    ph = np.exp(2j * np.pi * layer.phi_batch(X) / modulus(d))
    GPi = ph[:, None] * Pi                        # G Pi (G diagonal)
    return float(np.linalg.norm(GPi - Pi @ GPi))  # ||(I - Pi) G Pi||_F


# --------------------------------------------------------------------------- logical escape (transmission)
def axis_parity_report(code: Z4Code) -> dict:
    """For a UNIFORM single-qudit level-3 layer q with Delta^3 q = 4 (any carrier-type table), the
    induced logical F satisfies
        Delta_{2 e_i} Delta_{2 e_l} F  =  8 * #{t : a_{i,t}, a_{l,t} both odd}   (mod 16)
        Delta_{e_i}Delta_{e_j}Delta_{e_l} F = 4 * sum_t a_{i,t} a_{j,t} a_{l,t}      (mod 16)
    over the logical-X axes a_i.  Hence F ESCAPES A2 (universality transmitted) iff some axis has an
    odd number of odd entries or some pair of axes has an odd odd-overlap; and F has genuine order-4
    cubic content iff some triple (with repetition) has an odd all-odd overlap.  Pure code property."""
    axes = [np.asarray(a, dtype=np.int64) % code.d for a in code.logical_axes]
    odd = [(a % 2).astype(int) for a in axes]
    k = len(axes)
    odd_weights = [int(o.sum()) for o in odd]
    pair_overlaps = {(i, l): int((odd[i] * odd[l]).sum()) for i in range(k) for l in range(i, k)}
    escape = any(v % 2 for v in pair_overlaps.values())
    triple_odd = any(int((odd[i] * odd[j] * odd[l]).sum()) % 2
                     for i in range(k) for j in range(i, k) for l in range(j, k))
    return {"k": k, "odd_weights": odd_weights, "pair_odd_overlaps": {str(key): v for key, v in pair_overlaps.items()},
            "logical_escape_possible": escape, "logical_odd_cubic_possible": triple_odd}


def layer_from_poly_mod4(coeffs, mons, n: int) -> TransversalLayer:
    """A Phase-5 cubic phase i^{sum c * monomial} (coeffs aligned with cubic.monomials(n,3)) as a layer."""
    qt = [None] * n; pairs = []; trips = []
    g1 = np.arange(4); g2 = np.indices((4, 4)); g3 = np.indices((4, 4, 4))
    for c, m in zip(coeffs, mons):
        c = int(c) % 4
        if c == 0:
            continue
        idx = sorted(set(m))
        if len(idx) == 1:
            i = idx[0]
            t = (4 * ((c * g1 ** len(m)) % 4)) % 16
            qt[i] = t if qt[i] is None else (qt[i] + t) % 16
        elif len(idx) == 2:
            i, j = idx
            e = (m.count(i), m.count(j))
            pairs.append((i, j, (4 * ((c * g2[0] ** e[0] * g2[1] ** e[1]) % 4)) % 16))
        else:
            i, j, k = idx
            trips.append((i, j, k, (4 * ((c * g3[0] * g3[1] * g3[2]) % 4)) % 16))
    return TransversalLayer(n, 4, qt, pairs, trips, "poly_mod4")


def escape_rank(code: Z4Code, layer: TransversalLayer, k_max: int = 8) -> dict:
    """Number of universality-enabling LOGICAL outputs of a hosted layer.  For the induced logical table F
    on Z_d^k (d = 4) the sector-mixing content is the symmetric F_2 bilinear form
        B[i,l] = [ Delta_{2 e_i} Delta_{2 e_l} F  !=  0 ]   (values are 0 or 8 mod 16 for a level-3 F),
    and its F_2-rank r is the number of logical qudits that receive an escaping (T-like) gate after a logical
    Clifford change of basis (diagonal part = single-qudit escapers, hyperbolic pairs = 2-body escapers, both
    universality-enabling).  The universality-enabling overhead exponent is gamma_T = log(n / r) / log(d),
    NOT log(n / k) / log(d): a code whose k logical outputs carry only r escaping gates yields r T-like states
    per round (the other k - r outputs are non-universalizing level-3 magic).  Requires k <= k_max (dense F)."""
    from .lifts import binary_rank
    F = induced_logical(code, layer, k_max)
    k = F.ndim
    if k == 0:
        return {"k": 0, "rank": 0, "B": []}
    B = np.zeros((k, k), dtype=np.int64)
    for i in range(k):
        for l in range(k):
            D = PT.iterated_delta(F, [PT.unit(k, i, 2), PT.unit(k, l, 2)], code.d) % modulus(code.d)
            B[i, l] = int(np.any(D != 0))
    r = binary_rank(B) if k else 0
    return {"k": k, "rank": r, "B": B.tolist(), "escapes": bool(r > 0)}


# --------------------------------------------------------------------------- convenience
def carrier_layers(n: int, d: int = 4, strict_all: bool = False) -> list[TransversalLayer]:
    """Uniform layers for each strict level-3 single-qudit table (12 canonical carriers, or all 48)."""
    tabs = PT.strict_level3_single(d) if strict_all else PT.CARRIERS_D4
    return [TransversalLayer.uniform(n, np.array(q), d, label=f"uniform{list(q)}") for q in tabs]
