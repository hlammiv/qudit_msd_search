"""Cultivation circuits on the signed Z_d triangular colour code (Clifford proxy for T), simulated with zdsim.

Proxy (as Gidney–Shutty–Jones do with T -> S): the escaping T is replaced by the Clifford S^sigma (diag
zeta_{2d}^{sigma x^2}); the per-qudit check T^{tau} X T^{-tau} becomes X Z^{sigma tau} and the injected state
S^sigma|+> is a stabiliser state with the same check structure (cultivation.py: the logical check equals the
physical one).  Noise: DEP1(p) after 1-qudit gates, DEP2(p) after CSUM, X^j error p after resets / before
measurements.

Protocol (stage-1 cultivation, distance D colour code):
  inject  : unitary CSS encoder, the proxy-T on the input qudit q*            (fault distance 1)
  rounds  : stabiliser round, then `checks` T-checks, ...                     (all detectors post-selected)
  check   : Z_d cat state on n ancillas, data: T^{-tau} . CSUM(cat_i -> data_i) . T^{tau}; uncompute the cat,
            measure: cat_0 (Fourier) = check eigenvalue, cat_i = flags
  output  : ideal decode of the final data frame (lookup, weight-1) and test whether the residual logical
            error anticommutes with the check (= changes the cultivated state).
"""

from __future__ import annotations

import numpy as np

from .cultivation import color_code_plaquettes, honeycomb_signs
from .zdsim import Circuit, FrameSim, check_determinism


_SLOTS = {}


def stim_slots(distance: int) -> dict:
    """(plaquette, data) -> CX slot 0..5 of stim's color_code:memory_xyz schedule (conflict-free on every data qudit),
    in the plaquette / data order of color_code_plaquettes."""
    if distance not in _SLOTS:
        import stim
        c = stim.Circuit.generated("color_code:memory_xyz", distance=distance, rounds=2).flattened()
        tick, legs, anc, data = 0, {}, set(), []
        for ins in c:
            t = [x.value for x in ins.targets_copy()]
            if ins.name == "TICK":
                tick += 1
            elif ins.name == "CX":
                for a, b in zip(t[::2], t[1::2]):
                    legs.setdefault(b, []).append((tick, a))
            elif ins.name == "MY":
                data = t
            elif ins.name == "MR":
                anc |= set(t)
        idx = {q: i for i, q in enumerate(data)}
        out = {}
        for r, a in enumerate(sorted(anc)):
            t0 = min(tt for tt, _ in legs[a])
            for tt, q in legs[a]:
                if tt - t0 < 6:
                    out[(r, idx[q])] = tt - t0
        _SLOTS[distance] = out
    return _SLOTS[distance]


def _inv(a, d):
    return pow(int(a) % d, -1, d)


def pivot_form(H, d):
    """Row-reduce H over Z_d with unit pivots: (G, pivots), G[k, pivots[j]] = delta_kj."""
    G = H.copy() % d
    piv = []
    for r in range(G.shape[0]):
        cols = [c for c in range(G.shape[1]) if c not in piv and np.gcd(int(G[r, c]), d) == 1]
        c = cols[0]
        G[r] = (G[r] * _inv(G[r, c], d)) % d
        for k in range(G.shape[0]):
            if k != r and G[k, c]:
                G[k] = (G[k] - G[k, c] * G[r]) % d
        piv.append(c)
    return G, piv


class ColorCult:
    def __init__(self, distance: int = 3, d: int = 4, p: float = 1e-3, sigma: int = 1, inj: int | None = None,
                 checks_per_round: int = 2, rounds: int = 2, noisy_injection: bool = True, parallel: bool = False,
                 schedule=None, final_by_colour: bool = False, hooks: bool = False,
                 final_block: bool = False):
        self.parallel, self.schedule = parallel or schedule is not None, schedule
        self.final_by_colour, self.hooks, self.final_block = final_by_colour, hooks, final_block
        H, xy = color_code_plaquettes(distance)
        self.H, self.tau, self.d, self.n = H, honeycomb_signs(H, xy), d, H.shape[1]
        self.distance, self._tab = distance, None
        self.sigma = sigma
        self.inj = sigma if inj is None else inj
        self.MZ = (H * self.tau) % d
        self.p = p
        c = self.c = Circuit(d)
        self.data = c.q(self.n)
        self.xa = c.q(len(H)); self.za = c.q(len(H)); self.cat = c.q(self.n)
        self.cat2 = c.q(self.n) if self.parallel else None     # second cat register: consecutive checks overlap
        self.encode(noisy_injection)
        for _ in range(rounds):
            self.stab_round()
            for _ in range(checks_per_round):
                self.check()
        self.stab_round(final=True)

    # noisy primitives
    def R(self, qs):
        self.c.R(qs); self.c.XERR(qs, self.p)
    def M(self, qs):
        self.c.XERR(qs, self.p); return self.c.M(qs)
    def F(self, qs):
        self.c.F(qs); self.c.DEP1(qs, self.p)
    def FD(self, qs):
        self.c.FD(qs); self.c.DEP1(qs, self.p)
    def S(self, qs, a):
        self.c.S(qs, a); self.c.DEP1(qs, self.p)
    def CSUM(self, pairs, w=1):
        """hooks=True: a layer whose pairs share a qudit (a fan-out) is applied pair by pair, each followed by its DEP2, so
        faults between the legs (hook errors) are sampled; the default applies the whole layer, then DEP2 on every pair."""
        pairs = list(pairs)
        if getattr(self, "hooks", False):
            qs = [q for pq in pairs for q in pq]
            if len(set(qs)) < len(qs):
                ws = np.broadcast_to(np.asarray(w, dtype=np.int64), (len(pairs),))
                for pq, wq in zip(pairs, ws):
                    self.c.CSUM([pq], int(wq)); self.c.DEP2([pq], self.p)
                return
        self.c.CSUM(pairs, w); self.c.DEP2(pairs, self.p)

    def encode(self, noisy: bool):
        d, D = self.d, self.data
        p0 = self.p
        if not noisy:
            self.p = 0.0
        G, piv = pivot_form(self.H, d)
        XL = np.ones(self.n, dtype=np.int64)
        for k, pc in enumerate(piv):
            XL = (XL - XL[pc] * G[k]) % d
        qs = [i for i in range(self.n) if i not in piv and np.gcd(int(XL[i]), d) == 1][0]
        XL = (XL * _inv(XL[qs], d)) % d
        self.qstar = qs
        self.R(D)
        self.F([D[qs]]); self.S([D[qs]], self.inj)                   # proxy |T> = S^inj |+> on q*
        self.CSUM([(D[qs], D[i]) for i in range(self.n) if i != qs and XL[i]], [int(XL[i]) for i in range(self.n) if i != qs and XL[i]])
        for k, pc in enumerate(piv):
            self.F([D[pc]])
            tg = [i for i in range(self.n) if i != pc and G[k, i]]
            self.CSUM([(D[pc], D[i]) for i in tg], [int(G[k, i]) for i in tg])
        self.p = p0

    def _measure_plaquettes(self):
        """parallel schedule: all X plaquettes in one block, then all Z plaquettes (CSUMs into a common target, resp. out of a
        common control, commute, so each block measures its stabilisers exactly in any step order); the steps are greedy:
        each plaquette visits its data in the serial order, a leg waits only while its data qudit is busy.  Returns (m_x, m_z) per plaquette."""
        D, H, d = self.data, self.H, self.d
        xa, za = self.xa[:len(H)], self.za[:len(H)]
        if getattr(self, "schedule", None) is not None:
            return self._measure_plaquettes_slots(*self.schedule)
        steps = []
        busy_q = {}
        for r in range(len(H)):
            k = 0
            for i in np.nonzero(H[r])[0]:          # each plaquette keeps the serial leg order (its hooks stay the serial ones)
                while k in busy_q.setdefault(int(i), set()):
                    k += 1
                busy_q[int(i)].add(k)
                while len(steps) <= k:
                    steps.append([])
                steps[k].append((r, int(i)))
                k += 1
        self.R(xa + za); self.F(xa)
        for st in steps:
            self.CSUM([(xa[r], D[i]) for r, i in st], 1)
        for st in steps:
            self.CSUM([(D[i], za[r]) for r, i in st], [int(self.tau[i]) % d for r, i in st])
        self.FD(xa)
        m = self.M(xa + za)
        return m[:len(H)], m[len(H):]

    def _measure_plaquettes_slots(self, ax: int, bx: int, az: int, bz: int):
        """stim-slot schedule: the X leg of plaquette r on data q at tick ax s + bx, its Z leg at az s + bz (s = stim slot
        0..5); legs are emitted tick by tick (X before Z within a tick), so residual collisions are serialised by the
        ASAP scheduler; correctness is checked by detector determinism."""
        D, H, d = self.data, self.H, self.d
        xa, za = self.xa[:len(H)], self.za[:len(H)]
        slot = stim_slots(self.distance)
        ev = []
        for r in range(len(H)):
            for i in np.nonzero(H[r])[0]:
                sl = slot[(r, int(i))]
                ev.append((ax * sl + bx, 0, r, int(i))); ev.append((az * sl + bz, 1, r, int(i)))
        self.R(xa + za); self.F(xa)
        for t in sorted({e[0] for e in ev}):
            for kind in (0, 1):
                legs = [(r, i) for tt, k, r, i in ev if tt == t and k == kind]
                if not legs:
                    continue
                if kind == 0:
                    self.CSUM([(xa[r], D[i]) for r, i in legs], 1)
                else:
                    self.CSUM([(D[i], za[r]) for r, i in legs], [int(self.tau[i]) % d for r, i in legs])
        self.FD(xa)
        m = self.M(xa + za)
        return m[:len(H)], m[len(H):]

    def _plaquette_colours(self):
        """proper 3-colouring of the plaquette overlap graph (plaquettes of one colour share no data)."""
        H = self.H
        adj = [[r2 for r2 in range(len(H)) if r2 != r and np.any(H[r] * H[r2])] for r in range(len(H))]
        col = [-1] * len(H)

        def bt(r):
            if r == len(H):
                return True
            for c in range(3):
                if all(col[x] != c for x in adj[r]):
                    col[r] = c
                    if bt(r + 1):
                        return True
            col[r] = -1
            return False
        assert bt(0), "plaquettes not 3-colourable"
        return col

    def _measure_plaquettes_by_colour(self):
        """final-round schedule: one colour class at a time (X legs in stim slot order, then Z legs); the X checks of later
        classes catch the Z hooks of earlier ones as in the serial round, and the classes overlap under ASAP."""
        D, H, d = self.data, self.H, self.d
        xa, za = self.xa[:len(H)], self.za[:len(H)]
        slot = stim_slots(self.distance)
        col = self._plaquette_colours()
        m_x, m_z = [None] * len(H), [None] * len(H)
        for c in range(3):
            P = [r for r in range(len(H)) if col[r] == c]
            if not P:
                continue
            self.R([xa[r] for r in P] + [za[r] for r in P]); self.F([xa[r] for r in P])
            if self.final_by_colour == "sup":            # serial leg order (ascending data index) inside each class
                legs = [(k, r, int(i)) for r in P for k, i in enumerate(np.nonzero(H[r])[0])]
            else:
                legs = [(slot[(r, int(i))], r, int(i)) for r in P for i in np.nonzero(H[r])[0]]
            for sl in range(6):
                L = [(r, i) for s_, r, i in legs if s_ == sl]
                if L:
                    self.CSUM([(xa[r], D[i]) for r, i in L], 1)
            for sl in range(6):
                L = [(r, i) for s_, r, i in legs if s_ == sl]
                if L:
                    self.CSUM([(D[i], za[r]) for r, i in L], [int(self.tau[i]) % d for r, i in L])
            self.FD([xa[r] for r in P])
            m = self.M([xa[r] for r in P] + [za[r] for r in P])
            for j, r in enumerate(P):
                m_x[r], m_z[r] = m[j], m[len(P) + j]
        return m_x, m_z

    def _cat_tree(self, C, sign):
        """log-depth Z_d cat (sign = +1, C[0] carries the value) or its inverse (sign = -1)."""
        levels, k = [], 1
        while k < len(C):
            levels.append([(C[i], C[i + k]) for i in range(k) if i + k < len(C)]); k *= 2
        for lv in (levels if sign > 0 else levels[::-1]):
            self.CSUM(lv, sign % self.d)

    def stab_round(self, final: bool = False):
        """final: the last round before the ideal decode stays serial even in the parallel schedule -- there the X checks
        of the later plaquettes catch the Z hooks of the earlier ones, which the block schedule would leave undetected."""
        D, H, d = self.data, self.H, self.d
        if self.parallel and final and getattr(self, "final_by_colour", False):
            m_x, m_z = self._measure_plaquettes_by_colour()
            for ma, mb in zip(m_x, m_z):
                self.c.detector([(ma, 1)], "X-stab"); self.c.detector([(mb, 1)], "Z-stab")
            return
        if self.parallel and (not final or getattr(self, "final_block", False)):   # final_block: with final-syndrome postselection
            m_x, m_z = self._measure_plaquettes()
            for ma, mb in zip(m_x, m_z):
                self.c.detector([(ma, 1)], "X-stab"); self.c.detector([(mb, 1)], "Z-stab")
            return
        for r in range(len(H)):
            a, b = self.xa[r], self.za[r]
            sup = list(np.nonzero(H[r])[0])
            self.R([a, b]); self.F([a])
            self.CSUM([(a, D[i]) for i in sup], 1)
            self.CSUM([(D[i], b) for i in sup], [int(self.tau[i]) % d for i in sup])
            self.FD([a])
            ma, mb = self.M([a, b])
            self.c.detector([(ma, 1)], "X-stab"); self.c.detector([(mb, 1)], "Z-stab")

    def check(self):
        D, C, d, n = self.data, self.cat, self.d, self.n
        if self.parallel and self.cat2 is not None:     # alternate the two cat registers
            self.cat, self.cat2 = self.cat2, self.cat
        self.R(C); self.F([C[0]])
        if self.parallel:
            self._cat_tree(C, 1)
        else:
            self.CSUM([(C[0], C[i]) for i in range(1, n)], 1)
        a = [int(-self.sigma * t) % d for t in self.tau]
        self.S(D, a)
        self.CSUM([(C[i], D[i]) for i in range(n)], 1)
        self.S(D, [(-v) % d for v in a])
        if self.parallel:
            self._cat_tree(C, -1)
        else:
            self.CSUM([(C[0], C[i]) for i in range(1, n)], -1 % d)
        self.FD([C[0]])
        m = self.M(C)
        self.c.detector([(m[0], 1)], "T-check")
        for i in range(1, n):
            self.c.detector([(m[i], 1)], "flag")

    # ---------------------------------------------------------------- evaluation
    def _tables(self):
        """min-weight lookup (weight <= (D-1)/2) for X-type and Z-type errors."""
        import itertools as it
        if getattr(self, "_tab", None) is not None:
            return self._tab
        d, n = self.d, self.n
        t = (self.distance - 1) // 2
        tx, tz = {}, {}
        for w in range(1, t + 1):
            for sup in it.combinations(range(n), w):
                for vals in it.product(range(1, d), repeat=w):
                    e = np.zeros(n, dtype=np.int64); e[list(sup)] = vals
                    tx.setdefault(tuple((self.MZ @ e) % d), e)     # X-type error -> Z-stab syndrome
                    tz.setdefault(tuple((self.H @ e) % d), e)      # Z-type error -> X-stab syndrome
        self._tab = (tx, tz)
        return self._tab

    def logical_failures(self, x, z):
        """x,z: (shots, n) final data frame -> bool array: residual error changes the cultivated state."""
        d = self.d
        tx, tz = self._tables()
        zero = np.zeros(self.n, dtype=np.int64)
        sx = (x @ self.MZ.T) % d; sz = (z @ self.H.T) % d
        cx = np.array([zero if not s.any() else tx.get(tuple(s), None) if tuple(s) in tx else np.full(self.n, -1) for s in sx])
        cz = np.array([zero if not s.any() else tz.get(tuple(s), None) if tuple(s) in tz else np.full(self.n, -1) for s in sz])
        undecodable = (cx[:, 0] < 0) | (cz[:, 0] < 0)
        xr, zr = (x - np.where(cx < 0, 0, cx)) % d, (z - np.where(cz < 0, 0, cz)) % d
        chk_x = np.ones(self.n, dtype=np.int64); chk_z = (self.sigma * self.tau) % d
        sym = (xr @ chk_z - zr @ chk_x) % d
        return undecodable | (sym != 0), (sx.any(1) | sz.any(1))

    def simulate(self, shots: int, seed: int = 0, batch: int = 200_000):
        rng = np.random.default_rng(seed)
        acc = fail = fail_triv = acc_triv = 0
        done = 0
        while done < shots:
            b = min(batch, shots - done)
            s = FrameSim(self.c, b, rng).run()
            ok = ~np.any(s.detector_values(), axis=1)
            x, z = s.x[ok][:, self.data], s.z[ok][:, self.data]
            f, nontriv = self.logical_failures(x, z)
            acc += int(ok.sum()); fail += int(f.sum())
            acc_triv += int((~nontriv).sum()); fail_triv += int((f & ~nontriv).sum())
            done += b
        return {"shots": shots, "accepted": acc, "accept_rate": acc / shots, "logical_failures": fail,
                "p_L": fail / max(acc, 1), "failures_trivial_final_syndrome": fail_triv,
                "accepted_trivial_final_syndrome": acc_triv}      # final-syndrome postselection (the next stage detects)

    def single_fault_scan(self, final_postselect: bool = False):
        """Inject every single fault (each noise location, each non-identity Pauli); count undetected
        harmful ones (fault distance >= 2 iff 0).  final_postselect: a fault that leaves a nonzero final data syndrome
        counts as detected (the next stage, e.g. the escape's colour-region rounds, sees it)."""
        d, c = self.d, self.c
        faults = []
        for k, op in enumerate(c.ops):
            if op[0] in ("DEP1", "XERR", "ZERR"):
                for q in op[1]:
                    for v in range(1, d * d):
                        faults.append((k, [q], [(v % d, v // d)]))
            elif op[0] == "DEP2":
                for (a, b) in op[1]:
                    for v in range(1, d ** 4):
                        faults.append((k, [a, b], [(v % d, (v // d) % d), ((v // d ** 2) % d, v // d ** 3)]))
        bad = []
        B = 50_000
        for s0 in range(0, len(faults), B):
            chunk = faults[s0:s0 + B]
            sim = FrameSim(c, len(chunk), np.random.default_rng(1), gauge=True)
            forced = {}
            for j, (k, qs, P) in enumerate(chunk):
                forced.setdefault((k, tuple(qs)), []).append((j, P))
            # apply by op index: build per-op arrays
            per_op = {}
            for (k, qs), lst in forced.items():
                mask = np.zeros(len(chunk), dtype=bool); Pa = np.zeros((len(chunk), len(qs), 2), dtype=np.int64)
                for j, P in lst:
                    mask[j] = True; Pa[j] = np.array(P)
                per_op.setdefault(k, []).append((mask, list(qs), Pa))
            sim.run(scale=0.0, forced=per_op)
            ok = ~np.any(sim.detector_values(), axis=1)
            f, nontriv = self.logical_failures(sim.x[:, self.data], sim.z[:, self.data])
            if final_postselect:
                ok = ok & ~nontriv
            for j in np.nonzero(ok & f)[0]:
                bad.append(chunk[j])
        return {"single_faults": len(faults), "undetected_harmful": len(bad), "examples": bad[:10], "harmful_ops": sorted({k for k, _, _ in bad})}


def find_sigma_inj(distance=3, d=4):
    """(sigma, inj) choices for which every declared detector is deterministic."""
    out = []
    for sigma in range(1, d):
        for inj in range(1, d):
            cc = ColorCult(distance, d, 0.0, sigma, inj, checks_per_round=1, rounds=1)
            if not check_determinism(cc.c):
                out.append((sigma, inj))
    return out
