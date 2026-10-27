"""Escape of a cultivated state from the signed Z_d colour code into a surface code: Gidney et al.'s grafted code
(arXiv:2409.17595, d1 = 5 -> d2 = 15) translated to ring qudits.

Data (`data/graft_d5_z4.json`, produced by graft_extract.py + graft_signs_solve.py from the published stim circuit):
  * colour : the 19-qudit 6.6.6 colour code (9 X + 9 Z plaquettes), Z plaquettes signed by the honeycomb bipartition tau
             (the cultivation substrate of cultivation.py);
  * grafted: 228 checks on 229 data qudits = colour qudits + 105 fresh |0> + 105 fresh |+>; 14 colour plaquettes kept,
             boundary plaquettes extended into the surface region, interface checks of weight 2/3/5;
  * final  : the matchable code after the transition (8 colour-region Z plaquettes split into weight-2 pairs, 4 plaquettes
             dropped); 232 checks.
  Signs: solved over GF(2) so that every X/Z overlap sums to exactly zero (valid for every d = 2^k), kept/extended/
  refined checks agree with their parents up to a global unit, and the logicals Z_L = Z^tau (colour) x Z^v (14 fresh |0>
  qudits on a diagonal), X_L = X (colour) x X^u (14 fresh |+>) commute with every grafted and final check.

Circuit (zdsim): one ancilla per check, Z checks by CSUM(data -> anc), X checks by F CSUM(anc -> data) F^dag, all Z
checks of a round first, then all X checks (so interleaving never changes the measured operators); legs scheduled into
conflict-free layers, surface-region weight-4 checks with hooks perpendicular to the logical chains.  Noise as in
surface_zd (DEP1 idle during reset, XERR after reset / before measurement, DEP1 after F, DEP2 after CSUM).
Detectors: consecutive differences of identical checks; first-round singles for checks determined by the fresh
preparation; parent/child relations at the escape (extended plaquette vs colour plaquette) and at the transition
(split plaquette vs its parts); final data readout against the last round.  Postselection: detectors of the colour-region
plaquettes of one colour class per sector during the grafted rounds (Gidney's recipe; makes the region matchable).
"""
from __future__ import annotations

import itertools as it
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import surface_zd as SZ
from .surface_zd import DEM, PeelDecoder, _same_class, _val2, propagate_basis, Memory, P_BUILD
from .zdsim import Circuit, FrameSim, check_determinism

DATA = Path(__file__).parent / "data" / "graft_d5_z4.json"


@dataclass(frozen=True)
class Check:
    typ: str
    legs: tuple            # ((qudit, weight), ...) sorted by qudit
    anc_coord: tuple = ()

    @property
    def support(self):
        return frozenset(q for q, _ in self.legs)

    def weights(self):
        return dict(self.legs)


def load_graft(d: int = 4):
    J = json.load(open(DATA))
    coords = {int(k): tuple(v) for k, v in J["coords"].items()}

    def conv(lst):
        out = []
        for v in lst:
            legs = tuple(sorted((int(q), w % d) for q, w in v["legs"].items()))
            out.append(Check(v["type"], legs, coords[v["anc"]]))
        return out
    codes = {k: conv(J[k]) for k in ("colour", "grafted", "final")}
    fresh = {int(q): b for q, b in J["fresh"].items()}
    ZL = {int(q): w % d for q, w in J["ZL"].items()}
    XL = {int(q): w % d for q, w in J["XL"].items()}
    colours = {frozenset(s): c for s, c in J["colours"]}          # Z-plaquette support -> colour class
    return dict(coords=coords, codes=codes, fresh=fresh, ZL=ZL, XL=XL, colours=colours, tau={int(q): w for q, w in J["tau"].items()})


def unit_ratio(a: Check, b: Check, d: int):
    """sigma with legs(a) = sigma * legs(b) on the common support (None if not proportional)."""
    wa, wb = a.weights(), b.weights()
    common = set(wa) & set(wb)
    sig = None
    for q in common:
        for s in (1, d - 1):
            if (wa[q] - s * wb[q]) % d == 0:
                if sig is None:
                    sig = s
                elif sig != s:
                    return None
                break
        else:
            return None
    return sig


# ------------------------------------------------------------------------------------------------------------------
# scheduling

X_ORDER = [(-1, -1), (1, 1), (-1, 1), (1, -1)]      # hook (last two) perpendicular to the X_L diagonal (1, 1)
Z_ORDER = [(-1, 1), (1, -1), (-1, -1), (1, 1)]      # hook perpendicular to the Z_L diagonal (1, -1)


def schedule(checks, coords):
    """Assign each (check, leg) a layer: surface weight-4 checks by the fixed order, the rest greedily; no data qudit is
    touched twice in a layer.  Returns list of layers, each a list of (check index, qudit, weight)."""
    layers = []
    busy = []                                   # per layer: set of qudits

    def place(ci, q, w, lay):
        while len(layers) <= lay:
            layers.append([]); busy.append(set())
        if q in busy[lay]:
            return False
        layers[lay].append((ci, q, w)); busy[lay].add(q)
        return True
    pending = []
    for ci, ch in enumerate(checks):
        order = X_ORDER if ch.typ == "X" else Z_ORDER
        ax, ay = ch.anc_coord
        rel = {q: (coords[q][0] - ax, coords[q][1] - ay) for q, _ in ch.legs}
        if len(ch.legs) == 4 and all(r in order for r in rel.values()):
            for q, w in ch.legs:
                if not place(ci, q, w, order.index(rel[q])):
                    pending.append((ci, q, w))
        else:
            pending.extend((ci, q, w) for q, w in ch.legs)
    for ci, q, w in pending:
        lay = 0
        used = {l for l, L in enumerate(layers) for (c2, _, _) in L if c2 == ci}
        while lay in used or not place(ci, q, w, lay):
            lay += 1
    return layers


# ------------------------------------------------------------------------------------------------------------------
# circuit builder

@dataclass
class Phase:
    code: list                          # checks
    rounds: int
    noisy: bool = True
    fresh: dict = field(default_factory=dict)     # qudit -> 'R' | 'RX' prepared at phase start
    relations: list = field(default_factory=list)  # (parent Check from previous code, [(child index in code, sigma)])


def build(d: int, data: list, phases: list, basis: str, logical: dict, p: float, idle: bool = True,
          postselect=None, flag_weight: int = 0, flag_phases=None):
    """phases: list of Phase.  The first phase is the (noiseless) preparation of `data` in `basis` + one round of its code.
    logical: qudit -> weight of the logical measured at the end (Z_L for basis 'Z', X_L for 'X').
    postselect: function(check, phase_index, round) -> bool marking detectors to postselect.
    flag_weight: if > 0, every check of weight >= flag_weight in a noisy phase gets a flag qudit: two cancelling CSUMs
    between its ancilla and the flag, right after the check's first data CSUM and right before its last, so an ancilla
    fault that would spread to >= 2 data qudits (a hook) changes the flag; the flag measurement is a postselected
    detector.  X-type check (ancilla = CSUM control): hooks are X^a on the ancilla -> flag |0>, CSUM(anc->flag, +1/-1),
    measured in Z.  Z-type check (ancilla = CSUM target): hooks are Z^b -> flag |+>, CSUM(flag->anc, +1/-1), F^-1, Z."""
    c = Circuit(d)
    qidx = {q: i for i, q in enumerate(data)}
    dq = c.q(len(data))
    coords = None
    det_sector, det_post, det_tags = [], [], []
    noise_locs = []
    prepared = set()

    def noise(kind, qs, pr):
        if pr > 0:
            getattr(c, kind)(qs, pr)
            noise_locs.append((len(c.ops) - 1, kind, list(qs)))

    def det(terms, sector, post, tag):
        c.detector(terms, tag); det_sector.append(sector); det_post.append(post); det_tags.append(tag)

    prev = {}            # Check -> measurement index of its last measurement
    fresh_set, fresh_basis = set(), {}     # qudits prepared since the last measured round, and their basis
    for pi, ph in enumerate(phases):
        pr = p if ph.noisy else 0.0
        checks = ph.code
        anc = c.q(len(checks))
        xanc = [anc[i] for i, ch in enumerate(checks) if ch.typ == "X"]
        layers = schedule(checks, ph.coords) if hasattr(ph, "coords") else None
        # flags: check index -> flag qudit; couplings slotted between the layers of the check's own type
        flagged = [i for i, ch in enumerate(checks) if flag_weight and ph.noisy and len(ch.legs) >= flag_weight
                   and (flag_phases is None or pi in flag_phases)]
        flag_of = dict(zip(flagged, c.q(len(flagged))))
        zflags = [flag_of[i] for i in flagged if checks[i].typ == "Z"]
        early, late = {"Z": {}, "X": {}}, {"Z": {}, "X": {}}       # typ -> layer index -> [(control, target, w)]
        for i in flagged:
            typ = checks[i].typ
            lays = sorted(l for l, L in enumerate(LAYERS[pi][typ]) if any(ci == i for ci, _, _ in L))
            a, f = anc[i], flag_of[i]
            pair = (a, f) if typ == "X" else (f, a)
            early[typ].setdefault(lays[0], []).append((pair, 1))
            late[typ].setdefault(lays[-1], []).append((pair, -1))
        # preparation of data qudits at phase start
        if pi == 0:
            c.R(dq)
            if basis == "X":
                c.F(dq)
            fresh_set |= set(data); fresh_basis.update({q: ("R" if basis == "Z" else "RX") for q in data})
        if ph.fresh:
            fresh_set |= set(ph.fresh); fresh_basis.update(ph.fresh)
            fz = [dq[qidx[q]] for q, b in ph.fresh.items() if b == "R"]
            fx = [dq[qidx[q]] for q, b in ph.fresh.items() if b == "RX"]
            c.R(fz + fx); noise("XERR", fz + fx, pr)
            if fx:
                c.F(fx); noise("DEP1", fx, pr)
        for r in range(ph.rounds):
            if idle and (r > 0 or pi > 0):
                noise("DEP1", dq, pr)
            c.R(anc); noise("XERR", anc, pr)
            c.F(xanc); noise("DEP1", xanc, pr)
            if flagged:
                fl = list(flag_of.values())
                c.R(fl); noise("XERR", fl, pr)
                if zflags:
                    c.F(zflags); noise("DEP1", zflags, pr)
            for typ in ("Z", "X"):
                for li, L in enumerate(LAYERS[pi][typ]):
                    if li in late[typ]:
                        cp = late[typ][li]; c.CSUM([pq for pq, _ in cp], [w for _, w in cp]); noise("DEP2", [pq for pq, _ in cp], pr)
                    pairs = [((dq[qidx[q]], anc[ci]) if typ == "Z" else (anc[ci], dq[qidx[q]])) for ci, q, w in L]
                    ws = [w for ci, q, w in L]
                    c.CSUM(pairs, ws); noise("DEP2", pairs, pr)
                    if li in early[typ]:
                        cp = early[typ][li]; c.CSUM([pq for pq, _ in cp], [w for _, w in cp]); noise("DEP2", [pq for pq, _ in cp], pr)
            c.FD(xanc); noise("DEP1", xanc, pr)
            noise("XERR", anc, pr)
            m = c.M(anc)
            if flagged:
                if zflags:
                    c.FD(zflags); noise("DEP1", zflags, pr)
                noise("XERR", fl, pr)
                mf = c.M(fl)
                for j, i in enumerate(flagged):
                    ch = checks[i]
                    det([(mf[j], 1)], "X" if ch.typ == "X" else "Z", True, f"p{pi}r{r} flag {ch.typ}{ch.anc_coord}")
            for i, ch in enumerate(checks):
                sector = "X" if ch.typ == "Z" else "Z"
                post = bool(postselect and postselect(ch, pi, r))
                tag = f"p{pi}r{r} {ch.typ}{ch.anc_coord}"
                if ch in prev:
                    det([(m[i], 1), (prev[ch], -1)], sector, post, tag)
                elif r == 0:
                    # determined by preparation?  all legs prepared in the matching basis
                    legs = ch.support
                    if legs <= fresh_set and all(fresh_basis[q] == ("R" if ch.typ == "Z" else "RX") for q in legs):
                        det([(m[i], 1)], sector, post, tag)
                    else:
                        # parent / child relations
                        for parent, children in ph.relations:
                            if any(ci == i for ci, _ in children) and parent in prev:
                                if i == min(ci for ci, _ in children):
                                    det([(m[ci], s) for ci, s in children] + [(prev[parent], -1)], sector, post, tag)
                                break
                        else:
                            # extension of a prepared-state check: colour part equals a previous check, fresh part prepared
                            for pc, mi in list(prev.items()):
                                if pc.typ == ch.typ and pc.support <= legs and (legs - pc.support) <= fresh_set and \
                                        all(fresh_basis[q] == ("R" if ch.typ == "Z" else "RX") for q in legs - pc.support):
                                    sig = unit_ratio(ch, pc, d)
                                    if sig is not None:
                                        det([(m[i], 1), (mi, -sig)], sector, post, tag)
                                        break
            for i, ch in enumerate(checks):
                prev[ch] = m[i]
            fresh_set = set()
    # final readout
    if basis == "X":
        c.FD(dq); noise("DEP1", dq, p)
    noise("XERR", dq, p)
    md = c.M(dq)
    last = phases[-1].code
    for ch in last:
        if ch.typ != basis:
            continue
        sign = 1 if basis == "Z" else -1
        det([(md[qidx[q]], sign * w) for q, w in ch.legs] + [(prev[ch], -1)], "X" if ch.typ == "Z" else "Z", False, f"final {ch.typ}{ch.anc_coord}")
    obs = [(md[qidx[q]], w) for q, w in logical.items()]
    mem = Memory(c, d, 0, sum(ph.rounds for ph in phases), basis, np.array(det_sector), obs, noise_locs)
    mem.det_post = np.array(det_post, dtype=bool)
    mem.det_tags = det_tags
    return mem


LAYERS = {}


def build_dem_masked(mem, sector, mask):
    """surface_zd.build_dem with the faults that fire a postselected detector removed."""
    from .surface_zd import build_dem
    return build_dem(mem, sector, exclude=np.asarray(mask, dtype=int))


class Escape:
    """Colour code (d1 = 5) -> grafted code (r1 rounds) -> final matchable code (r2 rounds) -> readout, or the final
    code's own memory (mode='final'), with optional postselection of the colour-region plaquettes."""

    def __init__(self, d: int = 4, basis: str = "Z", r1: int = 5, r2: int = 5, mode: str = "escape", idle: bool = True,
                 post_colours=(None, None), pre_rounds: int = 1, post_phases=(1,), hyper_post: bool = True,
                 flag_weight: int = 0, flag_phases=(1,)):
        G = load_graft(d)
        self.hyper_post, self.flag_weight, self.flag_phases = hyper_post, flag_weight, flag_phases
        self.G, self.d, self.basis, self.mode = G, d, basis, mode
        codes, coords = G["codes"], G["coords"]
        data_all = sorted(set(q for ch in codes["grafted"] for q in ch.support))
        colour_q = sorted(set(q for ch in codes["colour"] for q in ch.support))
        self.sector = "X" if basis == "Z" else "Z"
        logical = G["ZL"] if basis == "Z" else G["XL"]
        # colour class of a check (by its colour-region support), for postselection
        colours = G["colours"]
        def colour_of(ch):
            s = ch.support & set(colour_q)
            for sup, cl in colours.items():
                if s == sup:
                    return cl
            return None
        cX, cZ = post_colours
        cX = set() if cX is None else ({cX} if isinstance(cX, int) else set(cX))
        cZ = set() if cZ is None else ({cZ} if isinstance(cZ, int) else set(cZ))
        self.post_phases = set(post_phases)
        def postselect(ch, pi, r):
            if mode != "escape" or pi not in self.post_phases:
                return False
            cl = colour_of(ch)
            return cl is not None and ((ch.typ == "X" and cl in cX) or (ch.typ == "Z" and cl in cZ))
        if mode == "escape":
            # transition relations: grafted parent = product of final children (sigma_j)
            rel = []
            for parent in codes["grafted"]:
                if parent in codes["final"]:
                    continue
                parts = [(j, ch) for j, ch in enumerate(codes["final"]) if ch.typ == parent.typ and ch.support <= parent.support]
                for k in (2, 3):
                    for combo in it.combinations(parts, k):
                        if sum(len(ch.support) for _, ch in combo) == len(parent.support) and frozenset().union(*(ch.support for _, ch in combo)) == parent.support:
                            rel.append((parent, [(j, unit_ratio(ch, parent, d)) for j, ch in combo]))
                            break
                    else:
                        continue
                    break
            phases = [Phase(codes["colour"], pre_rounds, noisy=False),
                      Phase(codes["grafted"], r1, fresh=G["fresh"]),
                      Phase(codes["final"], r2, relations=rel)]
            data = colour_q + sorted(G["fresh"])
            # fresh qudits not yet prepared in phase 0: the colour phase must only touch colour qudits (true)
        else:
            phases = [Phase(codes["final"], 0, noisy=False), Phase(codes["final"], r1 + r2)]
            data = data_all
        for pi, ph in enumerate(phases):
            ph.coords = coords
            lay = schedule(ph.code, coords)
            LAYERS[pi] = {"Z": [[e for e in L if ph.code[e[0]].typ == "Z"] for L in lay],
                          "X": [[e for e in L if ph.code[e[0]].typ == "X"] for L in lay]}
            LAYERS[pi] = {t: [L for L in LAYERS[pi][t] if L] for t in "ZX"}
        self.mem = build(d, data, phases, basis, logical, P_BUILD, idle, postselect, flag_weight=flag_weight,
                         flag_phases=flag_phases)
        bad = check_determinism(self.mem.circ, 64, 0)
        self.nondet = [self.mem.det_tags[i] for i in bad]
        if bad:
            raise RuntimeError(f"non-deterministic detectors: {self.nondet[:10]}")
        self.sel = np.nonzero(self.mem.det_sector == self.sector)[0]
        mask = np.nonzero(self.mem.det_post)[0]
        self.mask = mask
        self.dem = build_dem_masked(self.mem, self.sector, mask)
        self.post_idx = np.nonzero(self.mem.det_post[self.sel])[0]        # within the sector
        self.locations = sum(len(qs) for _, _, qs in self.mem.noise_locs)
        self.n_data = len(data)

    def colour_region_patterns(self, phases=(0, 1, 2)):
        """Kept single-fault patterns (dem.singles keys) that fire a detector of a check touching the colour-code qudits."""
        tags = [self.mem.det_tags[i] for i in self.sel]
        colour_q = set(q for ch in self.G["codes"]["colour"] for q in ch.support)
        codes = [self.G["codes"]["colour"], self.G["codes"]["grafted"], self.G["codes"]["final"]] if self.mode == "escape" \
            else [self.G["codes"]["final"]] * 2
        by_coord = [{ch.anc_coord: ch for ch in code} for code in codes]
        cdet = set()
        for j, t in enumerate(tags):
            if not t.startswith("p"):
                continue
            ph = int(t[1:t.index("r")])
            coord = eval(t.split(" ", 1)[1][1:])
            ch = by_coord[ph].get(coord)
            if ch is not None and ph in phases and (ch.support & colour_q):
                cdet.add(j)
        return [k for k in self.dem.singles if any(j in cdet for j, _ in k)]

    def certify_pairs(self, p: float, keys=None, n_jobs: int = 6, cache: str | None = None, scope: str = "colour",
                      degenerate_cap: float | None = None):
        """Second-order certification: decode every pair of colour-region single-fault patterns (scope 'colour') or every
        pair with at least one colour-region member (scope 'colour-all') and store the mis-decoded combined patterns in
        self.discard_pairs (used by run()).  degenerate_cap: a combined pattern that is itself a single-fault pattern
        costs first-order discard mass m1 when postselected; patterns are taken in decreasing order of gain m2/m1 until
        the summed m1 reaches the cap (per unit p; None = no cap).  cache: json path."""
        import json, os
        if cache and os.path.exists(cache):
            c = json.load(open(cache))
            bad = {tuple(tuple(x) for x in k): m for k, m in c["bad"]}
            self.pair_stats = c["stats"]
        else:
            keys = self.colour_region_patterns() if keys is None else keys
            keys_b = list(self.dem.singles) if scope == "colour-all" else None
            first = SZ.misdecoded_patterns(self.dem, PeelDecoder(self.dem, p))      # already postselected at first order
            bad, undet, checked = SZ.misdecoded_pairs(self.dem, p, keys, n_jobs, keys_b=keys_b, exclude=first)
            self.pair_stats = dict(candidates=len(keys), pairs=checked, bad_patterns=len(bad), bad_mass=sum(bad.values()),
                                   undetectable_pair_mass=undet, p=p, scope=scope)
            if cache:
                json.dump(dict(bad=[[list(map(list, k)), m] for k, m in bad.items()], stats=self.pair_stats), open(cache, "w"))
        # first-order cost of each bad pattern: mass of the single faults whose syndrome contains it (pattern_mask semantics)
        singles = self.dem.singles
        m1 = {}
        for sk, od in singles.items():
            ms = sum(od.values())
            for r_ in range(1, len(sk) + 1):
                for sub in it.combinations(sk, r_):
                    if sub in bad:
                        m1[sub] = m1.get(sub, 0.0) + ms
        self.pair_stats.update(degenerate=len(m1), degenerate_m1=sum(m1.values()),
                               degenerate_m2=sum(bad[k] for k in m1))
        if degenerate_cap is not None:
            order = sorted(m1, key=lambda k: -bad[k] / m1[k])
            keep, acc = set(), 0.0
            for k in order:
                if acc + m1[k] > degenerate_cap:
                    continue
                keep.add(k); acc += m1[k]
            bad = {k: m for k, m in bad.items() if k not in m1 or k in keep}
            self.pair_stats.update(degenerate_kept=len(keep), degenerate_kept_m1=acc, cap=degenerate_cap)
        self.discard_pairs = bad
        return self.pair_stats
        return self.pair_stats

    def make_decoder(self, p: float, decoder: str = "peel", **kw):
        if decoder == "peel":
            return PeelDecoder(self.dem, p)
        from .hyperdecoder import HyperDecoder
        return HyperDecoder(self.dem, p, **kw)

    def run(self, p: float, shots: int = 10_000, seed: int = 0, chunk: int = 5_000, decoder: str = "peel", **dec_kw) -> dict:
        dec = self.make_decoder(p, decoder, **dec_kw)
        gaps = []                                     # (gap, failed) per kept shot when the decoder reports a gap
        rng = np.random.default_rng(seed)
        d = self.d
        fails = kept = leftovers = done = hyper_disc = 0
        # first-order certification: single-fault patterns the decoder gets wrong are postselected away
        bad = SZ.misdecoded_patterns(self.dem, dec) if self.hyper_post else {}
        pairs = getattr(self, "discard_pairs", {})
        self.discard_patterns = list(bad) + [k for k in pairs if k not in bad]; bad_mass = sum(bad.values())
        while done < shots:
            n = min(chunk, shots - done)
            sim = FrameSim(self.mem.circ, n, rng).run(scale=p / P_BUILD)
            Dall = sim.detector_values()
            Dv = Dall[:, self.sel]
            O = np.zeros(n, dtype=np.int64)
            for i, cf in self.mem.obs:
                O += cf * sim.meas[:, i]
            O %= d
            # joint postselection: every postselected detector of either sector (the DEM excludes faults firing any of them)
            keep = ~np.any(Dall[:, self.mask] != 0, axis=1) if len(self.mask) else np.ones(n, dtype=bool)
            if self.discard_patterns:
                hm = SZ.pattern_mask(self.discard_patterns, Dv) & keep
                hyper_disc += int(hm.sum()); keep &= ~hm
            for s in np.nonzero(keep)[0]:
                corr, left = dec.decode(Dv[s])
                f = ((O[s] + corr) % d) != 0
                if left is None or isinstance(left, float):
                    if left is not None:
                        gaps.append((left, bool(f)))
                else:
                    leftovers += left
                fails += f
            kept += int(keep.sum()); done += n
        gap_curve = None
        if gaps:
            g = np.array([x for x, _ in gaps]); fl = np.array([y for _, y in gaps])
            gap_curve = [(float(t), float(1 - (kept * (g >= t).mean()) / shots), float(fl[g >= t].sum() / max((g >= t).sum(), 1)), int(fl[g >= t].sum()))
                         for t in np.percentile(g, [0, 10, 20, 30, 40, 50, 60, 70, 80, 90])]
        return dict(mode=self.mode, d=d, p=p, basis=self.basis, shots=shots, kept=kept, fails=int(fails), decoder=decoder, gap_curve=gap_curve,
                    p_fail=fails / max(kept, 1), discard=1 - kept / shots, leftover=int(leftovers), n_det=int(len(self.sel)),
                    edges=len(self.dem.edges), hyper_share=self.dem.hyper_prob / self.dem.total_prob,
                    undetectable_share=self.dem.undetectable_prob / self.dem.total_prob, conflicts=self.dem.conflicts,
                    discard_model=getattr(self.dem, 'discarded_prob', 0.0) * p, hyper_discard=hyper_disc / shots,
                    hyper_post=self.hyper_post, n_hyper=len(getattr(self.dem, 'hyperedges', [])),
                    n_bad_patterns=len(bad), bad_pattern_mass=bad_mass, n_pair_patterns=len(pairs),
                    pair_mass=sum(pairs.values()),
                    locations=self.locations, n_data=self.n_data, flag_weight=self.flag_weight,
                    n_qudits=self.mem.circ.nq)


if __name__ == "__main__":
    import sys, time
    mode = sys.argv[1]; p = float(sys.argv[2]); shots = int(sys.argv[3]); basis = sys.argv[4] if len(sys.argv) > 4 else "Z"
    cX = int(sys.argv[5]) if len(sys.argv) > 5 else None; cZ = int(sys.argv[6]) if len(sys.argv) > 6 else None
    t = time.time()
    ex = Escape(4, basis, mode=mode, post_colours=(cX, cZ))
    print("built", round(time.time() - t, 1), "s; dets", len(ex.mem.circ.detectors), "hyper", ex.dem.hyper_prob / ex.dem.total_prob)
    print(ex.run(p, shots))
