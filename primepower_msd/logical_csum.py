"""Logical CSUM between two signed Z_d rotated surface-code patches by a transversal CSUM layer (surface_zd model).

Two identical D-patches (control, target) are prepared in the memory basis, run `rounds` syndrome rounds, a transversal
CSUM(control_q -> target_q, w) layer with DEP2 noise is applied, then `rounds` more rounds and the data are read out.
Identical layouts and signs make the transversal CSUM a logical CSUM: X-check_c -> X-check_c X-check_t^w and
Z-check_t -> Z-check_t Z-check_c^{-w} are products of checks, so the detectors across the layer are
    m_t^Z(r) - m_t^Z(r-1) - w m_c^Z(r-1)   and   m_c^X(r) - m_c^X(r-1) + w m_t^X(r-1)
(the Z ancilla records +sum w x, the X ancilla -sum w z).  Logical action: X_c -> X_c X_t^w, Z_t -> Z_t Z_c^{-w}, so for
|0,0> and |+,+> both logical readouts are deterministic and the experiment fails when either is wrong.  A fault on the
control before the layer fires detectors in both patches (a weight-4 hyperedge) and is decomposed stim-style by
build_dem.  The observables are decoded with one DEM per observable (same graph, same matching, different observable
coefficients).  csum=False gives the baseline: two independent 2*rounds memories.

    python -m primepower_msd.logical_csum sweep d jobs [Dlist] [plist] [shots_low shots_high]   -> runs/surface/csum_d{d}.jsonl
    python -m primepower_msd.logical_csum fit d
"""
from __future__ import annotations

import dataclasses
import json
import math
import sys
import time
from collections import defaultdict

import numpy as np

from .surface_zd import P_BUILD, Memory, PeelDecoder, build_dem, propagate_basis, rotated_code
from .zdsim import Circuit, FrameSim, check_determinism


def csum_circuit(D: int, d: int, rounds_before: int, rounds_after: int, p: float, basis: str = "Z", idle: bool = True,
                 csum: bool = True, w: int = 1) -> Memory:
    checks = rotated_code(D, d)
    c = Circuit(d)
    data = [c.q(D * D), c.q(D * D)]
    anc = [c.q(len(checks)), c.q(len(checks))]
    xanc = [[a[i] for i, ch in enumerate(checks) if ch.typ == "X"] for a in anc]
    sector = {"X": "Z", "Z": "X"}
    det_sector, locs = [], []

    def noise(kind, qs, pr):
        if pr > 0:
            getattr(c, kind)(qs, pr)
            locs.append((len(c.ops) - 1, kind, list(qs)))

    alld, alla, allx = data[0] + data[1], anc[0] + anc[1], xanc[0] + xanc[1]
    c.R(alld); noise("XERR", alld, p)
    if basis == "X":
        c.F(alld); noise("DEP1", alld, p)
    prev, csum_op = None, None
    for r in range(rounds_before + rounds_after):
        if r == rounds_before:
            csum_op = len(c.ops)
            if csum:
                pairs = list(zip(data[0], data[1]))
                c.CSUM(pairs, w); noise("DEP2", pairs, p)
        if idle and r > 0:
            noise("DEP1", alld, p)
        c.R(alla); noise("XERR", alla, p)
        c.F(allx); noise("DEP1", allx, p)
        for step in range(4):
            pairs, ws = [], []
            for P in range(2):
                for i, ch in enumerate(checks):
                    for q, wq, s in ch.legs:
                        if s != step:
                            continue
                        pairs.append((data[P][q], anc[P][i]) if ch.typ == "Z" else (anc[P][i], data[P][q])); ws.append(wq)
            c.CSUM(pairs, ws); noise("DEP2", pairs, p)
        c.FD(allx); noise("DEP1", allx, p)
        noise("XERR", alla, p)
        m = [c.M(anc[0]), c.M(anc[1])]
        for P in range(2):
            for i, ch in enumerate(checks):
                if prev is None:
                    if ch.typ == basis:
                        c.detector([(m[P][i], 1)], f"P{P} {ch.typ}{i} r0"); det_sector.append(sector[ch.typ])
                    continue
                terms = [(m[P][i], 1), (prev[P][i], -1)]
                if r == rounds_before and csum:
                    if P == 1 and ch.typ == "Z":
                        terms.append((prev[0][i], -w))
                    if P == 0 and ch.typ == "X":
                        terms.append((prev[1][i], w))
                c.detector(terms, f"P{P} {ch.typ}{i} r{r}"); det_sector.append(sector[ch.typ])
        prev = m
    if basis == "X":
        c.FD(alld); noise("DEP1", alld, p)
    noise("XERR", alld, p)
    md = [c.M(data[0]), c.M(data[1])]
    sign = 1 if basis == "Z" else -1
    for P in range(2):
        for i, ch in enumerate(checks):
            if ch.typ == basis:
                c.detector([(md[P][q], sign * wq) for q, wq, _ in ch.legs] + [(prev[P][i], -1)], f"P{P} {ch.typ}{i} final")
                det_sector.append(sector[ch.typ])
    if basis == "Z":
        obs = [[(md[P][i * D], 1) for i in range(D)] for P in range(2)]      # logical Z = column 0
    else:
        obs = [[(md[P][j], 1) for j in range(D)] for P in range(2)]          # logical X = row 0
    mem = Memory(c, d, D, rounds_before + rounds_after, basis, np.array(det_sector), obs[0], locs)
    mem.obs_all, mem.csum_op, mem.data, mem.anc = obs, csum_op, data, anc
    return mem


class SequentialDecoder:
    """Pauli-frame (correlated) decoding of the transversal CSUM.  Errors of one sector flow one way through the layer
    (X: control -> target, Z: target -> control), so the source patch is decoded from its own detectors alone, the
    inferred faults are propagated through the layer onto the sink patch's detectors and observable (the shadow of
    every source edge, read off the basis-fault propagation; the most probable shadow when a source edge is shared by
    faults before and after the layer), and the sink is decoded with a DEM of sink-only faults."""

    def __init__(self, mem: Memory, sector: str, p: float):
        d = mem.d
        self.d = d
        src, snk = (0, 1) if sector == "X" else (1, 0)
        tags = mem.circ.det_tags
        patch = np.array([int(t[1]) for t in tags])
        sel = np.nonzero(mem.det_sector == sector)[0]
        self.pos_src = np.nonzero(patch[sel] == src)[0]
        self.pos_snk = np.nonzero(patch[sel] == snk)[0]
        ds = mem.det_sector.copy(); ds[sel[patch[sel] == snk]] = "N"
        dem_src = build_dem(dataclasses.replace(mem, det_sector=ds, obs=mem.obs_all[src]), sector)
        snk_q = set(mem.data[snk] + mem.anc[snk])
        locs = []
        for k, kind, qs in mem.noise_locs:
            keep = [q for q in qs if (set(q) if kind == "DEP2" else {q}) <= snk_q]
            if keep:
                locs.append((k, kind, keep))
        ds2 = mem.det_sector.copy(); ds2[sel[patch[sel] == src]] = "N"
        dem_snk = build_dem(dataclasses.replace(mem, det_sector=ds2, obs=mem.obs_all[snk], noise_locs=locs), sector)
        self.dec_src, self.dec_snk = PeelDecoder(dem_src, p), PeelDecoder(dem_snk, p)
        self.dem_src, self.dem_snk = dem_src, dem_snk
        # shadow of every source edge: sink detector pattern and sink-observable coefficient per unit multiplier of the edge's pattern
        faults, per_loc, V, O = propagate_basis(dataclasses.replace(mem, obs=mem.obs_all[snk]))
        Vs = V[:, sel]
        units = [lam for lam in range(1, d) if lam % 2]
        inv = {lam: next(mu for mu in units if (lam * mu) % d == 1) for lam in units}
        votes = {}
        for f in range(len(faults)):
            vs = Vs[f, self.pos_src]; vk = Vs[f, self.pos_snk]
            dets = np.nonzero(vs)[0]
            if len(dets) == 0 or len(dets) > 2:
                continue
            key = (int(dets[0]), int(dets[1]) if len(dets) == 2 else -1)
            e = dem_src.edges.get(key)
            if e is None:
                continue
            ca, cb, _ = e["pat"]
            lam = next((l for l in units if (l * ca - vs[key[0]]) % d == 0 and (key[1] == -1 or (l * cb - vs[key[1]]) % d == 0)), None)
            if lam is None:
                continue
            mu = inv[lam]
            shadow = tuple((int(j), int(mu * vk[j] % d)) for j in np.nonzero(vk)[0])
            o = int(mu * O[f] % d)
            vt = votes.setdefault(key, {})
            vt[(shadow, o)] = vt.get((shadow, o), 0) + 1
        self.shadow = {key: max(vt.items(), key=lambda kv: kv[1])[0] for key, vt in votes.items()}
        self.ambiguous = sum(1 for vt in votes.values() if len(vt) > 1)
        self.n_shadow = sum(1 for sh, o in self.shadow.values() if sh or o)

    def decode(self, sigma):
        """sigma over the sector detectors (both patches).  Returns (correction for the source observable, for the sink)."""
        d = self.d
        s_src = sigma[self.pos_src]; s_snk = sigma[self.pos_snk].copy()
        self.dec_src.trace = []
        o_src, _ = self.dec_src.decode(s_src)
        shift = 0
        for key, w in self.dec_src.trace:
            sh = self.shadow.get(key)
            if sh is None:
                continue
            vec, o = sh
            for j, v in vec:                      # the correction's shadow is minus the actual faults' shadow
                s_snk[j] = (s_snk[j] + w * v) % d
            shift += w * o
        o_snk, _ = self.dec_snk.decode(s_snk)
        return o_src % d, (o_snk + shift) % d


class CsumExperiment:
    """Transversal logical CSUM (or its no-gate baseline) between two D-patches; run() fails when either logical is wrong."""

    def __init__(self, D: int, d: int, basis: str = "Z", rounds: int | None = None, csum: bool = True, idle: bool = True,
                 decoder: str = "indep"):
        self.D, self.d, self.basis, self.csum, self.decoder = D, d, basis, csum, decoder
        self.rounds = D if rounds is None else rounds
        self.mem = csum_circuit(D, d, self.rounds, self.rounds, P_BUILD, basis, idle, csum)
        bad = check_determinism(self.mem.circ, 64, 0)
        assert bad == [], f"non-deterministic detectors {[self.mem.circ.det_tags[i] for i in bad[:8]]}"
        self.sector = "X" if basis == "Z" else "Z"
        self.sel = np.nonzero(self.mem.det_sector == self.sector)[0]
        self.dems = [build_dem(dataclasses.replace(self.mem, obs=o), self.sector) for o in self.mem.obs_all]
        self.locations = sum(len(qs) for _, _, qs in self.mem.noise_locs)

    def make_decoder(self, p: float):
        """Returns a function sigma -> (correction for observable 0 (control), observable 1 (target))."""
        if self.decoder == "seq":
            sq = SequentialDecoder(self.mem, self.sector, p)
            self.seq = sq
            if self.sector == "X":
                return lambda s: sq.decode(s)
            return lambda s: sq.decode(s)[::-1]
        decs = [PeelDecoder(dem, p) for dem in self.dems]
        return lambda s: (decs[0].decode(s)[0], decs[1].decode(s)[0])

    def run(self, p: float, shots: int = 10_000, seed: int = 0, chunk: int = 10_000) -> dict:
        dec = self.make_decoder(p)
        rng = np.random.default_rng(seed)
        d = self.d
        fails = [0, 0]; both = done = 0
        while done < shots:
            n = min(chunk, shots - done)
            sim = FrameSim(self.mem.circ, n, rng).run(scale=p / P_BUILD)
            Dv = sim.detector_values()[:, self.sel]
            O = []
            for o in self.mem.obs_all:
                v = np.zeros(n, dtype=np.int64)
                for i, cf in o:
                    v += cf * sim.meas[:, i]
                O.append(v % d)
            for s in range(n):
                corr = dec(Dv[s])
                bad = [((O[k][s] + corr[k]) % d) != 0 for k in range(2)]
                fails[0] += bad[0]; fails[1] += bad[1]; both += any(bad)
            done += n
        dem = self.dems[0]
        return dict(D=self.D, d=d, p=p, rounds=self.rounds, basis=self.basis, csum=self.csum, decoder=self.decoder, shots=shots, fails=int(both),
                    fails_control=int(fails[0]), fails_target=int(fails[1]), p_fail=both / shots, n_det=int(len(self.sel)),
                    edges=len(dem.edges), hyper_share=dem.hyper_prob / dem.total_prob, conflicts=dem.conflicts,
                    undetectable_share=dem.undetectable_prob / dem.total_prob, locations=self.locations)


def point(D, d, basis, p, shots, csum, seed, decoder="indep"):
    t = time.time()
    try:
        r = CsumExperiment(D, d, basis, csum=csum, decoder=decoder).run(p, shots, seed)
    except MemoryError as e:                       # keep the sweep alive under ulimit -v
        r = dict(D=D, d=d, p=p, basis=basis, csum=csum, decoder=decoder, shots=0, fails=0, error=repr(e))
    r["seconds"] = round(time.time() - t, 1)
    return r


def load(d):
    agg = defaultdict(lambda: [0, 0])
    for l in open(f"runs/surface/csum_d{d}.jsonl"):
        r = json.loads(l)
        key = ("seq" if r.get("decoder", "indep") == "seq" else ("csum" if r["csum"] else "base"))
        a = agg[(key, r["basis"], r["D"], r["p"])]; a[0] += r["fails"]; a[1] += r["shots"]
    return agg


def fit(pts):
    """log p_fail = log A + ((D+1)/2) (log p - log p_th) over sub-threshold points (p_fail < 0.1, >= 5 events)."""
    X, y, w = [], [], []
    for (D, p), (f, s) in pts.items():
        if f < 5 or f / s >= 0.1:
            continue
        X.append([1.0, (D + 1) / 2]); y.append(math.log(f / s) - (D + 1) / 2 * math.log(p)); w.append(math.sqrt(f))
    if len(X) < 4:
        return None
    X, y, w = np.array(X), np.array(y), np.array(w)
    sol, *_ = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)
    resid = y - X @ sol
    return dict(A=math.exp(sol[0]), p_th=math.exp(-sol[1]), n=len(y), rms=float(np.sqrt(np.mean(resid ** 2))))


def main():
    cmd = sys.argv[1]; d = int(sys.argv[2])
    if cmd == "sweep":
        from joblib import Parallel, delayed
        jobs = int(sys.argv[3])
        Ds = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [3, 5, 7, 9]
        ps = [float(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else [1e-3, 2e-3, 3e-3, 4e-3, 6e-3]
        lo, hi = (int(sys.argv[6]), int(sys.argv[7])) if len(sys.argv) > 7 else (100_000, 40_000)
        decoder = sys.argv[8] if len(sys.argv) > 8 else "indep"
        css = (True, False) if decoder == "indep" else (True,)
        tasks = [(D, basis, p, lo if p <= 2e-3 else hi, cs) for D in Ds for basis in "ZX" for p in ps for cs in css]
        tasks.sort(key=lambda t: -t[0])
        out = f"runs/surface/csum_d{d}.jsonl"
        gen = Parallel(n_jobs=jobs, verbose=0, return_as="generator_unordered")(
            delayed(point)(D, d, b, p, s, cs, hash((D, b, p, cs, decoder)) % 2**31, decoder) for D, b, p, s, cs in tasks)
        for r in gen:                                  # one line per finished point, so a crash loses nothing done
            if "error" in r:
                print("ERROR", r["D"], r["basis"], r["p"], r["error"], flush=True)
                continue
            with open(out, "a") as f:
                f.write(json.dumps(r) + "\n")
            print(r["D"], r["basis"], r["p"], r["decoder"], "csum" if r["csum"] else "base", r["fails"], r["fails_control"], r["fails_target"],
                  r["shots"], "hyper %.4f" % r["hyper_share"], "s", r["seconds"], flush=True)
        return
    agg = load(d)
    for cs in ("csum", "seq", "base"):
        for basis in "ZX":
            pts = {(D, p): tuple(v) for (c_, b, D, p), v in agg.items() if c_ == cs and b == basis}
            if not pts:
                continue
            print(cs, basis, fit(pts))
            for D in sorted({D for D, _ in pts}):
                row = ["D=%d" % D]
                for p in sorted({p for _, p in pts}):
                    if (D, p) in pts:
                        f, s = pts[(D, p)]; row.append("%.0e: %.2e (%d)" % (p, f / s, f))
                print("   ", "  ".join(row))


if __name__ == "__main__":
    main()
