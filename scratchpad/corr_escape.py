"""Correlated (two-sector) decoding of the f = 1 escape (prototype, 2026-11-27).

The escape decodes each basis run on one detector sector (Z-basis run: the X-sector detectors of the Z checks).  The other
sector's detectors are in the same circuit, and a depolarizing fault X^a Z^b fires both, so the other sector says where the
faults were.  Here (Fowler-style correlated matching, as for qubit surface codes):
  1. decode the other sector with a plain 2-adic peel decoder and take the edge set it traverses (trace);
  2. for every main-sector edge e sharing faults with a traversed other-sector edge f, replace its probability by
       P'_e = P_e - p J(e, f) + J(e, f) / P_f          (J, P_f = sum_lev P_f at unit location probability; P_e absolute)
     i.e. the share of e that comes with f is conditioned on f having happened;
  3. decode the main sector (and the closed-window syndrome decoder) with the complementary gap on the reweighted graphs
     (edge weights replaced in place in the pymatching graphs, restored after the shot).
Joint table: every Pauli of every noise location is enumerated (as build_dem), its main-sector and other-sector projections are
kept when both are graphlike edges of the respective models; postselected faults (excluded detectors) are left out.

    PYTHONPATH=.:scratchpad python scratchpad/corr_escape.py p shots jobs seed [d] [out]
"""
import itertools as it
import json
import math
import os
import sys
import time

import numpy as np

from primepower_msd import surface_zd as SZ
from primepower_msd.graft_run import VARIANTS
from primepower_msd.graft_zd import Escape, build_dem_masked
from primepower_msd.surface_zd import PeelDecoder, _val2, propagate_basis
from primepower_msd.zdsim import FrameSim
from primepower_msd.surface_zd import P_BUILD


def _edge_key(vrow, d):
    """(key, level) of a graphlike sector pattern, or None (no detector -> 'empty', > 2 detectors or mixed levels -> None)."""
    nz = np.nonzero(vrow)[0]
    if len(nz) == 0:
        return "empty"
    if len(nz) > 2:
        return None
    lev = _val2(int(vrow[nz[0]]), d)
    if len(nz) == 2 and _val2(int(vrow[nz[1]]), d) != lev:
        return None
    key = (int(nz[0]), -1) if len(nz) == 1 else (int(nz[0]), int(nz[1]))
    return key, lev


def joint_table(mem, sel_main, sel_other, exclude, main_edges, other_edges, drop=None, prop=None):
    """J[(main key, lev)] -> {other key: joint probability at unit location probability}."""
    d = mem.d
    faults, per_loc, V, O = prop if prop is not None else propagate_basis(mem)
    sel_main = np.asarray(sel_main); sel_other = np.asarray(sel_other)
    exclude = np.asarray(exclude, dtype=int) if exclude is not None else np.zeros(0, int)
    bridge = None
    if drop is not None and len(drop):
        Vs = V[:, sel_main]
        bridge = (V[:, np.asarray(drop, dtype=int)] != 0).any(axis=1) & (Vs != 0).any(axis=1)
    J = {}
    for kind, ids in per_loc:
        m = len(ids)
        ids = np.array(ids)
        Vm = V[ids][:, sel_main].astype(np.int64)
        if bridge is not None:
            Vm = Vm * (~bridge[ids])[:, None]
        Vo = V[ids][:, sel_other].astype(np.int64)
        if not Vm.any() or not Vo.any():
            continue
        Vx = V[ids][:, exclude].astype(np.int64) if len(exclude) else None
        cm = np.nonzero(Vm.any(axis=0))[0]; co = np.nonzero(Vo.any(axis=0))[0]
        U = np.array(list(it.product(range(d), repeat=m)), dtype=np.int64)[1:]
        Sm = U @ Vm[:, cm] % d; So = U @ Vo[:, co] % d
        ok = np.ones(len(U), bool)
        if Vx is not None:
            cx = np.nonzero(Vx.any(axis=0))[0]
            if len(cx):
                ok = ~((U @ Vx[:, cx] % d) != 0).any(axis=1)
        both = ok & (Sm != 0).any(axis=1) & (So != 0).any(axis=1)
        unit = 1.0 / (d ** m - 1)
        for r in np.nonzero(both)[0]:
            vm = np.zeros(len(sel_main), np.int64); vm[cm] = Sm[r]
            vo = np.zeros(len(sel_other), np.int64); vo[co] = So[r]
            km = _edge_key(vm, d); ko = _edge_key(vo, d)
            if km in (None, "empty") or ko in (None, "empty"):
                continue
            if km[0] not in main_edges or ko[0] not in other_edges:
                continue
            dct = J.setdefault(km, {})
            dct[ko[0]] = dct.get(ko[0], 0.0) + unit
    return J


class CorrGap:
    """Gap peel decoder of `dem_main` reweighted by the other sector's decoded edges."""

    def __init__(self, dem_main, dem_other, J, p, cap=0.45):
        self.base = PeelDecoder(dem_main, p, gap=True)
        self.other = PeelDecoder(dem_other, p)
        self.p, self.cap, self.n = p, cap, dem_main.n_det
        self.labels = [SZ.boundary_potential(dem_main, lev)[1] for lev in range(self.base.k)]
        self.edges = dem_main.edges
        Pf = {k: sum(e["P"]) for k, e in dem_other.edges.items()}     # unit location probability: P(e-share | f) = p J / (p Pf)
        self.inv = {}                      # other key -> [(main key, lev, J/P_f - p J)]
        for (mk, lev), dct in J.items():
            for ok, j in dct.items():
                if Pf.get(ok, 0) > 0:
                    self.inv.setdefault(ok, []).append((mk, lev, j / Pf[ok] - p * j))
        self.n_changed = []

    def _w(self, P):
        P = min(max(P, 1e-300), self.cap)
        return math.log((1 - P) / P)

    def _set(self, lev, key, P):
        MG = self.base.gap_matchers[lev]
        a, b = key
        w = self._w(P)
        if b == -1:
            MG.add_edge(a, self.n + self.labels[lev][key], weight=w, merge_strategy="replace")
        else:
            MG.add_edge(a, b, weight=w, merge_strategy="replace")

    def decode(self, sig_main, sig_other):
        self.other.trace = []
        self.other.decode(sig_other)
        used = set(k for k, _ in self.other.trace)
        self.other.trace = None
        delta = {}
        for ok in used:
            for mk, lev, dj in self.inv.get(ok, ()):
                delta[(mk, lev)] = delta.get((mk, lev), 0.0) + dj
        self.n_changed.append(len(delta))
        for (mk, lev), dj in delta.items():
            self._set(lev, mk, self.p * self.edges[mk]["P"][lev] + dj)
        try:
            out = self.base.decode(sig_main)
            self.last_level_gaps = list(self.base.last_level_gaps)
        finally:
            for (mk, lev) in delta:
                self._set(lev, mk, self.p * self.edges[mk]["P"][lev])
        return out


def setup(d, basis, p, variant="one+gap-tick1-closed"):
    ex = Escape(d, basis, **VARIANTS[variant])
    other = "Z" if ex.sector == "X" else "X"
    sel_o = np.nonzero(ex.mem.det_sector == other)[0]
    dem_o = build_dem_masked(ex.mem, other, ex.mask)
    prop = propagate_basis(ex.mem)
    J = joint_table(ex.mem, ex.sel, sel_o, ex.mask, ex.dem.edges, dem_o.edges, prop=prop)
    fin = np.array([i for i, t in enumerate(ex.mem.det_tags) if t.startswith("final")], dtype=int)
    J2 = joint_table(ex.mem, ex.sel_syn, sel_o, ex.mask, ex.dem_syn.edges, dem_o.edges, drop=fin, prop=prop)
    return ex, sel_o, dem_o, J, J2


_ctx = {}


CACHE = os.environ.get("CORR_CACHE", "runs/surface/corr_cache")


def cached_setup(d, basis, p):
    """setup() once per (d, basis), pickled (the basis-fault propagation needs ~3 GB; workers only load the result)."""
    import pickle
    fn = f"{CACHE}/corr_d{d}_{basis}.pkl"
    if os.path.exists(fn):
        return pickle.load(open(fn, "rb"))
    os.makedirs(CACHE, exist_ok=True)
    obj = setup(d, basis, p)
    pickle.dump(obj, open(fn + ".tmp", "wb")); os.replace(fn + ".tmp", fn)
    return obj


def run_chunk(d, basis, p, n, seed):
    if (d, basis) not in _ctx:
        _ctx[(d, basis)] = cached_setup(d, basis, p)
    if (d, basis, p) not in _ctx:                       # decoders + first-order certified patterns, once per worker
        ex, sel_o, dem_o, J, J2 = _ctx[(d, basis)]
        dec = PeelDecoder(ex.dem, p, gap=True)
        _ctx[(d, basis, p)] = (dec, PeelDecoder(ex.dem_syn, p, gap=True), CorrGap(ex.dem, dem_o, J, p),
                               CorrGap(ex.dem_syn, dem_o, J2, p), SZ.misdecoded_patterns(ex.dem, dec))
    ex, sel_o, dem_o, J, J2 = _ctx[(d, basis)]
    dec, dec2, cdec, cdec2, bad = _ctx[(d, basis, p)]
    cdec.n_changed = []
    rng = np.random.default_rng([seed, 0])
    sim = FrameSim(ex.mem.circ, n, rng).run(scale=p / P_BUILD)
    Dall = sim.detector_values()
    O = np.zeros(n, np.int64)
    for i, cf in ex.mem.obs:
        O += cf * sim.meas[:, i]
    O %= d
    keep = ~np.any(Dall[:, ex.mask] != 0, axis=1) if len(ex.mask) else np.ones(n, bool)
    cert = SZ.pattern_mask(list(bad), Dall[:, ex.sel]) & keep if bad else np.zeros(n, bool)
    rows = []
    for s in np.nonzero(keep)[0]:
        Dv, Ds, Do = Dall[s, ex.sel], Dall[s, ex.sel_syn], Dall[s, sel_o]
        c0, g0 = dec.decode(Dv); _, s0 = dec2.decode(Ds)
        c1, g1 = cdec.decode(Dv, Do); _, s1 = cdec2.decode(Ds, Do)
        rows.append((int(s), int((O[s] + c0) % d), float(s0), int((O[s] + c1) % d), float(s1), bool(cert[s]), float(g0), float(g1)))
    return n, rows, float(np.mean(cdec.n_changed)) if cdec.n_changed else 0.0


def summary(res, d):
    """Joint (b' - a' mod d) failures vs min(Z, X) syndrome-gap cut, base vs correlated, on the shots kept by both runs."""
    if len(res["Z"]) == 0 or len(res["X"]) == 0:
        return {}
    iz = {int(r[0]): r for r in res["Z"]}; ix = {int(r[0]): r for r in res["X"]}
    com = sorted(set(iz) & set(ix))
    A = np.array([iz[s] for s in com]); B = np.array([ix[s] for s in com])
    out = {"kept_both": len(com)}
    for name, f, g in (("base", 1, 2), ("corr", 3, 4)):
        bad = ((B[:, f] - A[:, f]) % d) != 0
        gm = np.minimum(A[:, g], B[:, g])
        out[name] = {c: [int((gm >= c).sum()), int((bad & (gm >= c)).sum())] for c in range(0, 14)}
    return out


def main():
    p = float(sys.argv[1]); shots = int(sys.argv[2]); jobs = int(sys.argv[3]); seed = int(sys.argv[4])
    d = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    out = sys.argv[6] if len(sys.argv) > 6 else f"runs/surface/corr_escape_d{d}_p{p:g}_s{seed}.npz"
    from joblib import Parallel, delayed
    CH = 5000
    t = time.time()
    for basis in "ZX":
        t1 = time.time(); cached_setup(d, basis, p); print("setup", basis, round(time.time() - t1), "s", flush=True)
    tasks = [(b, k, min(CH, shots - i), seed * 100000 + k) for k, i in enumerate(range(0, shots, CH)) for b in "ZX"]
    res, done = {"Z": [], "X": []}, set()
    if os.path.exists(out):                      # resume: chunks already saved
        z = np.load(out, allow_pickle=True)
        res = {"Z": list(map(tuple, z["Z"])), "X": list(map(tuple, z["X"]))}; done = set(map(tuple, z["done"]))
        print("resume", len(done), "/", len(tasks), "chunks", flush=True)
    todo = [tk for tk in tasks if (tk[0], tk[1]) not in done]
    B = 2 * jobs
    for i in range(0, len(todo), B):
        batch = todo[i:i + B]
        R = Parallel(n_jobs=jobs)(delayed(run_chunk)(d, b, p, n, sd) for b, k, n, sd in batch)
        for (b, k, n, sd), (_, rows, _) in zip(batch, R):
            res[b] += [(k * CH + r[0],) + r[1:] for r in rows]
            done.add((b, k))
        np.savez_compressed(out + ".tmp.npz", Z=np.array(res["Z"], float), X=np.array(res["X"], float), done=np.array(sorted(done), dtype=object),
                            p=p, shots=shots, d=d, seed=seed, cols="shot, flip_base, gapsyn_base, flip_corr, gapsyn_corr, cert, gap_base, gap_corr")
        os.replace(out + ".tmp.npz", out)
        sm = summary(res, d)
        el = time.time() - t; nd = i + len(batch)
        line = " ".join(f"c{c}: {sm['base'][c][1]}/{sm['base'][c][0]} vs {sm['corr'][c][1]}/{sm['corr'][c][0]}" for c in (0, 4, 6, 7, 8, 10)) if sm else ""
        print(f"progress {len(done)}/{len(tasks)} chunks, {el/60:.1f} min, eta {el/nd*(len(todo)-nd)/60:.1f} min | fails/kept base vs corr {line}", flush=True)
    sm = summary(res, d)
    json.dump(dict(p=p, d=d, seed=seed, shots=shots, **{k: v for k, v in sm.items()}), open(out[:-4] + ".json", "w"))
    print("saved", out, flush=True)


if __name__ == "__main__":
    main()
