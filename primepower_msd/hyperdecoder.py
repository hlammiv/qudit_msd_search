"""Hypergraph decoder for the Z_{2^k} detector error models of surface_zd / graft_zd: BP + LSD (ldpc) per 2-adic level.

Level j works on the binary syndrome ((sigma >> j) & 1).  Its columns are the distinct level-j bit patterns of the kept
single-fault patterns (DEM.singles) whose values have minimal 2-adic valuation j -- hyperedges included, nothing is
decomposed -- plus, for j >= 1, one "unit-flip" column per level-(j-1) column: a fault decoded at level j-1 with the
wrong unit (Z vs Z^-1 for d = 4) leaves residue 2v, i.e. the pattern's level-(j-1) bits at level j.  Their priors are set
per shot (0.4 for the columns chosen at the previous level, ~0 otherwise).  After each level the chosen patterns'
exact Z_d values (the key whose values best match the syndrome) are subtracted, so the next level sees an even residue.
With gap=True an extra parity row carries the observable bit of every column; the level is decoded with that bit
forced to 0 and to 1, the lighter solution wins and the weight difference is the level's decoder gap (the shot's gap
is the minimum over levels) -- the complementary-gap postselection of the qubit cultivation literature.
Written for d = 4 (only the lambda = -1 unit flip is modelled); other d decode but without the extra unit alternatives.
"""
from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp


def _val2(x: int, d: int) -> int:
    k = int(round(math.log2(d)))
    if x % d == 0:
        return k
    v = 0
    while (x >> v) & 1 == 0:
        v += 1
    return v


class HyperDecoder:
    def __init__(self, dem, p: float, gap: bool = False, max_iter: int = 30, lsd_order: int = 0, bp_method: str = "minimum_sum",
                 ms_scaling_factor: float = 1.0, flip_prior: float = 0.4, tiny: float = 1e-12):
        from ldpc.bplsd_decoder import BpLsdDecoder
        self.d, self.n, self.gap, self.tiny, self.flip_prior = dem.d, dem.n_det, gap, tiny, flip_prior
        d, n = self.d, self.n
        self.k = int(round(math.log2(d)))
        # group singles by (level, bit pattern)
        cols = [dict() for _ in range(self.k)]          # level -> bits(tuple) -> {'mass', 'keys': {key: (mass, obs)}}
        for key, od in dem.singles.items():
            lev = min(_val2(v, d) for _, v in key)
            bits = tuple(j for j, v in key if (v >> lev) & 1)
            c = cols[lev].setdefault(bits, {"mass": 0.0, "keys": {}})
            m = sum(od.values()); c["mass"] += m
            c["keys"][key] = (m, max(od, key=od.get))
        self.levels = []
        prev_base = 0
        for lev in range(self.k):
            base = list(cols[lev].items())
            nb = len(base)
            rows, cc, prior, obs, keys = [], [], [], [], []
            for ci, (bits, c) in enumerate(base):
                rows.extend(bits); cc.extend([ci] * len(bits))
                prior.append(min(0.5, p * c["mass"]))
                dom = max(c["keys"], key=lambda kk: c["keys"][kk][0])
                obs.append(c["keys"][dom][1]); keys.append(c["keys"])
            # unit-flip columns for the previous level's base columns
            nflip = prev_base if lev >= 1 else 0
            if nflip:
                pb = self.levels[lev - 1]
                for fi in range(nflip):
                    bits = pb["bits"][fi]
                    rows.extend(bits); cc.extend([nb + fi] * len(bits))
                    prior.append(tiny)
                    obs.append((-2 * pb["obs"][fi]) % d)       # observable change of lambda = -1
                    keys.append(None)
            ncol = nb + nflip
            Hrows = rows + ([] if not gap else [n] * 0)
            H = sp.csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cc)), shape=(n, ncol))
            if gap:
                ob = np.array([((o >> lev) & 1) for o in obs], dtype=np.uint8)
                H = sp.vstack([H, sp.csr_matrix(ob.reshape(1, -1))]).tocsr()
            dec = BpLsdDecoder(H, error_channel=list(prior), max_iter=max_iter, bp_method=bp_method,
                               ms_scaling_factor=ms_scaling_factor, lsd_order=lsd_order)
            self.levels.append(dict(H=H, dec=dec, prior=np.array(prior), obs=obs, keys=keys, nb=nb, nflip=nflip,
                                    bits=[b for b, _ in base]))
            prev_base = nb

    def _weight(self, lev, sol):
        pr = self.levels[lev]["cur_prior"]
        idx = np.nonzero(sol)[0]
        return float(np.sum(np.log((1 - pr[idx]) / pr[idx])))

    def decode(self, sigma: np.ndarray):
        """sigma: detector charges mod d (one shot) -> (observable correction mod d, gap or None)."""
        d, n = self.d, self.n
        sigma = sigma.copy() % d
        obs = 0
        gap_min = math.inf
        chosen_prev = []
        for lev in range(self.k):
            L = self.levels[lev]
            syn = ((sigma >> lev) & 1).astype(np.uint8)
            pr = L["prior"].copy()
            if L["nflip"]:
                for ci in chosen_prev:
                    pr[L["nb"] + ci] = self.flip_prior
                L["dec"].update_channel_probs(list(pr))
            L["cur_prior"] = pr
            if not syn.any():
                chosen_prev = []
                continue
            if self.gap:
                sols = []
                for bit in (0, 1):
                    s = np.concatenate([syn, [bit]]).astype(np.uint8)
                    sol = np.array(L["dec"].decode(s), dtype=np.uint8)
                    sols.append((self._weight(lev, sol), sol))
                w0, w1 = sols[0][0], sols[1][0]
                sol = sols[0][1] if w0 <= w1 else sols[1][1]
                gap_min = min(gap_min, abs(w0 - w1))
            else:
                sol = np.array(L["dec"].decode(syn), dtype=np.uint8)
            idx = np.nonzero(sol)[0]
            chosen = []
            for ci in idx:
                ci = int(ci)
                if ci < L["nb"]:
                    # exact key: best value match with sigma on its detectors
                    best, bm = None, -1
                    for key, (m, o) in L["keys"][ci].items():
                        match = sum(1 for j, v in key if sigma[j] % d == v % d)
                        if match > bm or (match == bm and m > L["keys"][ci][best][0]):
                            best, bm = key, match
                    for j, v in best:
                        sigma[j] = (sigma[j] - v) % d
                    obs -= L["keys"][ci][best][1]              # correction = minus the fault's observable
                    chosen.append(ci)
                else:
                    fi = ci - L["nb"]
                    bits = self.levels[lev - 1]["bits"][fi]
                    for j in bits:
                        sigma[j] = (sigma[j] - (1 << lev)) % d     # residue 2v at those detectors (d = 4: value 2)
                    obs -= L["obs"][ci]
            chosen_prev = chosen
        return obs % d, (None if not self.gap else (gap_min if gap_min < math.inf else 0.0))

    # PeelDecoder-compatible interface for surface_zd.misdecoded_patterns / misdecoded_pairs
    def decode_corr(self, sigma):
        return self.decode(sigma)[0]
