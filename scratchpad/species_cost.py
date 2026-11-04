"""Species accounting for p = 7 protocols (U^mu magic states, mu in the three classes of F_7^*/{+-1}).

Part 1 checks, by dense state-vector arithmetic, the gate-teleportation conversion
|U^a> + |U^b> -> |U^(a+b)> (Clifford correction, Z errors add), and the error composition.
Part 2 checks that rescaling charges maps a code on species t to the same code on species s*t.
Part 3 recomputes protocol costs in two scenarios by a forward Pareto search over (species, eps, cost):
   multi  : raw noisy states of all three species are available at eps0, cost 1 each
   single : only raw U at eps0; other species come from conversions (or from towers fed by conversions)
A protocol round takes every input at error <= eps_in.  Candidates are run twice: with f_hi, pacc_lo
(a proved upper bound on cost) and with the leading-order surrogate c eps^d, (1-eps)^n (an estimate).
Writes runs/weighted/species_cost.json.
"""
import json, math, sys
from collections import Counter

import mpmath
import numpy as np

from qmsd.doubling import paper_codes
from qmsd.objective import standard_library, DEFAULT_GRID
from qmsd.protocol import best_protocol
from qmsd.rigorous_curve import from_census
from qmsd.triorthogonal import logical_lambdas
from qmsd.weighted import quadric, weighted_support

P = 7
W = np.exp(2j * np.pi / P)
CUBE = [pow(x, 3, P) for x in range(P)]
SPECIES = {1: 0, 6: 0, 3: 1, 4: 1, 2: 2, 5: 2}      # discrete log base 3 mod 3; multiplication adds
COST_NAME = {0: "U (class 1,6)", 1: "U^3 (class 3,4)", 2: "U^2 (class 2,5)"}


def magic(a):
    return np.array([W ** (a * CUBE[x] % P) for x in range(P)]) / math.sqrt(P)


def part1():
    """Teleport U^b into a data qudit holding psi, using ancilla |U^b> Z^c, data error Z^e.
    Map |s>|x> -> |s>|s-x>, measure y = s-x: the data keeps phase w^(b (s-y)^3 + c (s-y)).
    Correct with the quadratic phase w^(3 b s^2 y - 3 b s y^2) and the Z power w^(c y) is global."""
    worst = 0.0
    rng = np.random.default_rng(1)
    for a in range(1, P):
        for b in range(1, P):
            if (a + b) % P == 0:
                continue
            for c in range(P):
                for e in (0, 3):
                    psi = magic(a) * np.array([W ** (e * s % P) for s in range(P)])
                    anc = magic(b) * np.array([W ** (c * x % P) for x in range(P)])
                    for y in range(P):
                        out = np.zeros(P, dtype=complex)
                        for s in range(P):
                            out[s] = psi[s] * anc[(s - y) % P]
                        out /= np.linalg.norm(out)
                        corr = np.array([W ** ((3 * b * s * s * y - 3 * b * s * y * y) % P) for s in range(P)])
                        out = out * corr
                        tgt = magic(a + b) * np.array([W ** ((e + c) * s % P) for s in range(P)])
                        worst = max(worst, 1 - abs(np.vdot(tgt, out)) ** 2)
    # error composition: eps1, eps2 with uniform nonzero Z power
    e1, e2 = 1e-2, 3e-2
    d1 = np.full(P, e1 / (P - 1)); d1[0] = 1 - e1
    d2 = np.full(P, e2 / (P - 1)); d2[0] = 1 - e2
    conv = np.array([sum(d1[a] * d2[(z - a) % P] for a in range(P)) for z in range(P)])
    uniform = float(np.ptp(conv[1:]))
    formula = compose(e1, e2)
    return {"teleport_worst_infidelity": float(worst), "composition_nonuniformity": uniform,
            "composition_formula_error": abs(1 - conv[0] - formula)}


def compose(e1, e2):
    """Error of Z^a Z^b when a, b are 0 w.p. 1 - e, uniform nonzero w.p. e (again uniform on nonzero)."""
    return e1 + e2 - e1 * e2 * P / (P - 1)


def part2():
    """Scaling: triple sums of the [[20,1,4]]_7 stabilisers vanish for c * coeffs and lambda scales by c."""
    out = {}
    for label, G0, H1, co in paper_codes(P):
        G0 = np.asarray(G0, dtype=np.int64) % P
        for c in (1, 2, 3):
            cc = (c * np.asarray(co, dtype=np.int64)) % P
            worst = 0
            r = G0.shape[0]
            for i in range(r):
                for j in range(i, r):
                    for l in range(j, r):
                        worst = max(worst, int(np.sum(cc * G0[i] * G0[j] * G0[l]) % P))
            out[f"{label} c={c}"] = {"max_triple_sum": worst, "lambda": [int(x) for x in logical_lambdas(H1, len(H1), P, cc)]}
    return out


class Tagged:
    """A curve with its input species counts relative to the output class (t -> count)."""

    def __init__(self, curve, counts):
        self.c, self.counts = curve, counts
        assert sum(counts.values()) == curve.n, (curve.name, counts, curve.n)
        self.name = curve.name


def tag_library():
    co = {label.split()[0]: (co, H1) for label, _, H1, co in paper_codes(P)}
    tagged = []
    for c in standard_library(P):
        key = c.name.split()[0]
        if c.name.endswith("dbl") and key in co:
            cf, H1 = co[key]
            lam = logical_lambdas(H1, len(H1), P, cf)[0] % P
            rel = Counter(SPECIES[int(x) % P * pow(int(lam), -1, P) % P] for x in cf)
            tagged.append(Tagged(c, dict(rel)))
        else:
            tagged.append(Tagged(c, {0: c.n}))
    return tagged


class Pareto:
    def __init__(self):
        self.pts = []       # sorted by eps descending, cost ascending

    def add(self, eps, cost, how):
        eps = float(eps) if eps > 1e-300 else 1e-300
        for e, c, _ in self.pts:
            if e <= eps and c <= cost:
                return False
        self.pts = [t for t in self.pts if not (eps <= t[0] and cost <= t[1])]
        self.pts.append((eps, cost, how))
        self.pts.sort(key=lambda t: -t[0])
        return True

    def lookup(self, eps):
        best = None
        for e, c, h in self.pts:
            if e <= eps * (1 + 1e-12) and (best is None or c < best[0]):
                best = (c, e, h)
        return best

    def thin(self, target, cap=240):
        """Keep the target-reaching end exact; thin the rest by log eps (valid: only removes options)."""
        if len(self.pts) <= cap:
            return
        logs = np.array([math.log10(t[0]) for t in self.pts])
        keep, last = [], None
        step = (logs.max() - logs.min()) / cap
        for t, l in zip(self.pts, logs):
            if last is None or last - l >= step:
                keep.append(t); last = l
        keep.append(self.pts[-1])
        self.pts = sorted(set(keep), key=lambda t: -t[0])


def search(tagged, eps0, target, scenario, sweeps=6, mode="upper"):
    """Forward search: returns (cost, how) for species 0 at error <= target.  Every recorded point
    is a real protocol: inputs supplied at <= eps_in, output error bound f_hi(eps_in)."""
    F = [Pareto() for _ in range(3)]
    F[0].add(eps0, 1.0, "raw")
    if scenario == "multi":
        F[1].add(eps0, 1.0, "raw")
        F[2].add(eps0, 1.0, "raw")
    cache = {}

    def fp(t, e):
        k = (t.name, e)
        if k not in cache:
            c = t.c
            if mode == "upper" or not hasattr(c, "f_surrogate"):    # library curves are exact in both modes
                f, pa = (c.f_hi, c.pacc_lo)
                cache[k] = (float(f(e)), float(pa(e)))
            else:                                                   # leading-order estimate, not a bound
                cache[k] = (float(c.f_surrogate(e)), float((1 - e) ** c.n))
        return cache[k]

    pairs = [(1, 1, 2), (1, 2, 3), (2, 2, 4), (1, 3, 4), (3, 3, 6), (2, 3, 5), (1, 4, 5)]
    for _ in range(sweeps):
        changed = False
        grid = sorted({t[0] for f in F for t in f.pts}, reverse=True)
        grid = [e for e in grid if e > target * 1e-3]
        if scenario == "single":            # conversions, any pair of frontier points
            for a, b, ab in pairs:
                sa, sb, s = SPECIES[a], SPECIES[b], SPECIES[ab]
                for ea, ca, _ in list(F[sa].pts):
                    for eb, cb, _ in list(F[sb].pts):
                        changed |= F[s].add(compose(ea, eb), ca + cb, f"conv({a}+{b})@{ea:.3g},{eb:.3g}")
        for e in grid:
            for t in tagged:
                fe, pa = fp(t, e)
                if pa <= 0 or fe >= e:
                    continue
                for out in range(3):        # scaled copy with output species ``out``
                    cost, ok = 0.0, True
                    for rel, cnt in t.counts.items():
                        L = F[(out + rel) % 3].lookup(e)
                        if L is None:
                            ok = False
                            break
                        cost += cnt * L[0]
                    if ok:
                        changed |= F[out].add(fe, cost / (t.c.k * pa), f"{t.name}@{e:.3g}")
        for f in F:
            f.thin(target)
        if not changed:
            break
    return F


def answer(F, target):
    L = F[0].lookup(target)
    return None if L is None else {"cost": L[0], "eps": L[1], "last": L[2]}


def main():
    r = {"part1": part1()}
    print(r["part1"], flush=True)
    r["part2"] = part2()
    print(r["part2"], flush=True)

    tagged = tag_library()
    r["library_species"] = {t.name: t.counts for t in tagged}
    print(r["library_species"], flush=True)

    GX, X, mu = weighted_support(P, 4, quadric(P, 4, (1, 1, 1, 1)), 7)
    cnt = Counter(SPECIES[int(x)] for x in mu)
    r["sphere_class_sizes"] = dict(cnt)
    print("sphere classes", cnt, flush=True)
    assert cnt[0] == cnt[1] == cnt[2] == 672

    cands = {}
    for path, kind in [("runs/weighted/p7_m4_sphere_1746_270_7_census.json", "sphere"),
                       ("runs/weighted/p7_m4_sphere_1746_270_6_census.json", "sphere"),
                       ("runs/weighted/p7_m4_sphere_1756_260_7_census.json", "sphere"),
                       ("runs/weighted/p7_m4_unit_2093_308_5_census.json", "unit")]:
        c = from_census(path)
        counts = {0: c.n} if kind == "unit" else {0: 672 - c.k, 1: 672, 2: 672}
        cands[c.name] = Tagged(c, counts)
        weights = {0: 1, 1: 3, 2: 2}
        n_eff = sum(weights[s] * v for s, v in counts.items())
        print(c.name, counts, "U per round", n_eff, "gamma", round(c.gamma, 3), "gamma_U-count", round(math.log(n_eff / c.k) / math.log(c.d), 3), flush=True)
        r.setdefault("gamma", {})[c.name] = {"gamma": c.gamma, "n_U": n_eff, "gamma_U_count": math.log(n_eff / c.k) / math.log(c.d)}

    rows = []
    base_raw = standard_library(P)
    all_cands = list(cands.values())
    for eps0, target in DEFAULT_GRID:
        ref = best_protocol(base_raw, eps0, target, max_depth=5)
        row = {"eps0": eps0, "target": target, "best_protocol_base": ref.cost}
        for sc in ("multi", "single"):
            b = answer(search(tagged, eps0, target, sc), target)
            row[sc] = {"base": b}
            print(f"{eps0:g}->{target:g} {sc} base {b['cost']:.2f} last {b['last']} (best_protocol {ref.cost:.2f})", flush=True)
            for name, t in list(cands.items()) + [("all", None)]:
                pool = tagged + (all_cands if t is None else [t])
                up = answer(search(pool, eps0, target, sc, mode="upper"), target)
                lo = answer(search(pool, eps0, target, sc, mode="surrogate"), target)
                row[sc][name] = {"proved": up, "estimate": lo, "saving_proved": b["cost"] / up["cost"], "saving_estimate": b["cost"] / lo["cost"]}
                print(f"   {name}: proved {up['cost']:.2f} ({b['cost'] / up['cost']:.2f}x) [{up['last']}]  estimate {lo['cost']:.2f} ({b['cost'] / lo['cost']:.2f}x) [{lo['last']}]", flush=True)
        rows.append(row)
    r["rows"] = rows
    json.dump(r, open("runs/weighted/species_cost.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
