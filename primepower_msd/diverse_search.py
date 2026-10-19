"""Driver for the diverse-gate ring-MSD search: evaluate candidate Z_4 codes for hosting of
(A) single-qudit level-3 layers (universality via sector escape of the INDUCED logical) and
(B) entangling layers (in-block ring-CS/CCZ, across-block CCZ / 2-body tables), with distance, gamma,
screens, fault-tolerance (protocol) flag and provenance; write a JSON catalog.

Experiments (see plan): E2 lifts, E3 Galois-ring GR(4,m) graph/evaluation codes, E4 constrained sampler
(X-generators drawn from the q-weight-0 set, S2/S3 enforced), E5 across-block ring-CCZ codes.

    python -m primepower_msd.diverse_search --exp E2 --out runs/diverse/catalog_E2.json
    python -m primepower_msd.diverse_search --exp E4 --n 12 --k 2 --trials 20000 --jobs 8

Resource hygiene: distance MITM capped (`--dmax`), left-table guard, `--mem-gb` RLIMIT_AS cap.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import time

import numpy as np

from . import phase_tables as PT
from . import rgt
from .symplectic import Z4Code, blocks
from .z4_distance import code_distance


# --------------------------------------------------------------------------- evaluation
def _layer_summary(code: Z4Code, layer: rgt.TransversalLayer, mode: str, block_of=None) -> dict:
    v = rgt.hosts(code, layer, mode)
    out = {"label": layer.label, "hosted": bool(v.hosted), "exact": bool(v.exact), "mode": mode,
           "physical_escapes": layer.escapes_physically()}
    out.update({f"ft_{k}": val for k, val in layer.degree_report(block_of).items()})
    if v.hosted and code.k_genuine <= 6 and code.k_genuine >= 1:
        rep = rgt.induced_logical_report(code, layer)
        out["logical"] = {k: rep[k] for k in ("level", "escapes_A2", "escapes_all_sectors", "entangling",
                                               "odd_third_difference", "genuine_T", "mixed_parts_mod_d")}
        out["logical"]["cs_terms"] = len(rep["cs_terms"]); out["logical"]["ccz_terms"] = len(rep["ccz_terms"])
        out["logical"]["table"] = rep["table"] if len(rep["table"]) <= 256 else None
    return out


def evaluate(code: Z4Code, d_max: int = 4, single_tables=None, entangling: bool = True,
             across_block: bool = True, provenance: str = "") -> dict:
    """Full report for one code."""
    t0 = time.time()
    rep = {"label": code.label, "provenance": provenance, "n": code.n, "d_ring": code.d,
           "k": code.k, "k_genuine": code.k_genuine, "soft": code.soft, "free_quotient_ok": code.free_quotient_ok,
           "is_css": code.is_css, "isotropic": code.isotropic}
    try:
        dmin, cert, det = code_distance(code, d_max=d_max)
    except MemoryError as e:
        dmin, cert, det = None, False, {"error": str(e)}
    rep.update({"distance": dmin, "distance_certified": cert, "distance_detail": det})
    rep["gamma"] = (math.log(code.n / code.k) / math.log(dmin)) if (dmin and dmin > 1 and code.k > 0) else None
    rep["axis_parity"] = rgt.axis_parity_report(code)
    # (A) single-qudit layers
    tabs = single_tables if single_tables is not None else PT.strict_level3_single(code.d)
    A = []
    s1_pass = 0
    for q in tabs:
        scr = rgt.uniform_screens(code, q)
        s1_pass += scr["S1_qweight"]
        if not scr["all"]:
            continue
        A.append(_layer_summary(code, rgt.TransversalLayer.uniform(code.n, np.array(q), code.d, f"uniform{list(q)}"),
                                "support" if code.is_css else "cosets"))
    rep["single_qudit"] = {"n_tables": len(tabs), "S1_pass": s1_pass, "screens_pass": len(A),
                           "hosted": sum(a["hosted"] for a in A),
                           "hosted_escaping_logical": sum(1 for a in A if a["hosted"] and a.get("logical", {}).get("escapes_all_sectors")),
                           "layers": A}
    # (B) entangling across-block layers on copies of the code
    B = []
    if entangling and across_block and code.n <= 40:
        tri = blocks(code, 3)
        L = rgt.TransversalLayer.across_block_ccz(code.n, 1, code.d)
        block_of = [0] * code.n + [1] * code.n + [2] * code.n
        if tri.k_genuine <= 6:
            B.append(_layer_summary(tri, L, "support" if code.is_css else "cosets", block_of))
        # across-block 2-body escaper (the m=2 catalog's canonical entangling escaper)
        two = blocks(code, 2)
        T2 = np.array([[0, 0, 0, 0], [0, 2, 0, 10], [0, 0, 8, 8], [0, 10, 8, 10]])
        L2 = rgt.TransversalLayer.across_block_pairs(code.n, T2, code.d, "across2_escaper")
        if two.k_genuine <= 6:
            B.append(_layer_summary(two, L2, "support" if code.is_css else "cosets", [0] * code.n + [1] * code.n))
    rep["across_block"] = B
    rep["seconds"] = round(time.time() - t0, 2)
    return rep


def best_protocols(reports: list[dict]) -> list[dict]:
    """Rows that are PROTOCOLS: hosted + transversal/across-block + escaping logical + certified d >= 2."""
    out = []
    for r in reports:
        if not r.get("distance") or r["distance"] < 2:
            continue
        for a in r["single_qudit"]["layers"]:
            if a["hosted"] and a["ft_protocol"] and a.get("logical", {}).get("escapes_all_sectors"):
                out.append({"code": r["label"], "n": r["n"], "k": r["k"], "d": r["distance"], "gamma": r["gamma"],
                            "gate": a["label"], "kind": "T-like (escaping logical)"})
        for b in r["across_block"]:
            if b["hosted"] and b["ft_protocol"] and b.get("logical", {}).get("level") == 3:
                out.append({"code": r["label"], "n": r["n"], "k": r["k"], "d": r["distance"], "gamma": r["gamma"],
                            "gate": b["label"], "kind": "entangling across-block",
                            "logical_escapes": b["logical"]["escapes_all_sectors"]})
    return sorted(out, key=lambda x: (x["gamma"] if x["gamma"] else 99))


# --------------------------------------------------------------------------- candidate sources
def candidates_E2(rng) -> list[tuple[Z4Code, str]]:
    from .lifts import lift_families, z_refinements
    out = []
    for c in lift_families(m_values=(4, 5), rng=rng):
        out.append((c, "lift"))
        for r in z_refinements(c, random.Random(int(rng.integers(1 << 30))), n_samples=6):
            out.append((r, "lift+Zref"))
    return out


def candidates_E3(rng, m: int = 3, n_points: int = 4, trials: int = 300, k_min: int = 2) -> list[tuple[Z4Code, str]]:
    from .galoisring import gr, random_C0, random_symmetric_matrix
    R = gr(m)
    prng = random.Random(int(rng.integers(1 << 30)))
    out = []
    for t in range(trials):
        C0 = random_C0(R, n_points, R.m * n_points - 2 * k_min, prng)
        M = random_symmetric_matrix(R, n_points, prng)
        rows, N = R.graph_code(C0, M)
        if rows.shape[0] == 0:
            continue
        code = Z4Code.from_symplectic(rows, N, 4, label=f"GR(4,{m}) graph n={n_points} t={t}")
        if not code.isotropic or code.k_genuine < k_min:
            continue
        out.append((code, "GR graph"))
    return out


def constrained_sample(n: int, k_target: int, q, rng, max_tries: int = 2000) -> Z4Code | None:
    """E4: draw X-generators with q-weight 0 (S1), accept only if S2/S3 hold; grow the support
    V >= M_X with random vectors keeping Delta_a Phi = 0 (uniform layer q) until |V| = |M_X| 4^k;
    set M_Z = V^perp."""
    from .ringlinalg import module_size, right_kernel_howell, in_span
    q = np.asarray(q, dtype=np.int64)
    layer = rgt.TransversalLayer.uniform(n, q)
    n_gens = max(1, n - k_target - 1)
    for _ in range(max_tries):
        gens = []
        ok = True
        for _ in range(n_gens):
            for _ in range(200):
                a = np.array([rng.randrange(4) for _ in range(n)], dtype=np.int64)
                if not np.any(a) or rgt.q_weight(q, a) != 0:
                    continue
                if in_span(a, gens, 4) if gens else False:
                    continue
                trial = gens + [list(map(int, a))]
                if all(rgt.pair_defect(q, trial[i], trial[j]) == 0 for i in range(len(trial)) for j in range(i, len(trial))) and \
                   all(rgt.triple_defect(q, trial[i], trial[j], trial[l]) == 0 for i in range(len(trial)) for j in range(i, len(trial)) for l in range(j, len(trial))):
                    gens = trial
                    break
            else:
                ok = False
                break
        if not ok or not gens:
            continue
        MX = np.array(gens, dtype=np.int64)
        # grow support
        V = [list(map(int, g)) for g in gens]
        target = module_size(V, 4) * 4 ** k_target
        grown = 0
        for _ in range(400):
            if module_size(V, 4) >= target:
                break
            x = np.array([rng.randrange(4) for _ in range(n)], dtype=np.int64)
            if in_span(x, V, 4):
                continue
            cand = V + [list(map(int, x))]
            if all(rgt.delta_vanishes_on_coset(layer, a, np.zeros(n, dtype=np.int64), cand)[0] for a in gens):
                V = cand; grown += 1
        if module_size(V, 4) < target:
            continue
        MZ = right_kernel_howell(np.array(V, dtype=np.int64), 4)
        MZ = np.array([z for z in MZ if np.any(np.array(z) % 4)], dtype=np.int64).reshape(-1, n)
        try:
            code = Z4Code.from_css(MX, MZ, n, 4, label=f"constrained n={n} k={k_target}")
        except AssertionError:
            continue
        if code.k_genuine >= 1:
            return code
    return None


def worker_E4(seed: int, trials: int, n: int = 12, k_target: int = 2, d_max: int = 3) -> list[dict]:
    rng = random.Random(seed)
    out = []
    for q in (PT.CARRIER_Q,):
        for _ in range(trials):
            code = constrained_sample(n, k_target, q, rng, max_tries=20)
            if code is None:
                continue
            rep = evaluate(code, d_max=d_max, single_tables=[q], entangling=False, provenance=f"E4 seed={seed}")
            out.append(rep)
    return out


# --------------------------------------------------------------------------- CLI
def run(exp: str, out: str, seed: int = 0, jobs: int | None = None, trials: int = 200, n: int = 12,
        k: int = 2, d_max: int = 4, mem_gb: float | None = 4.0) -> dict:
    if mem_gb:
        try:
            import resource
            cap = int(mem_gb * (1 << 30))
            resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
        except Exception:
            pass
    rng = np.random.default_rng(seed)
    t0 = time.time()
    reports = []
    if exp == "E2":
        for code, prov in candidates_E2(rng):
            reports.append(evaluate(code, d_max=d_max, provenance=prov))
    elif exp == "E3":
        for code, prov in candidates_E3(rng, m=3, n_points=n, trials=trials, k_min=k):
            reports.append(evaluate(code, d_max=min(d_max, 3), provenance=prov, entangling=False))
    elif exp == "E4":
        from .parallel import run_chunks, default_jobs
        nj = jobs or default_jobs()
        chunks = run_chunks(lambda s, t: worker_E4(s, t, n, k, min(d_max, 3)), nj, max(1, trials // nj), n_jobs=nj,
                            base_seed=seed * 7919 + 1)
        for ch in chunks:
            reports.extend(ch)
    else:
        raise ValueError(exp)
    catalog = {"exp": exp, "seed": seed, "seconds": round(time.time() - t0, 1), "n_codes": len(reports),
               "protocols": best_protocols(reports), "reports": reports}
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w") as f:
        json.dump(catalog, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return catalog


def main(argv=None):
    ap = argparse.ArgumentParser(description="diverse-gate ring MSD search")
    ap.add_argument("--exp", required=True, choices=["E2", "E3", "E4"])
    ap.add_argument("--out", default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--k", type=int, default=2)
    ap.add_argument("--dmax", type=int, default=4)
    ap.add_argument("--mem-gb", type=float, default=4.0)
    a = ap.parse_args(argv)
    out = a.out or f"runs/diverse/catalog_{a.exp}.json"
    cat = run(a.exp, out, a.seed, a.jobs, a.trials, a.n, a.k, a.dmax, a.mem_gb)
    print(f"{a.exp}: {cat['n_codes']} codes in {cat['seconds']}s -> {out}")
    for p in cat["protocols"][:10]:
        print("  PROTOCOL", p)


if __name__ == "__main__":
    main()
