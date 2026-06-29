"""M2b — general CSS codes (M_X != M_Z): does decoupling X- and Z-stabilizers admit a distance>=2
single-ququart distillation code?

The M2 (weakly-self-dual) null came from M_X = M_Z forcing self-orthogonality, which fights distance.
In a GENERAL CSS code the X-stabilizers (which control transversality) and the Z-stabilizers (which
provide distance) decouple.  A diagonal gate D is transversal when M_X lies in the gate's translation
stabilizer Stab(phi); the Z-stabilizers M_Z subseteq M_X^perp are then free to supply distance.

We first report a sharp STRUCTURAL fact (which gates can be transversal on ANY nontrivial code), then
search Stab(phi)-compatible general CSS codes for a distance>=2, k=1, cyclic, level-3 witness.

Run:  python -m primepower_msd.m2b_search
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from primepower_msd.ring_css import build_css_general, logical_group_invariants
from primepower_msd.ring_transversal import extract_phase, stab1, stab_elements, certify_distillation_code
from primepower_msd.hierarchy_search import search as level3_family


def structural_report(d: int):
    """For each level-3 gate, report stab1 (single-coordinate translation stabilizer) — its size bounds
    whether the gate can be transversal on any nontrivial code."""
    print("=" * 78)
    print(f"STRUCTURAL: single-coordinate translation stabilizer stab1(phi) of each level-3 gate (d={d})")
    print("=" * 78)
    can, cannot = [], []
    for (a, b, sig, U) in level3_family(d):
        phi, N = extract_phase(U, d)
        s1 = stab1(phi, N, d)
        nontrivial = len(s1) > 1
        (can if nontrivial else cannot).append((a, b, tuple(sorted(s1.keys())), N, sig))
    print(f"  gates with NONTRIVIAL stab1 (CAN be globally transversal): {len(can)}/{d*(d-1)}")
    for a, b, keys, N, sig in can[:6]:
        print(f"    antidiff(S^{a} Z^{b}): stab1={keys} (N={N}), phase exps(2pi/{4*d})={sig}")
    print(f"  gates with TRIVIAL stab1={{0}} (CANNOT be globally transversal on any nontrivial code): {len(cannot)}")
    for a, b, keys, N, sig in cannot[:6]:
        print(f"    antidiff(S^{a} Z^{b}): phase exps(2pi/{4*d})={sig}")
    return can, cannot


def gc_worker(seed: int, trials: int):
    """One parallel worker: `trials` random Stab(phi)-guided general-CSS codes per (gate, n).

    Top-level (picklable) so joblib/loky can ship it to a process.  Returns a picklable summary.
    """
    import numpy as _np
    from primepower_msd.ringlinalg import dual_module
    from primepower_msd.ring_css import build_css_general
    from primepower_msd.ring_transversal import extract_phase, stab_elements, certify_distillation_code
    from primepower_msd.hierarchy_search import search as _level3_family

    d = 4
    B_PER_A = 8                       # sample several Z-stabilizer sets per X-stabilizer set (reuse dual A)
    rng = _np.random.default_rng(seed)
    fam = _level3_family(d)
    best = None
    trivial = 0
    scanned = 0
    for gi, (ga, gb, sig, U) in enumerate(fam):
        phi, N = extract_phase(U, d)
        for n in (5, 7):
            stab = stab_elements(phi, N, d, n)
            if not stab or len(stab) <= 1:
                continue
            stab_arr = [_np.array(s) for s in stab if any(s)]
            if not stab_arr:
                continue
            for _ in range(trials):
                rx = int(rng.integers(1, min(4, len(stab_arr) + 1)))
                A = _np.array([stab_arr[i] for i in rng.choice(len(stab_arr), size=rx, replace=False)])
                mxperp = dual_module(A, d)               # computed once per A, reused for all B
                if len(mxperp) <= 1:
                    continue
                mxperp_arr = [_np.array(y) for y in mxperp if any(y)]
                for _ in range(B_PER_A):
                    rz = int(rng.integers(1, min(n, len(mxperp_arr) + 1)))
                    B = _np.array([mxperp_arr[i] for i in rng.choice(len(mxperp_arr), size=rz, replace=False)])
                    code = build_css_general(A, B, d, mxperp_pre=mxperp)   # k1-prune + dual(A) reuse
                    scanned += 1
                    if code is None or code.k != 1 or not code.cyclic:
                        continue
                    cert = certify_distillation_code(code, U)
                    if not cert.get("is_distillation_code"):
                        continue
                    if code.distance < 2:
                        trivial += 1
                    elif best is None or code.distance > best["distance"]:
                        best = {"n": n, "A": A.tolist(), "B": B.tolist(), "gi": gi,
                                "distance": code.distance, "cert": cert}
    return {"best": best, "trivial": trivial, "scanned": scanned}


def main() -> int:
    import time
    from primepower_msd.parallel import run_chunks, default_jobs

    d = 4
    structural_report(d)
    print()
    n_jobs = default_jobs()
    n_chunks, trials = n_jobs, 1_200
    t0 = time.time()
    print(f"  parallel general-CSS search: {n_jobs} workers x {n_chunks} chunks x {trials} trials/(gate,n) ...",
          flush=True)
    results = run_chunks(gc_worker, n_chunks=n_chunks, trials_per_chunk=trials, n_jobs=n_jobs)
    trivial = sum(r["trivial"] for r in results)
    scanned = sum(r["scanned"] for r in results)
    bests = [r["best"] for r in results if r["best"]]
    best = max(bests, key=lambda b: b["distance"]) if bests else None

    print("=" * 78)
    print("VERDICT (M2b — general CSS)")
    print("=" * 78)
    print(f"  scanned {scanned} general-CSS k=1 cyclic codes in {time.time()-t0:.1f}s "
          f"({n_jobs} cores; Stab(phi)-guided, transversality guaranteed).")
    print(f"  trivial (distance-1) transversal hits: {trivial}")
    if best:
        fam = level3_family(d)
        ga, gb, _, U = fam[best["gi"]]
        code = build_css_general(best["A"], best["B"], d)
        cert = best["cert"]
        print(f"\n  NONTRIVIAL WITNESS: [[{code.n},1,{best['distance']}]]_{d} general CSS, "
              f"transversal antidiff(S^{ga} Z^{gb})")
        print(f"    X-stabilizers (Z_{d}): {best['A']}")
        print(f"    Z-stabilizers (Z_{d}): {best['B']}")
        print(f"    logical orders in M_Z^perp/M_X: {logical_group_invariants(code)} cyclic={code.cyclic}")
        print(f"    induced logical level: {cert['logical_level']} strict-L3: {cert['logical_strict_level3']}")
        print(f"    -> NON-TRIVIAL single-ququart distillation code via GENERAL CSS CONFIRMED.")
        return 0
    print("\n  No distance>=2 witness found via Stab(phi)-guided general CSS.")
    print("  Combined with the structural fact above: the genuinely-quadratic level-3 gates cannot be")
    print("  globally transversal, and the transversable (qubit-like) gates meet only distance-1 codes.")
    print("  Next: per-coset search for the quadratic gates (R1 ring-triorthogonality); M3 distance at larger n.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
