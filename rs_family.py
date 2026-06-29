"""T2 — the optimal small gamma<1 Reed-Solomon (m=1) qudit MSD code family, p >= 23.

At m=1 the quantum code G0^perp is a *punctured MDS* code, so its distance is the closed form
    d = r_max - k + 2 ,   r_max = floor((p-2)/3) ,   n = p - k ,   k logical qudits,
INDEPENDENT of which k columns are punctured (MDS puncture-invariance). Hence the entire gamma<1
family is exactly enumerable with no search and no distance computation. This script enumerates it,
verifies the closed form against the qmsd triorthogonal + MacWilliams tooling, and reports the
per-prime Pareto-optimal codes (min gamma / smallest block n / min cost C).

Result recorded in NEW_CODES.md ("p >= 23 (Reed-Solomon, m=1) ..."). Run: PYTHONPATH=. python rs_family.py
"""
import math
from qmsd.triorthogonal import build_triorthogonal_code
from qmsd.weightdist import exact_distance_and_Ad
from qmsd.reedmuller import rm_generator
from qmsd import distillation as dist

DIN = 1e-2
PRIMES = [23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
rmax = lambda p: (p - 2) // 3


def verify_closed_form(p=23):
    """Confirm d = r_max - k + 2 (and that A_d is finite) against qmsd for a few codes."""
    rm = rmax(p); G = rm_generator(rm, 1, p)
    out = []
    for k in (4, 5, 6):
        b = build_triorthogonal_code(p, 1, rm, list(range(1, k + 1)), G=G)
        r = exact_distance_and_Ad(b['X_stab'], p, max_words=10**7)
        out.append((p - k, k, rm - k + 2, r['distance'], r['B_d']))
    return out


def family(p):
    """All gamma<1 codes [[p-k, k, r_max-k+2]] for prime p, with cost C."""
    rm = rmax(p); fam = []
    for k in range(1, rm + 1):
        n, d = p - k, rm - k + 2
        if d < 2:
            continue
        g = math.log(n / k) / math.log(d)
        if g < 1.0:
            C = dist.cost(n, dist.nbar_T(n, k, p, DIN))
            fam.append((g, n, k, d, C))
    return fam


if __name__ == "__main__":
    print("SANITY (closed-form d vs qmsd, p=23):")
    for n, k, dcf, dq, Ad in verify_closed_form(23):
        print(f"  [[{n},{k},{dcf}]]: qmsd d={dq} A_d={Ad}  {'OK' if dq == dcf else 'MISMATCH'}")
    print("\nMin-gamma gamma<1 RS code per prime:")
    for p in PRIMES:
        fam = family(p)
        if not fam:
            continue
        g, n, k, d, C = min(fam, key=lambda x: x[0])
        print(f"  p={p:>2}: [[{n},{k},{d}]]_{p}  gamma={g:.3f}  C={C:.1f}  ({len(fam)} gamma<1 codes total)")
