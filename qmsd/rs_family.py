"""Closed-form exact performance of the Reed-Solomon (m = 1) triorthogonal family.

For m = 1 the triorthogonal space is RM(r, 1) over F_p = the extended Reed-Solomon code
[p, r+1, p-r] (evaluations of polynomials of degree <= r at every field element), triorthogonal
iff 3r < p - 1.  Puncturing k points gives the magic code [[p-k, k, r-k+2]]_p and every space the
performance engine needs is again MDS, so its weight enumerator is closed-form:

* A = rowspan(G0) = RS shortened at the k punctures: MDS [n, r+1-k, p-r], n = p - k;
* span(A, h_i) (words vanishing at the other k-1 punctures, restricted to the n physical
  coordinates) = RS shortened at k-1 points and punctured at puncture i: MDS [n, r+2-k, p-r-1];
  the k scalar cosets A + t h_i share one histogram, (W_span - W_A)/(p-1);
* H = span(A, h_1..h_k) = RS punctured at all k punctures: MDS [n, r+1, p-r-k].

With the MDS weight enumerator  A_w = C(n,w) sum_{j=0}^{w-d} (-1)^j C(w,j) (q^{w-d+1-j} - 1)
these feed ``profile_from_histograms`` directly, so eps_out, P_acc and the threshold are exact
for every member of the family at every prime (no enumeration; validated against the enumerating
engine in tests/test_rs_family.py).
"""
from __future__ import annotations

from math import comb

from .performance import DistillationProfile, profile_from_histograms


def r_max_rs(p: int) -> int:
    """Largest r with 3r < p - 1 (triorthogonality of RM(r, 1))."""
    return (p - 2) // 3


def mds_weight_enumerator(n: int, K: int, q: int) -> list[int]:
    """Weight distribution of an [n, K, n-K+1]_q MDS code (K = 0 gives the zero code)."""
    if K <= 0:
        return [1] + [0] * n
    d = n - K + 1
    A = [0] * (n + 1)
    A[0] = 1
    for w in range(d, n + 1):
        s = 0
        for j in range(0, w - d + 1):
            s += (-1) ** j * comb(w, j) * (q ** (w - d + 1 - j) - 1)
        A[w] = comb(n, w) * s
    assert sum(A) == q ** K, "MDS enumerator does not sum to q^K"
    return A


def rs_profile(p: int, k: int, r: int | None = None, block: bool = True) -> DistillationProfile:
    """Exact DistillationProfile of the RS magic code [[p-k, k, r-k+2]]_p (default r = r_max)."""
    r = r_max_rs(p) if r is None else r
    assert 3 * r < p - 1, "RM(r,1) is triorthogonal only for 3r < p - 1"
    assert 1 <= k <= r, "need 1 <= k <= r (d = r - k + 2 >= 2)"
    n = p - k
    Aw = mds_weight_enumerator(n, r + 1 - k, p)
    Wspan = mds_weight_enumerator(n, r + 2 - k, p)
    coset = [(s - a) // (p - 1) for s, a in zip(Wspan, Aw)]
    assert all((s - a) % (p - 1) == 0 for s, a in zip(Wspan, Aw))
    Hw = mds_weight_enumerator(n, r + 1, p) if block and r + 1 <= n else None
    return profile_from_histograms(p, r + 1 - k, Aw, cosets=[coset] * k, Hw=Hw)


def rs_parameters(p: int, k: int, r: int | None = None) -> tuple[int, int, int]:
    """(n, k, d) of the RS magic code."""
    r = r_max_rs(p) if r is None else r
    return p - k, k, r - k + 2
