"""Magic-state distillation cost and output-error formulas (NOTES eqs 38/39).

A triorthogonal CSS code [[n,k,d]]_p consumes n noisy magic states with input
depolarizing error ``delta_in`` and post-selects on the X-syndrome to output k
higher-fidelity magic states.  The single-round figures of merit are:

  * survival / yield        nbar_T = (1 - (p-1)/p*delta_in)**n * k   (NOTES eq 39)
  * distillation cost       C      = n / nbar_T                       (NOTES sec 7)
  * per-output infidelity   delta_out^(i) = A_d^(i)/((p-1)p^(d-1)) * delta_in**d
                                                                       (NOTES eq 38)
  * averaged infidelity     delta_out  = A_d/(nbar_T (p-1) p^(d-1))
                                          * delta_in**d
                                          * (1 - (p-1)/p*delta_in)**(n-d)
                                                                       (NOTES eq 39)

Load-bearing correction (NOTES sec 7, "Critical correction"): the nbar_T
exponent is **n**, not n/k; the in-eq39 survival factor exponent is **(n-d)**.
"""
from __future__ import annotations


def nbar_T(n, k, p, delta_in) -> float:
    """Expected number of accepted outputs (NOTES eq 39).

        nbar_T = (1 - (p-1)/p * delta_in)**n * k

    The survival probability is raised to the power n (the number of consumed
    input states), NOT n/k -- this is the corrected exponent that reproduces the
    paper's worked numbers (NOTES sec 7, "Critical correction").
    """
    assert p >= 2, "p must be a prime >= 2"
    assert n >= 0 and k >= 0, "n, k must be non-negative"
    assert 0.0 <= delta_in <= 1.0, "delta_in must lie in [0, 1]"
    survival = 1.0 - (p - 1) / p * delta_in
    return survival ** n * k


def cost(n, nbar_T_value) -> float:
    """Distillation cost C = n / nbar_T -- inputs per accepted output (NOTES sec 7).

    Suggested as a better single-round search objective than gamma once both
    delta_in and delta_out are specified.
    """
    assert nbar_T_value > 0.0, "nbar_T must be positive to define a cost"
    return n / nbar_T_value


def delta_out_avg(n, k, d, A_d, p, delta_in) -> float:
    """Averaged output infidelity per output qudit (NOTES eq 39).

        delta_out ~ A_d / (nbar_T (p-1) p^(d-1))
                    * delta_in**d
                    * (1 - (p-1)/p * delta_in)**(n-d)

    with nbar_T = (1 - (p-1)/p*delta_in)**n * k.  The survival factor here carries
    exponent (n - d), the leading-order correction for the d failed locations.
    """
    assert p >= 2, "p must be a prime >= 2"
    assert d >= 1, "distance d must be >= 1"
    assert 0.0 <= delta_in <= 1.0, "delta_in must lie in [0, 1]"
    nbar = nbar_T(n, k, p, delta_in)
    assert nbar > 0.0, "nbar_T must be positive"
    survival = 1.0 - (p - 1) / p * delta_in
    denom = nbar * (p - 1) * p ** (d - 1)
    return A_d / denom * delta_in ** d * survival ** (n - d)


def delta_out_per_output(A_d_i, d, p, delta_in) -> float:
    """Output infidelity for a single output qudit i (NOTES eq 38).

        delta_out^(i) = A_d^(i) / ((p-1) p^(d-1)) * delta_in**d ,   i = 1,...,k

    where A_d^(i) is the number of weight-d logical-Z operators acting on output
    i (the per-output subset of the total A_d used in eq 39).
    """
    assert p >= 2, "p must be a prime >= 2"
    assert d >= 1, "distance d must be >= 1"
    assert 0.0 <= delta_in <= 1.0, "delta_in must lie in [0, 1]"
    denom = (p - 1) * p ** (d - 1)
    return A_d_i / denom * delta_in ** d
