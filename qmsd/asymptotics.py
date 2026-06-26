"""Asymptotic gamma analysis via saddle point (NOTES sec 8, Appendix A; eqs 37/45/47).

The asymptotic distillation yield gamma_0(theta) of the distinguished puncture
family (m = 3*alpha, r = alpha(p-1)-1) is governed by a generalized p-ary
entropy H_p(theta) = lim_m (1/m) log_p [m, theta*m]_p.  A steepest-descent
analysis (Appendix A) evaluates H_p via the dominant real saddle xi(theta) of
the contour integral for [m, w]_p, giving the closed forms below.
"""
from __future__ import annotations

import math

import scipy.optimize


def _theta_of_xi(xi, p):
    """theta as a function of the saddle xi (NOTES eq 45, solved for theta).

        theta(xi) = ( (p(xi-1) - xi) xi^p + xi ) / ( (xi-1)(xi^p - 1) ).

    Monotone increasing on (0, inf): theta(0)=0, theta(1)=(p-1)/2 (removable
    singularity), theta(inf)=p-1.
    """
    if abs(xi - 1.0) < 1e-7:
        # Removable 0/0 singularity; the saddle xi=1 corresponds to the
        # entropy peak theta = (p-1)/2.
        return (p - 1) / 2.0
    num = (p * (xi - 1.0) - xi) * xi ** p + xi
    den = (xi - 1.0) * (xi ** p - 1.0)
    return num / den


def saddle_xi(theta, p) -> float:
    """Dominant real root xi>0 of eq 45 (closed forms for p=2,3; else brentq) (NOTES eq 45).

    theta(xi) is monotone increasing from 0 (xi->0) to p-1 (xi->inf), so for
    any theta in (0, p-1) there is a unique positive saddle xi.
    """
    assert p >= 2, "p must be a prime >= 2"
    assert 0.0 < theta < p - 1, "saddle requires 0 < theta < p-1"

    if p == 2:
        # eq 48: xi = theta/(1-theta).
        return theta / (1.0 - theta)
    if p == 3:
        # eq 50: xi = (theta + sqrt(3(2-theta)theta + 1) - 1) / (2(2-theta)).
        return (theta + math.sqrt(3.0 * (2.0 - theta) * theta + 1.0) - 1.0) / (
            2.0 * (2.0 - theta)
        )

    median = (p - 1) / 2.0
    if abs(theta - median) < 1e-12:
        return 1.0

    f = lambda xi: _theta_of_xi(xi, p) - theta

    if theta < median:
        # root in (0, 1)
        lo, hi = 1e-14, 1.0 - 1e-12
    else:
        # root in (1, inf); grow upper bracket until residual turns positive
        lo, hi = 1.0 + 1e-12, 2.0
        while f(hi) < 0.0:
            hi *= 2.0
            if hi > 1e12:
                raise RuntimeError("saddle_xi: failed to bracket root")
    return scipy.optimize.brentq(f, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=200)


def H_p(theta, p) -> float:
    """Entropy H_p(theta) = log_p((xi^p-1)/(xi-1)) - theta*log_p(xi) (NOTES eq 47).

    Coefficient on log_p(xi) is theta (the printed eq 47 'theta/3' is a typo).
    For p=2 this is the binary entropy expressed in bits (log base 2).
    Boundary values: H_p(0)=0, H_p((p-1)/2)=1 (peak), H_p->0 as theta->p-1.
    """
    assert p >= 2, "p must be a prime >= 2"
    if theta <= 0.0:
        return 0.0
    if theta >= p - 1:
        return 0.0
    xi = saddle_xi(theta, p)
    lnp = math.log(p)
    # (xi^p - 1)/(xi - 1) = 1 + xi + ... + xi^{p-1}; stable at xi->1.
    geom = sum(xi ** j for j in range(p))
    return (math.log(geom) - theta * math.log(xi)) / lnp


def gamma0(theta, p) -> float:
    """Asymptotic gamma_0(theta), two-branch split at theta=(p-1)/6 (NOTES eq 37).

        gamma_0(theta) = 3*(1 - H_p(theta))                 if theta < (p-1)/6
                       = 3*(1 - H_p(theta)) / H_p(3*theta)   if theta > (p-1)/6

    Both branches agree at the boundary theta=(p-1)/6 where 3*theta hits the
    entropy peak (p-1)/2 and H_p(3*theta)=1.  Requires 3*theta < p-1.
    """
    assert p >= 2, "p must be a prime >= 2"
    assert 3.0 * theta < p - 1, "gamma0 requires 3*theta < p-1"
    assert theta > 0.0, "gamma0 requires theta > 0"
    boundary = (p - 1) / 6.0
    top = 3.0 * (1.0 - H_p(theta, p))
    if theta < boundary:
        return top
    return top / H_p(3.0 * theta, p)


def optimal_gamma(p) -> tuple:
    """Minimize gamma0 over theta in (0,(p-1)/3); return (gamma0_min, t0=theta_opt/(p-1)) (NOTES sec 8, Table 1)."""
    assert p >= 2, "p must be a prime >= 2"
    upper = (p - 1) / 3.0
    res = scipy.optimize.minimize_scalar(
        lambda theta: gamma0(theta, p),
        bounds=(1e-9, upper - 1e-12),
        method="bounded",
        options={"xatol": 1e-12},
    )
    theta_opt = float(res.x)
    return (float(res.fun), theta_opt / (p - 1))
