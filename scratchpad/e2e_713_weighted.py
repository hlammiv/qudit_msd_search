"""[[7,1,3]] s = 1 check gadget: leading-order (W = 3) output-error coefficient with a NON-uniform input error distribution.

meas_sim.leading_order assumes every input / T-gate error is Z^e with e uniform on 1..d-1.  The escape's measured
failures are not uniform (MEAS_DISTILL 5k): with the cultivation proxy convention a residual frame (a', b') shifts the
check eigenvalue by e = b' - sigma a' (mod d), sigma = 1.  Here every fault value e carries weight w(e) = (d-1) pi(e)
(pi = measured distribution of e given a failure; uniform pi gives w = 1 and reproduces leading_order).  The T gates of
the check are teleported from the same inputs, so their faults have the same pi.
  input term : input Z^e, outcome cancelled by a layer-1 Z^a (a = -u e) + layer-2 Z^-a on one inner qudit: w(e) w(a) w(-a)
  inner term : weight-3 inner logical z split over the two T layers with zero outcome shift: prod w(z_i)
c_W = (input + inner) / (d-1)^3, so p_out = c_W p^3 with p the total input error probability.
    python scratchpad/e2e_713_weighted.py d pi_1,...,pi_{d-1}      (or 'uniform')
"""
import sys

import numpy as np

from primepower_msd.meas_distill import signed_rm
from primepower_msd.meas_sim import Protocol, leading_order, min_weight_logicals


def coef(d, pi):
    S, tau = signed_rm(1, 3)
    P = Protocol(S, tau, d, np.array([[0]]))
    pi = np.asarray(pi, float); pi = pi / pi.sum()
    w = np.zeros(d); w[1:] = (d - 1) * pi
    u = int(P.u[0])
    A = sum(w[e] * w[(-u * e) % d] * w[(u * e) % d] for e in range(1, d)) * P.n
    _, zs = min_weight_logicals(P, 3)
    B = 0.0
    for z in zs:
        supp = np.nonzero(z)[0]; vals = z[supp]
        for m in range(1 << len(supp)):
            if sum(int(vals[i]) for i in range(len(supp)) if m >> i & 1) % d == 0:
                B += np.prod([w[int(v)] for v in vals])
    return dict(d=d, input=A, inner=B, c=(A + B) / (d - 1) ** 3, ref=leading_order(P, 3)["per_output_coef"])


if __name__ == "__main__":
    d = int(sys.argv[1])
    pi = np.ones(d - 1) if sys.argv[2] == "uniform" else np.array([float(x) for x in sys.argv[2].split(",")])
    print(coef(d, pi))
