"""U^a + U^b -> U^(a+b) by gate teleportation at p = 7 (WEIGHTED_SUPPORTS.md, Species accounting)."""
import math

import numpy as np

P = 7
W = np.exp(2j * np.pi / P)


def magic(a, z=0):
    return np.array([W ** ((a * pow(x, 3, P) + z * x) % P) for x in range(P)]) / math.sqrt(P)


def test_teleport_adds_exponents_and_z_errors():
    for a in range(1, P):
        for b in range(1, P):
            if (a + b) % P == 0:
                continue
            for z_data in (0, 3):
                for z_anc in range(P):
                    psi, anc = magic(a, z_data), magic(b, z_anc)
                    for y in range(P):
                        out = np.array([psi[s] * anc[(s - y) % P] for s in range(P)])
                        out /= np.linalg.norm(out)
                        out *= np.array([W ** ((3 * b * s * s * y - 3 * b * s * y * y) % P) for s in range(P)])
                        target = magic(a + b, z_data + z_anc)
                        assert abs(abs(np.vdot(target, out)) - 1) < 1e-12


def test_error_composition_stays_uniform():
    e1, e2 = 1e-2, 3e-2
    d1 = np.full(P, e1 / (P - 1)); d1[0] = 1 - e1
    d2 = np.full(P, e2 / (P - 1)); d2[0] = 1 - e2
    conv = np.array([sum(d1[a] * d2[(z - a) % P] for a in range(P)) for z in range(P)])
    assert np.ptp(conv[1:]) < 1e-15
    assert abs((1 - conv[0]) - (e1 + e2 - e1 * e2 * P / (P - 1))) < 1e-15
