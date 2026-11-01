import numpy as np
from primepower_msd.host_census import gadget, census, faults_depolarizing, faults_logical, schedule, floor, q_round, q_gate


def test_depolarizing_census_matches_cliff_census():
    P, C, inputs = gadget(4)
    N, NZ, W = census(C, inputs, 4, faults_depolarizing(C, 4))
    assert abs(sum(NZ.values()) / P.n_in - 0.357) < 0.002
    assert abs(sum(N.values()) / P.n_in - 2.92) < 0.02


def test_logical_census_exposures():
    P, C, inputs = gadget(4)
    F, T_end = faults_logical(C, inputs, 4, 0)
    N, NZ, W = census(C, inputs, 4, F)
    assert T_end > 20 and abs(W["gate"] - 2 * sum(1 for op in C.ops if op[0] == "CSUM")) < 1e-6
    assert abs(NZ["gate"] - 4.667) < 0.01 and abs(NZ["idle"] - 5.333) < 0.01
    assert NZ["fourier"] == 0 and NZ["meas"] == 0 and NZ["prep"] == 0   # only Z-type gate and idle faults survive the check
    assert N["idle_x"] > NZ["idle"]                                        # X-type faults are counted by the frame model but projected


def test_floor_monotone():
    f = [floor(1e-3, D) for D in (15, 19, 23, 27)]
    assert all(a > b for a, b in zip(f, f[1:])) and floor(1e-3, 23) > floor(1e-3, 23, gate_model="round")
    assert q_gate(1e-3, 23) > q_round(1e-3, 23)
