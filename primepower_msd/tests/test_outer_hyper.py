import numpy as np

from primepower_msd.outer_hyper import (brute_sensitive, certify, check_bits_from, construct2, random_biregular)


def test_certifier_matches_bruteforce_random():
    rng = np.random.default_rng(3)
    n_cases = 0
    for _ in range(200):
        n = int(rng.choice([12, 15, 18])); s = int(rng.choice([2, 3])); w = 3
        if (n * s) % w:
            continue
        bc, m = random_biregular(n, s, w, rng)
        if any(len(set(r)) < s for r in bc):
            continue
        cb = check_bits_from(bc, m, w)
        for D in (5, 7):
            assert (certify(bc, cb, D) < 0) == brute_sensitive([list(r) for r in bc], D)
            n_cases += 1
    assert n_cases > 20


def test_constructed_outer_is_sensitive():
    bc, m, _ = construct2(40, 3, 4, 7, seed=0, max_rounds=60)
    assert certify(bc, check_bits_from(bc, m, 4), 7) < 0
    assert brute_sensitive([list(r) for r in bc], 7)
