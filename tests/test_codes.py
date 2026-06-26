"""Tests for qmsd.codes: the analytic Manhattan family and explicit-puncture reconstruction."""
import math

import pytest

from qmsd.codes import code_from_manhattan, code_from_puncture, gamma, Code
from qmsd.oracle import load_oracle
from tests import ground_truth as gt


# --- Analytic engine: reproduce ALL of Table 2 exactly (closes the coverage gap) ---
@pytest.mark.parametrize("row", gt.TABLE2_SMALLEST_SUBLOG,
                         ids=[f"p{r['p']}m{r['m']}" for r in gt.TABLE2_SMALLEST_SUBLOG])
def test_code_from_manhattan_reproduces_table2(row):
    c = code_from_manhattan(row["p"], row["m"], row["w"])
    assert (c.n, c.k, c.d) == (row["n"], row["k"], row["d"])
    assert c.gamma < 1.0  # every Table-2 row is sublogarithmic
    assert c.d_certified is True


# --- Explicit engine: reconstruct the small oracle codes from their puncture columns ---
SMALL_ORACLE = [oc for oc in load_oracle() if oc.n <= 130]


@pytest.mark.parametrize("oc", SMALL_ORACLE, ids=[oc.label for oc in SMALL_ORACLE])
def test_code_from_puncture_small_oracle(oc):
    c = code_from_puncture(oc.p, oc.m, oc.puncture_columns_1indexed, max_distance=oc.d + 1)
    assert (c.n, c.k, c.d) == (oc.n, oc.k, oc.d)
    assert c.full_rank is True
    assert c.d_certified is True
    row = next(t for t in gt.TABLE3_CODES if (t["p"], t["n"], t["k"], t["d"]) ==
               (oc.p, oc.n, oc.k, oc.d))
    assert c.A_d == row["A_d"]
    assert round(c.gamma, 2) == row["gamma"]


def test_gamma_function_and_property_agree():
    c = code_from_manhattan(5, 16, 16)  # a Table-2 ququint code
    assert math.isclose(gamma(c.n, c.k, c.d), c.gamma)
    # headline ququint code value
    assert round(gamma(519, 106, 5), 2) == 0.99


def test_gamma_nan_when_distance_unknown():
    c = Code(p=5, n=20, k=5, d=None, d_certified=False)
    assert math.isnan(c.gamma)
