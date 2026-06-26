"""Unit tests for qmsd.field (NOTES sec 2: power-sum identity)."""
import galois
import pytest

from qmsd.field import GFp, field_power_sum


@pytest.mark.parametrize("p", [3, 5, 7])
def test_power_sum_matches_brute_force(p):
    """field_power_sum(p,a) equals the literal sum over F_p of x**a mod p."""
    for a in range(0, 3 * (p - 1) + 2):
        brute = sum(pow(x, a, p) for x in range(p)) % p  # 0^0 == 1 (pow(0,0,p)==1)
        assert field_power_sum(p, a) == brute, (p, a)


@pytest.mark.parametrize("p", [3, 5, 7])
def test_power_sum_identity(p):
    """-1 iff a>0 and (p-1)|a; a==0 -> 0; else 0 (NOTES sec 2)."""
    assert field_power_sum(p, 0) == 0
    for a in range(1, 4 * (p - 1) + 1):
        if a % (p - 1) == 0:
            assert field_power_sum(p, a) == p - 1
        else:
            assert field_power_sum(p, a) == 0


@pytest.mark.parametrize("p", [2, 3, 5, 7, 11])
def test_GFp_is_cached_field(p):
    """GFp returns galois.GF(p) and caches the same class object."""
    F = GFp(p)
    assert F is GFp(p)
    assert F.order == p
    assert issubclass(F, galois.FieldArray)
