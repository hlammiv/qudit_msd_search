"""Method 2: rigorous structured-enumeration validator for the m=7 cap qutrit code.

d(code) = min over nonzero c in RM_3(9,7) of |supp(c) \\ S|.  This package lower-bounds
d by enumerating structured low-weight codeword families and minimising |supp \\ S| via
an exact additive coset reduction (numba hot loop), distributed across blocks/nodes.
"""
from .capcode import load_cap, points, verify_build
from .enum_families import FAMILIES, family_min_punct, merge_block_maxima

__all__ = [
    "load_cap", "points", "verify_build",
    "FAMILIES", "family_min_punct", "merge_block_maxima",
]
