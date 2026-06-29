"""Structured RM_3(9,7) low-weight codeword families and their rigorous min |supp \\ S|.

A "family" is a set of codewords whose supports are unions of ``j`` parallel ``w``-flats
(j distinct cosets of a common dim-``w`` direction V), all genuine RM_3(9,7) codewords:

  - ``w18``  : w=2, j=2, weight 18 -- the MINIMUM-weight codewords (Kasami-Tokura /
               DGM: support = two parallel 2-flats inside a common 3-flat).  This family
               is the complete weight-18 class.
  - ``w27``  : w=3, j=1, weight 27 -- indicator of an affine 3-flat (degree 8 <= 9).

For such a family |supp(c) cap S| is additive over the j cosets, so

    max over the family of |supp cap S|  =  max over directions V of (top-j coset
    occupancies of S under V)            =  max_topj(...).max()

and the rigorous minimum punctured weight over the family is

    min |supp \\ S|  =  weight  -  max_topj.max().

The full enumeration size (e.g. 99463 * C(243,2) = 2,924,210,889 weight-18 supports) is
never materialised: the additive reduction is exact and visits only the
gaussian_binomial(7,w) directions.  A brute path over explicit supports exists in
``selftest`` purely as an independent cross-check on a sample.

Codewords NOT of pure parallel-flat form (the Leducq second/third-weight classes,
18 < weight <= 45 with non-flat support) are documented as the extension and are covered
probabilistically by Method 1 (Stern).  See RUNBOOK.md.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import subspaces as sub
from .kernels import max_topj

P = 3


@dataclass(frozen=True)
class Family:
    name: str
    w: int        # dimension of the direction subspace V
    j: int        # number of parallel cosets unioned
    weight: int   # codeword weight = j * 3**w

    @property
    def codim(self) -> int:
        return sub.M - self.w

    @property
    def n_directions(self) -> int:
        return sub.gaussian_binomial(sub.M, self.w)


FAMILIES = {
    "w18": Family("w18", w=2, j=2, weight=18),
    "w27": Family("w27", w=3, j=1, weight=27),
}


def directions(fam: Family) -> np.ndarray:
    """Syndrome maps (n_directions, codim, 7) for the family's direction dimension."""
    return sub.syndrome_maps(fam.w)


def family_block_max(fam: Family, Spts: np.ndarray, Hs: np.ndarray,
                     start: int = 0, count: int | None = None) -> dict:
    """Rigorous max |supp cap S| over the family, restricted to directions
    [start, start+count).  Returns the block max + witness (argmax direction index).

    Splitting on ``start/count`` is the unit of work distributed across nodes/threads;
    the global result is the max over all blocks.
    """
    nsub = Hs.shape[0]
    if count is None:
        count = nsub - start
    stop = min(start + count, nsub)
    chunk = np.ascontiguousarray(Hs[start:stop])
    out = max_topj(chunk, Spts, fam.codim, fam.j)
    loc = int(np.argmax(out))
    return {
        "family": fam.name,
        "weight": fam.weight,
        "start": start,
        "count": stop - start,
        "block_max_inter": int(out[loc]),       # max |supp cap S| in this block
        "witness_dir": start + loc,             # global direction index
        "block_min_punct": fam.weight - int(out[loc]),
    }


def family_min_punct(fam: Family, Spts: np.ndarray, Hs: np.ndarray | None = None) -> dict:
    """Whole-family rigorous min |supp \\ S| (single process, all directions)."""
    if Hs is None:
        Hs = directions(fam)
    res = family_block_max(fam, Spts, Hs, 0, Hs.shape[0])
    return {
        "family": fam.name,
        "weight": fam.weight,
        "n_directions": int(Hs.shape[0]),
        "n_supports": int(fam.n_directions) * _coset_pairs(fam),
        "max_inter": res["block_max_inter"],
        "min_punct": fam.weight - res["block_max_inter"],
        "witness_dir": res["witness_dir"],
    }


def _coset_pairs(fam: Family) -> int:
    """Number of coset-choices per direction (for reporting the full support count)."""
    from math import comb
    ncoset = P ** fam.codim
    return comb(ncoset, fam.j)


def merge_block_maxima(blocks: list[dict]) -> dict:
    """Merge per-block results (Method-2 merge rule = MAX of block maxima ->
    MIN of punctured weights).  Returns the family-wide rigorous result."""
    assert blocks, "no blocks to merge"
    fam = blocks[0]["family"]
    weight = blocks[0]["weight"]
    best = max(blocks, key=lambda b: b["block_max_inter"])
    return {
        "family": fam,
        "weight": weight,
        "max_inter": best["block_max_inter"],
        "min_punct": weight - best["block_max_inter"],
        "witness_dir": best["witness_dir"],
        "n_blocks": len(blocks),
        "directions_covered": sum(b["count"] for b in blocks),
    }
