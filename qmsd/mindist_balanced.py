"""Memory-rebalanced meet-in-the-middle minimum distance (additive sibling of qmsd.mindist).

WHY THIS MODULE EXISTS
----------------------
``qmsd.mindist.min_dependent_columns`` certifies the exact minimum distance of the
punctured-RM dual codes by a meet-in-the-middle (MITM) column-dependency search, splitting
a weight-``d`` codeword into a LEFT ``a``-subset (hashed) and a disjoint RIGHT ``b``-subset
(streamed), ``a + b = d``.  In the original split the LEFT half carries the FULL nonzero
coefficient freedom ``(p-1)**a`` and the RIGHT half fixes its leading coefficient to 1
(the global F_p* scaling redundancy).  The LEFT table therefore has

    (p-1)**a * C(n, a)   entries.

For the headline pin -- the exact distance of ``[[237,52,>=6]]_17`` (p=17, n=237) -- the
``d = 6`` step uses ``a = b = 3`` and the LEFT table balloons to

    (p-1)**3 * C(237,3) = 4096 * 2.19e6 = 9.0e9 entries  (~72 GB codes, ~287 GB with
    supports)  ->  OOM even on a 125 GB box.

THE REBALANCE
-------------
The global F_p* scaling lets us fix exactly ONE coefficient.  We move that fixed coefficient
from the RIGHT half to the LEFT half: the LEFT ``a``-subset fixes its smallest column's
coefficient to 1 (free ``(p-1)**(a-1)`` for the rest) and the RIGHT ``b``-subset now carries
the FULL ``(p-1)**b`` freedom.  The LEFT table shrinks to

    (p-1)**(a-1) * C(n, a) = (p-1)**2 * C(237,3) = 256 * 2.19e6 = 5.6e8 entries
    (~4.5 GB codes + supports ~ 18 GB total -- fits a 125 GB box),

and the RIGHT stream ``(p-1)**b * C(n,b) = 9.0e9`` is enumerated and parallelised across
processes.  Coverage is unchanged: every weight-6 codeword has a non-empty 3-column left
block, so scaling its left block's smallest column to coefficient 1 is a valid canonical
form that covers it.  The collision logic is identical to ``mindist``: a syndrome match
``s_L == -s_R`` on DISJOINT supports is a genuine weight-``d`` codeword, given no weight
``< d`` codeword exists (which the increasing-``d`` search establishes first).

WHAT IS / ISN'T REBALANCED
--------------------------
Only the ``d = 6`` step (``a = 3``) has the giant LEFT table, so only that step is
rebalanced here.  The ``d <= 5`` steps use ``a <= 2`` whose ORIGINAL LEFT table is already
small (``(p-1)**2 * C(n,2) ~ 7e6`` at the pin), so they are delegated verbatim to
``mindist._weight_d_exists_parallel`` (which already memmaps its small left table and
parallelises the b=3 right stream).  The end-to-end result is identical to
``mindist.min_dependent_columns``; this module only changes HOW the d=6 step is laid out in
memory so it actually fits.

MEMORY SAFETY
-------------
* The ~18 GB LEFT table is built ONCE in the main process and handed to workers as a joblib
  argument; joblib auto-memmaps ndarrays > 1 MB, so all workers mmap a single read-only copy
  rather than each copying it.
* ``n_jobs`` is bounded so ``n_jobs * (per-worker right batch)`` stays well under available
  RAM.
* A RAM GUARD refuses to run (``MemoryError``) and LOGS the projection if the projected LEFT
  table exceeds ``ram_fraction`` of available RAM -- so a 12 GB laptop will decline the 18 GB
  pin instead of OOM-crashing, and the heavy run is steered to the 125 GB box.
* HARD-CAPPED at ``d = 6`` (half size 3): it never constructs a weight-4 ``C(n,4)`` table.

This module is ADDITIVE: it imports ``qmsd.mindist`` helpers but never modifies them.
"""
from __future__ import annotations

import os
import sys
from math import comb

import numpy as np

# Reuse the (read-only) primitives from mindist -- never edit them.
from .mindist import (
    _build_left,  # noqa: F401  (kept for symmetry / potential reuse)
    _encode,
    _pairs_after,
    _powers,
    _weight_d_exists,
    _weight_d_exists_parallel,
)

HARD_CAP = 6  # halves of size <= 3 -> d <= 6; never builds a weight-4 table.
_SUPP_DTYPE = np.int32  # column indices; int32 covers n up to ~2.1e9.
_BYTES_PER_ENTRY = 32  # conservative resident estimate: int64 code (8) + supports + slack.


# ---------------------------------------------------------------------------
# RAM accounting / projection logging
# ---------------------------------------------------------------------------
def _log(msg: str) -> None:
    """Emit to stderr (captured by the detached logfile) so projections are always visible."""
    print(msg, file=sys.stderr, flush=True)


def available_ram_bytes() -> int:
    """Best-effort available RAM in bytes (0 if unknown -> guard is skipped, logged)."""
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        pass
    # Windows (and any other platform without /proc/meminfo or os.sysconf):
    # psutil works cross-platform, including Windows.
    try:
        import psutil
        return psutil.virtual_memory().available
    except ImportError:
        return 0

def left_table_entries(p: int, a: int, n: int) -> int:
    """Number of rebalanced LEFT entries: (p-1)**(a-1) * C(n, a)."""
    return (p - 1) ** (a - 1) * comb(n, a)


def left_table_bytes(p: int, a: int, n: int) -> int:
    """Projected resident bytes of the rebalanced LEFT table (conservative)."""
    return left_table_entries(p, a, n) * _BYTES_PER_ENTRY


def _bound_n_jobs(n_jobs, r: int, n: int, avail: int) -> int:
    """Clamp n_jobs so n_jobs * (per-worker right batch) stays well under RAM.

    A worker's transient peak is dominated by a handful of (r x C(n,2)) int64 syndrome
    buffers for the widest (smallest-i) right slice.  We budget ~6 such buffers and cap the
    pool at 40% of available RAM (the read-only left table is memmapped, not counted here).
    """
    cpu = os.cpu_count() or 1
    if n_jobs in (-1, None):
        n_jobs = cpu
    n_jobs = max(1, min(int(n_jobs), cpu))
    if avail:
        per_worker = max(1, r * comb(n, 2) * 8 * 6)
        cap = max(1, int(0.40 * avail / per_worker))
        n_jobs = min(n_jobs, cap)
    return max(1, n_jobs)


# ---------------------------------------------------------------------------
# Rebalanced LEFT table (a = 3, smallest-column coefficient fixed to 1)
# ---------------------------------------------------------------------------
def _build_left_fixed3(H, p, powers):
    """All 3-subsets (i<j<k) with the SMALLEST column's coeff fixed to 1, c2,c3 free.

    Returns (codes_sorted, supports_sorted): ``codes`` is encode(s_L) (int64), ``supports``
    is an (N, 3) int32 array of sorted column indices, both reordered by ``codes`` so a
    target can be range-looked-up with searchsorted.  N = (p-1)**2 * C(n,3).
    """
    r, n = H.shape
    nz = range(1, p)
    N = left_table_entries(p, 3, n)
    codes = np.empty(N, dtype=np.int64)
    supps = np.empty((N, 3), dtype=_SUPP_DTYPE)
    off = 0
    for i in range(n - 2):
        j, k = _pairs_after(i + 1, n)
        if j.size == 0:
            continue
        Hi = H[:, i] % p  # leading column, coefficient fixed to 1
        block_supp = np.stack(
            [
                np.full(j.size, i, dtype=_SUPP_DTYPE),
                j.astype(_SUPP_DTYPE),
                k.astype(_SUPP_DTYPE),
            ],
            axis=1,
        )
        for c2 in nz:
            Hj = (c2 * H[:, j]) % p
            for c3 in nz:
                S = (Hi[:, None] + Hj + c3 * H[:, k]) % p
                cc = _encode(S, powers)
                m = cc.size
                codes[off : off + m] = cc
                supps[off : off + m] = block_supp
                off += m
    assert off == N, (off, N)
    order = np.argsort(codes, kind="stable")
    return codes[order], supps[order]


# ---------------------------------------------------------------------------
# Rebalanced RIGHT stream (b = 3, FULL nonzero coefficients) + collision
# ---------------------------------------------------------------------------
def _right_full3_collision(H, p, powers, left_codes, left_supps, i_list):
    """Worker: does any weight-3 right-subset (first column in ``i_list``, FULL coeffs)
    collide with the shared read-only LEFT table on disjoint supports?

    Mirrors mindist._right_i_collision but with the extra leading-coefficient loop, because
    in the rebalance it is the RIGHT half (not the left) that carries full coefficients.
    """
    r, n = H.shape
    nz = range(1, p)
    for i in i_list:
        j, k = _pairs_after(i + 1, n)
        if j.size == 0:
            continue
        supp = np.stack([np.full(j.size, i), j, k], axis=1)
        for c1 in nz:
            Hi = (c1 * H[:, i]) % p
            for c2 in nz:
                Hj = (c2 * H[:, j]) % p
                for c3 in nz:
                    S = (Hi[:, None] + Hj + c3 * H[:, k]) % p
                    target = _encode((p - S) % p, powers)  # encode(-s_R)
                    pos = np.searchsorted(left_codes, target, side="left")
                    in_range = pos < left_codes.size
                    hit = np.zeros(target.shape, dtype=bool)
                    hit[in_range] = left_codes[pos[in_range]] == target[in_range]
                    if not hit.any():
                        continue
                    for ridx in np.nonzero(hit)[0]:
                        t = target[ridx]
                        lo = np.searchsorted(left_codes, t, side="left")
                        hi = np.searchsorted(left_codes, t, side="right")
                        rcols = set(int(x) for x in supp[ridx])
                        for li in range(lo, hi):
                            if rcols.isdisjoint(int(x) for x in left_supps[li]):
                                return True
    return False


def _weight6_exists_balanced(H, p, powers, n_jobs, ram_fraction) -> bool:
    """True iff H has a weight-6 column dependency, via the rebalanced (a=3) MITM.

    Assumes (as the increasing-d search guarantees) that no weight < 6 dependency exists, so
    a disjoint-support syndrome collision is a genuine weight-6 codeword.
    """
    r, n = H.shape
    if 6 > n:
        return False

    entries = left_table_entries(p, 3, n)
    proj = left_table_bytes(p, 3, n)
    avail = available_ram_bytes()
    _log(
        f"[mindist_balanced d=6] rebalanced LEFT table projection: "
        f"(p-1)^2 * C({n},3) = {(p - 1) ** 2} * {comb(n, 3):,} = {entries:,} entries "
        f"~ {proj / 1e9:.2f} GB resident; available RAM ~ {avail / 1e9:.2f} GB "
        f"(guard fraction {ram_fraction:.2f})."
    )
    if avail and proj > ram_fraction * avail:
        raise MemoryError(
            f"projected rebalanced LEFT table ~{proj / 1e9:.2f} GB exceeds "
            f"{ram_fraction:.2f} * available {avail / 1e9:.2f} GB; refusing to build it "
            f"(run the pin on a larger-RAM box, e.g. lenore 125 GB). Entries={entries:,}."
        )

    left_codes, left_supps = _build_left_fixed3(H, p, powers)
    if left_codes.size == 0:
        return False

    n_jobs = _bound_n_jobs(n_jobs, r, n, avail)
    _log(
        f"[mindist_balanced d=6] LEFT table built: {left_codes.size:,} entries, "
        f"{left_codes.nbytes / 1e9:.2f} GB codes + {left_supps.nbytes / 1e9:.2f} GB supports; "
        f"streaming RIGHT (p-1)^3 * C({n},3) = {(p - 1) ** 3 * comb(n, 3):,} across "
        f"n_jobs={n_jobs}."
    )

    if n_jobs == 1:
        return _right_full3_collision(
            H, p, powers, left_codes, left_supps, list(range(n - 2))
        )

    from joblib import Parallel, delayed

    # round-robin i-slices balance load (small i has the widest (j,k) range -> most work).
    chunks = [list(range(w, n - 2, n_jobs)) for w in range(n_jobs)]
    res = Parallel(n_jobs=n_jobs)(
        delayed(_right_full3_collision)(H, p, powers, left_codes, left_supps, ch)
        for ch in chunks
        if ch
    )
    return any(res)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def min_dependent_columns_balanced(
    H, p, d_max=None, n_jobs=-1, *, ram_fraction=0.5
) -> int:
    """Minimum number of F_p-linearly-dependent columns of H (= distance of ker H).

    Identical RESULT to ``mindist.min_dependent_columns`` (exact, certified up to
    ``HARD_CAP = 6``), but the ``d = 6`` step uses the memory-rebalanced layout above so it
    fits in RAM at the ``[[237,52,>=6]]_17`` pin where the original ``d=6`` left table OOMs.

    Steps d <= 5 are delegated to ``mindist._weight_d_exists_parallel`` (small left table,
    already memmapped + parallel).  Step d == 6 uses ``_weight6_exists_balanced`` with the
    shared-memmap rebalanced left table, a parallel right stream, and the RAM guard.

    Parameters
    ----------
    H : array_like, shape (r, n), entries reduced mod p (the X_stab / G0 matrix).
    p : int prime.
    d_max : int or None.  Search cap (default r+1), clamped to ``min(n, HARD_CAP)``.
    n_jobs : int.  Worker processes for the b=3 streams; -1 = all cores (then clamped to a
        RAM-safe bound for the d=6 step).
    ram_fraction : float.  Refuse the d=6 step if the projected left table exceeds this
        fraction of available RAM.

    Raises
    ------
    MemoryError : if the projected d=6 left table is too big for this machine (logged).
    ValueError : if no dependency of size <= cap exists (distance exceeds the MITM cap).
    """
    H = np.asarray(H, dtype=np.int64) % p
    r, n = H.shape
    if r == 0:
        return 1  # empty parity check: every single column is a codeword
    powers = _powers(p, r)
    cap = min(d_max if d_max is not None else (r + 1), n, HARD_CAP)
    if n_jobs in (-1, None):
        n_jobs = os.cpu_count() or 1

    for d in range(1, cap + 1):
        if d == 1:
            if (H == 0).all(axis=0).any():
                return 1
            continue
        if d <= 5:
            # a <= 2: ORIGINAL split's left table is already small -> reuse mindist verbatim.
            if n_jobs == 1:
                found = _weight_d_exists(H, p, d, powers)
            else:
                found = _weight_d_exists_parallel(H, p, d, powers, n_jobs)
            if found:
                return d
        else:  # d == 6 (a = b = 3): the rebalanced, memory-safe step.
            if _weight6_exists_balanced(H, p, powers, n_jobs, ram_fraction):
                return 6

    raise ValueError(
        f"no dependent column set of size <= {cap} found; minimum distance exceeds {cap} "
        f"(MITM certifies up to {HARD_CAP}) -- not certified"
    )
