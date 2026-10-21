"""Order-D outer codes with s checks per input (HHPW Lemma 4.4 condition  2|Mv| + |v| >= D  for every Z_d).

Inputs = bits with s distinct checks; checks have degree w.  For an error set V, |Mv| >= U(V) := #checks touched
exactly once (a leaf check sees one nonzero value), and with adversarial Z_d values every other check can cancel,
so the alphabet-independent condition is  2 U(V) + |V| >= D  for all 1 <= |V| <= D-1.

Certifier (complete for violators).  DFS from a root bit b0 (= min index of V): keep the checks touched once
("open"); take the smallest open check c and either FREEZE it (it stays unique; at most u_max(t) =
floor((D-1-t)/2) frozen) or add one of c's other bits (index > b0).  When nothing is open, V is closed up to the
frozen checks and is tested.  Any violator V* yields along the consistent path either V* or a closed V subset of
V* with the same frozen checks and smaller t, which is itself a violator, so "no violator found" is a proof.

Construction: random biregular configuration model, then edge switches that break every violator found
(local re-certification around the switched bits), then a full certification.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange


@njit
def _dfs(b0, V, t, deg, frozen, nfrozen, bit_checks, check_bits, D, require_min, out):
    # V[:t] current bits; deg[c] = #V-bits on check c; frozen[c] = 1 if frozen unique
    s = bit_checks.shape[1]
    # find smallest open check (deg == 1, not frozen) among checks of V
    best = -1
    U = 0
    O = 0
    for i in range(t):
        for j in range(s):
            c = bit_checks[V[i], j]
            if deg[c] == 1:
                U += 1
                if frozen[c] == 0:
                    O += 1
                    if best == -1 or c < best:
                        best = c
    if best == -1:
        if 2 * U + t < D:
            for i in range(t):
                out[i] = V[i]
            out[D - 1] = t
            return True
        return False
    if t == D - 1:
        return False
    umax_t = (D - 1 - t) // 2
    # each further bit closes at most s open checks; remaining freezes close one each
    if O > s * (D - 1 - t) + max(0, umax_t - nfrozen):
        return False
    # branch 1: freeze
    if nfrozen + 1 <= umax_t:
        frozen[best] = 1
        if _dfs(b0, V, t, deg, frozen, nfrozen + 1, bit_checks, check_bits, D, require_min, out):
            frozen[best] = 0
            return True
        frozen[best] = 0
    # branch 2: add a bit through `best`
    w = check_bits.shape[1]
    for k in range(w):
        b = check_bits[best, k]
        if b < 0:
            continue
        if require_min and b <= b0:
            continue
        inV = False
        for i in range(t):
            if V[i] == b:
                inV = True
                break
        if inV:
            continue
        ok = True
        for j in range(s):                       # cannot touch a frozen check
            if frozen[bit_checks[b, j]] == 1:
                ok = False
        if not ok:
            continue
        V[t] = b
        for j in range(s):
            deg[bit_checks[b, j]] += 1
        found = _dfs(b0, V, t + 1, deg, frozen, nfrozen, bit_checks, check_bits, D, require_min, out)
        for j in range(s):
            deg[bit_checks[b, j]] -= 1
        if found:
            return True
    return False


@njit
def violator_at(b0, bit_checks, check_bits, D, require_min):
    m = check_bits.shape[0]
    deg = np.zeros(m, np.int32)
    frozen = np.zeros(m, np.int8)
    V = np.zeros(D, np.int64)
    out = -np.ones(D, np.int64)
    V[0] = b0
    for j in range(bit_checks.shape[1]):
        deg[bit_checks[b0, j]] += 1
    _dfs(b0, V, 1, deg, frozen, 0, bit_checks, check_bits, D, require_min, out)
    return out


@njit
def certify(bit_checks, check_bits, D):
    """-1 if sensitive, else a root bit of a violator."""
    for b0 in range(bit_checks.shape[0]):
        o = violator_at(b0, bit_checks, check_bits, D, True)
        if o[D - 1] > 0:
            return b0
    return -1


def check_bits_from(bit_checks, m, w):
    cb = -np.ones((m, w), np.int64)
    fill = np.zeros(m, np.int64)
    for b, cs in enumerate(bit_checks):
        for c in cs:
            cb[c, fill[c]] = b; fill[c] += 1
    return cb


def brute_sensitive(bit_checks, D):
    """exhaustive over all bit subsets of size <= D-1 (tiny instances only)."""
    import itertools as it
    n = len(bit_checks)
    for t in range(1, D):
        for V in it.combinations(range(n), t):
            cnt = {}
            for b in V:
                for c in bit_checks[b]:
                    cnt[c] = cnt.get(c, 0) + 1
            U = sum(1 for v in cnt.values() if v == 1)
            if 2 * U + t < D:
                return False
    return True


def swap_edges(bc, cb, b, j, b2, j2):
    """swap check bc[b,j] <-> bc[b2,j2] and update check_bits in place."""
    c, c2 = bc[b, j], bc[b2, j2]
    cb[c, np.nonzero(cb[c] == b)[0][0]] = b2
    cb[c2, np.nonzero(cb[c2] == b2)[0][0]] = b
    bc[b, j], bc[b2, j2] = c2, c


def random_biregular(n_bits, s, w, rng):
    m = n_bits * s // w
    assert m * w == n_bits * s
    stubs = np.repeat(np.arange(m), w)
    rng.shuffle(stubs)
    bc = stubs.reshape(n_bits, s).copy()
    return bc, m


def construct(n_bits, s, w, D, seed=0, max_switch=2_000_000, verbose=False):
    """Random (s,w) graph + violator-breaking switches.  Returns (bit_checks, m, switches) or None."""
    rng = np.random.default_rng(seed)
    bc, m = random_biregular(n_bits, s, w, rng)
    # first remove repeated checks within a bit
    def fix_repeats():
        for b in range(n_bits):
            while len(set(bc[b])) < s:
                j = rng.integers(s); b2 = rng.integers(n_bits); j2 = rng.integers(s)
                bc[b, j], bc[b2, j2] = bc[b2, j2], bc[b, j]
    fix_repeats()
    cb = check_bits_from(bc, m, w)
    switches = 0
    dirty = list(range(n_bits))
    while switches < max_switch:
        bad = -1
        while dirty:
            b = dirty.pop()
            o = violator_at(b, bc, cb, D, False)
            if o[D - 1] > 0:
                bad = o[:o[D - 1]]; break
        if isinstance(bad, int):
            r = certify(bc, cb, D)
            if r < 0:
                return bc, m, switches
            o = violator_at(r, bc, cb, D, True); bad = o[:o[D - 1]]
        # switch one edge of a random violating bit with a random edge elsewhere
        b = int(bad[rng.integers(len(bad))]); j = int(rng.integers(s))
        b2 = int(rng.integers(n_bits)); j2 = int(rng.integers(s))
        c, c2 = bc[b, j], bc[b2, j2]
        if c == c2 or c2 in bc[b] or c in bc[b2]:
            continue
        swap_edges(bc, cb, b, j, b2, j2)
        switches += 1
        dirty.extend([b, b2] + [int(x) for x in cb[c] if x >= 0] + [int(x) for x in cb[c2] if x >= 0])
        if verbose and switches % 1000 == 0:
            print("switches", switches, flush=True)
    return None


@njit
def short_cycle_bit(b, bit_checks, check_bits):
    """True if bit b lies on a Tanner 4-cycle (shares >= 2 checks with a bit) or 6-cycle (bits x, y adjacent to b
    via different checks that share a check other than those)."""
    s = bit_checks.shape[1]; w = check_bits.shape[1]
    for j1 in range(s):
        c1 = bit_checks[b, j1]
        for k1 in range(w):
            x = check_bits[c1, k1]
            if x < 0 or x == b:
                continue
            for j in range(s):                         # 4-cycle: x shares another check with b
                cx = bit_checks[x, j]
                if cx != c1:
                    for jj in range(s):
                        if bit_checks[b, jj] == cx:
                            return True
            for j2 in range(j1 + 1, s):
                c2 = bit_checks[b, j2]
                for k2 in range(w):
                    y = check_bits[c2, k2]
                    if y < 0 or y == b or y == x:
                        continue
                    for jx in range(s):                # 6-cycle: x and y share a check not in {c1, c2}
                        cx = bit_checks[x, jx]
                        if cx == c1 or cx == c2:
                            continue
                        for jy in range(s):
                            if bit_checks[y, jy] == cx:
                                return True
    return False


def remove_short_cycles(bc, m, w, rng, max_switch=10_000_000, verbose=False):
    """Edge switches until the Tanner graph has girth >= 8."""
    n, s = bc.shape
    cb = check_bits_from(bc, m, w)
    bad = [b for b in range(n) if short_cycle_bit(b, bc, cb)]
    sw = 0
    while bad and sw < max_switch:
        b = bad[rng.integers(len(bad))]
        if not short_cycle_bit(b, bc, cb):
            bad.remove(b); continue
        j = int(rng.integers(s)); b2 = int(rng.integers(n)); j2 = int(rng.integers(s))
        c, c2 = bc[b, j], bc[b2, j2]
        if c == c2 or c2 in bc[b] or c in bc[b2]:
            continue
        old = short_cycle_bit(b2, bc, cb)
        swap_edges(bc, cb, b, j, b2, j2)
        sw += 1
        if short_cycle_bit(b, bc, cb) or (not old and short_cycle_bit(b2, bc, cb)):
            swap_edges(bc, cb, b, j, b2, j2)             # reject: did not help / created a new one
            continue
        for x in set([b2] + [int(v) for v in cb[c] if v >= 0] + [int(v) for v in cb[c2] if v >= 0]):
            if short_cycle_bit(x, bc, cb) and x not in bad:
                bad.append(x)
        if verbose and sw % 5000 == 0:
            print("  girth-8 switches", sw, "bad", len(bad), flush=True)
    return bc, sw, len(bad)


@njit(parallel=True)
def scan_roots(bit_checks, check_bits, D):
    """bool per bit: a violator with this bit as minimum exists (full certificate iff none)."""
    n = bit_checks.shape[0]
    bad = np.zeros(n, np.bool_)
    for b0 in prange(n):
        o = violator_at(b0, bit_checks, check_bits, D, True)
        bad[b0] = o[D - 1] > 0
    return bad


def construct2(n_bits, s, w, D, seed=0, max_rounds=200, verbose=False):
    """Random (s,w) graph; repeat: parallel scan -> for every violating root, switch an edge of a bit in its
    violator until that root is clean (local re-check of the two touched bits); stop when a scan is clean."""
    rng = np.random.default_rng(seed)
    bc, m = random_biregular(n_bits, s, w, rng)
    for b in range(n_bits):
        while len(set(bc[b])) < s:
            j = rng.integers(s); b2 = rng.integers(n_bits); j2 = rng.integers(s)
            bc[b, j], bc[b2, j2] = bc[b2, j2], bc[b, j]
    cb = check_bits_from(bc, m, w)
    switches = 0
    for rnd in range(max_rounds):
        bad = np.nonzero(scan_roots(bc, cb, D))[0]
        if verbose:
            print(f"round {rnd}: violating roots {len(bad)}, switches so far {switches}", flush=True)
        if len(bad) == 0:
            return bc, m, switches
        for b0 in rng.permutation(bad):
            for _ in range(200):
                o = violator_at(int(b0), bc, cb, D, False)
                if o[D - 1] <= 0:
                    break
                V = o[:o[D - 1]]
                b = int(V[rng.integers(len(V))]); j = int(rng.integers(s))
                b2 = int(rng.integers(n_bits)); j2 = int(rng.integers(s))
                c, c2 = bc[b, j], bc[b2, j2]
                if c == c2 or c2 in bc[b] or c in bc[b2]:
                    continue
                swap_edges(bc, cb, b, j, b2, j2)
                switches += 1
                if violator_at(b2, bc, cb, D, False)[D - 1] > 0:     # do not create a violator at b2
                    swap_edges(bc, cb, b, j, b2, j2); switches -= 1
    return None
