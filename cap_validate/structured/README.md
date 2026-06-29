# cap_validate.structured — Method 2 (structured_enum)

Rigorous structured-enumeration validator for the m=7 cap qutrit code `[[1968,219]]_3`.

`d = min_{c in RM_3(9,7), c != 0} |supp(c) \ S|`, S = the 219 cap puncture points.
We lower-bound `d` by minimising `|supp(c) \ S|` over the structured low-weight
codeword families of RM_3(9,7) (Kasami-Tokura / DGM: unions of `j` parallel `w`-flats).

Key idea (makes the 2.9e9-support enumeration trivial and *exact*): for a family whose
support is a union of `j` parallel `w`-flats, `|supp ∩ S|` is **additive** over the `j`
cosets, so the family-wide extremum is, per direction subspace `V`, the **top-j coset
occupancy of S**. We only visit `gaussian_binomial(7,w)` directions (numba `prange` hot
loop in `kernels.max_topj`), never the full support list.

| file | role |
|------|------|
| `capcode.py` | load cap S, F_3^7 points, rebuild/verify G0 (55×1968) |
| `subspaces.py` | enumerate dim-`w` directions (RREF) → syndrome/coset maps |
| `kernels.py` | numba hot loops: `max_topj` (additive), `inter_counts` (brute cross-check) |
| `enum_families.py` | families `w18`,`w27`; block max; merge rule |
| `distribute.py` | block split / per-node worker / merge (no network) |
| `mock_cluster.py` | LOCAL two-node split == single-process (+resume) |
| `selftest.py` | build, structure, count, kernel==brute, minima |
| `run_all.py` | single-machine Method-2 verdict |
| `RUNBOOK.md` | parameterized 2-machine (local + lenore:60022) launch |

## Result

`w18` (complete min-weight class, 2,924,510,589 supports) → **min |supp\S| = 10**;
`w27` (affine 3-flats) → 18. Rigorous min over enumerated families = **10**, consistent
with `d >= 10` (γ < 1). The weight-18 class is exact/complete; non-flat Leducq classes
(18<wt≤45) are covered probabilistically by Method 1 (Stern) — see RUNBOOK note.

```bash
cd /home/hlamm/Desktop/QC/prime_msd
python -m cap_validate.structured.selftest
python -m cap_validate.structured.mock_cluster
python -m cap_validate.structured.run_all
```
