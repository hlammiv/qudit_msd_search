"""Single-machine driver: rigorous min |supp \\ S| over all enumerated structured
families of RM_3(9,7), for the m=7 cap qutrit code.  Prints the Method-2 verdict.

Run:  python -m cap_validate.structured.run_all
"""
from __future__ import annotations

import json

from . import capcode as cc
from . import enum_families as ef


def run():
    cap = cc.load_cap()
    Spts = cap["cap_pts"]
    results = {}
    for name, fam in ef.FAMILIES.items():
        results[name] = ef.family_min_punct(fam, Spts)
    rigorous_min = min(r["min_punct"] for r in results.values())
    verdict = {
        "method": "structured_enum",
        "code": "[[1968,219]]_3 (m=7 cap qutrit)",
        "families": results,
        "rigorous_min_punct_over_enumerated": rigorous_min,
        "claim": (
            f"No RM_3(9,7) codeword in the enumerated structured families punctures "
            f"below {rigorous_min}. The complete weight-18 class gives exactly "
            f"{results['w18']['min_punct']} (= boundary for d>=10). Non-flat Leducq "
            f"classes (18<wt<=45) are covered probabilistically by Method 1 (Stern)."
        ),
        "consistent_with_d_ge_10": rigorous_min >= 10,
    }
    print(json.dumps(verdict, indent=2, default=int))
    return verdict


if __name__ == "__main__":
    run()
