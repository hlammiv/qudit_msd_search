"""Extract the colour / grafted / final stabilizer codes of Gidney et al.'s end-to-end cultivation circuit
(Zenodo 10.5281/zenodo.13777072, circuits.zip, c=end2end-inplace-distillation,...,d1=5,...,d2=15.stim) by back-propagating
every measured ancilla through its round's Clifford tableau.  Writes graft_d5.json (unsigned); graft_signs_solve.py then
assigns the Z_{2^k} signs -> data/graft_d5_z4.json.

python -m primepower_msd.graft_extract path/to/circuit.stim out.json  [tick_colour tick_graft tick_final]
Round windows for the published d1=5 circuit: colour (59, 68), grafted (89, 98), final (139, 148)."""
from __future__ import annotations

import collections
import json
import sys

import stim

NOISE = {"DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR", "Z_ERROR", "DETECTOR", "OBSERVABLE_INCLUDE", "QUBIT_COORDS", "SHIFT_COORDS"}


def ticks_of(circuit):
    ticks = [[]]
    for op in circuit:
        if op.name == "TICK":
            ticks.append([])
        else:
            ticks[-1].append(op)
    return ticks


def targets(op):
    return [t.value for t in op.targets_copy()]


def round_stabs(ticks, t_reset, t_meas):
    """Stabilizers measured by the ancillas measured at tick t_meas, as data Pauli strings right after the resets of
    tick t_reset.  Ancilla components must be consistent with the reset basis (asserted)."""
    sub = stim.Circuit()
    for t in range(t_reset, t_meas):
        for op in ticks[t]:
            if op.name in ("CX", "S", "S_DAG", "H", "CZ", "SQRT_X", "SQRT_X_DAG"):
                sub.append(op)
            elif op.name not in NOISE and op.name not in ("R", "RX"):
                raise RuntimeError(f"unexpected {op.name} at tick {t}")
    inv = stim.Tableau.from_circuit(sub, ignore_noise=True).inverse()
    resets = {q: op.name for op in ticks[t_reset] if op.name in ("R", "RX") for q in targets(op)}
    meas = {a: op.name for op in ticks[t_meas] if op.name in ("M", "MX") for a in targets(op)}
    out = []
    for a, kind in meas.items():
        P = inv.z_output(a) if kind == "M" else inv.x_output(a)
        s = str(P)[1:]
        supp = {q: ch for q, ch in enumerate(s) if ch != "_"}
        for q, ch in supp.items():
            if q in meas:
                assert q in resets and ((resets[q] == "R" and ch == "Z") or (resets[q] == "RX" and ch == "X")), (a, q, ch)
        out.append(dict(kind=kind, anc=a, data={q: ch for q, ch in supp.items() if q not in meas}))
    return out, resets


def main():
    path, out = sys.argv[1], sys.argv[2]
    w = [int(x) for x in sys.argv[3:9]] if len(sys.argv) > 3 else [59, 68, 89, 98, 139, 148]
    c = stim.Circuit.from_file(path).flattened()
    ticks = ticks_of(c)
    coords = {int(k): list(v) for k, v in c.get_final_qubit_coordinates().items()}
    cc, _ = round_stabs(ticks, w[0], w[1])
    g, rg = round_stabs(ticks, w[2], w[3])
    fn, _ = round_stabs(ticks, w[4], w[5])
    data_g = set(q for v in g for q in v["data"])
    fresh = {int(q): rg[q] for q in data_g if q in rg}
    for name, lst in (("colour", cc), ("grafted", g), ("final", fn)):
        wts = collections.Counter(("".join(sorted(set(v["data"].values()))), len(v["data"])) for v in lst)
        print(name, len(lst), sorted(wts.items()))
    print("fresh", collections.Counter(fresh.values()))
    json.dump(dict(coords=coords, colour=cc, grafted=g, final=fn, fresh=fresh), open(out, "w"))


if __name__ == "__main__":
    main()
