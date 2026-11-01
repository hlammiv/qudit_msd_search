"""Scaling table from results.jsonl: seconds per (instance, n, w) and solver.  '>T' = time limit hit,
'OOM' = MITM left table exceeded the 4 GB cap.  usage: summarize.py [results.jsonl]"""
import sys, json, os, collections
f = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "results.jsonl")
rows = [json.loads(l) for l in open(f) if l.strip()]
tab = collections.defaultdict(dict)
solvers = ["mitm", "cpsat-onehot", "cpsat-elem", "cpsat-elem-dfs", "sat-cadical195", "cpsat-sys"]
for r in rows:
    if r.get("tag") == "validate" or r.get("instance") not in ("A", "B", "C"):
        continue
    key = (r["instance"], r["n"], r["w"])
    st = r["status"]
    if st in ("INFEASIBLE", "OPTIMAL") or (st == "FEASIBLE" and r["solver"] == "mitm") or st == "FEASIBLE":
        cell = f"{r['seconds']:.2f}" + (" (word)" if st in ("FEASIBLE", "OPTIMAL") else "")
    else:
        cell = f">{r['seconds']:.0f}"
    name = r["solver"] + (" lin0" if r.get("tag") == "lin0" else "")
    if name not in solvers:
        solvers.append(name)
    prev = tab[key].get(name)
    if prev is None or (prev.startswith(">") and not cell.startswith(">")):
        tab[key][name] = cell
used = [s for s in solvers if any(s in v for v in tab.values())]
print("| inst | n | w | " + " | ".join(used) + " |")
print("|---|---|---|" + "---|" * len(used))
for key in sorted(tab, key=lambda k: (k[0], k[2], k[1])):
    print(f"| {key[0]} | {key[1]} | {key[2]} | " + " | ".join(tab[key].get(s, "") for s in used) + " |")
