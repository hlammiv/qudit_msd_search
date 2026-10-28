"""Summarise runs/weighted/plane_sweep/*.jsonl: per k, the d_plane histogram, survivor count
(d_plane >= 7 and no class-block word of weight <= 5), survivor seeds, and timing."""
import json, glob, collections, sys
out = {}
for f in sorted(glob.glob("runs/weighted/plane_sweep/plane_sweep_k*.jsonl")):
    k = int(f.split("_k")[1][:3])
    recs = [json.loads(l) for l in open(f) if l.strip()]
    e = out.setdefault(k, {"draws": 0, "d_plane": collections.Counter(), "survivors": [], "plane_pass_block_fail": [], "cpu_min": 0.0})
    for r in recs:
        e["draws"] += 1; e["d_plane"][r["d_plane"]] += 1; e["cpu_min"] += (r["t_plane"] + r.get("t_block", 0)) / 60
        if r["d_plane"] >= 7:
            (e["survivors"] if r.get("block_d_le5") is None else e["plane_pass_block_fail"]).append(r["seed"])
print("| k | n | draws | d_plane 3 / 4 / 5 / 6 / 7 / 8 | plane pass, block fail | survivors | rate | CPU h |")
print("|---|---|---|---|---|---|---|---|")
for k, e in sorted(out.items()):
    h = " / ".join(str(e["d_plane"].get(d, 0)) for d in range(3, 9))
    print(f"| {k} | {2016 - k} | {e['draws']} | {h} | {len(e['plane_pass_block_fail'])} | {len(e['survivors'])} | {len(e['survivors'])/e['draws']:.0%} | {e['cpu_min']/60:.1f} |")
json.dump({k: {"draws": e["draws"], "d_plane_hist": dict(sorted(e["d_plane"].items())), "survivor_seeds": sorted(e["survivors"]),
              "plane_pass_block_fail_seeds": sorted(e["plane_pass_block_fail"]), "cpu_hours": round(e["cpu_min"] / 60, 2)} for k, e in sorted(out.items())},
          open("runs/weighted/plane_sweep/summary.json", "w"), indent=1)
