"""total_cost.py at p = 7 with the three single-species paired-cap codes as PROVED curves (exact weight-8 census +
paired_heavier.py weight-9/10 bounds + generic tail, as in proved_paired_cost.py) next to the proved sphere curve.
Gates per round for the paired codes from env GATES:
  grid : stage_gates(7, 4, 6) (one-shot 4D grid encoder; over-count, only slices x1 = 5, 6 used)
  circuit : the explicit, simulated and verified encoder of paired_encoder.py (runs/coverage/paired_encoder.json "best":
            RM(6,3) on slice 5 + RM(5,3) on slice 6 + 343 CSUMs = 6 058)
usage: GATES=grid|circuit PYTHONPATH=.:scratchpad python scratchpad/total_cost_paired_proved.py OUT.json"""
import sys, os
src = open("scratchpad/total_cost.py").read()
old_big = src[src.index("        big = [(LeadingOrderCurve(\"[[1746"):src.index("    if p == 11:")]
new_big = ("        import json as _j\n        from qmsd.rigorous_curve import from_census, RigorousBoundCurve\n        from w8_proved import upgrade\n"
           "        from paired_heavier import heavier_bounds\n"
           "        gm = os.environ.get('GATES', 'grid')\n"
           "        g6 = stage_gates(7, 4, 6)[0] if gm == 'grid' else _j.load(open('runs/coverage/paired_encoder.json'))['best']\n"
           "        print('gates per round: sphere', g7, ' paired', gm, g6, flush=True)\n"
           "        big = [(upgrade(from_census('runs/weighted/p7_m4_sphere_1746_270_7_census.json'), 402)[0], g7, 'rigorous')]\n"
           "        for _f in ('Y47a', 't45', 't40'):\n"
           "            _c = _j.load(open(f'runs/weighted/pc_census_w8_{_f}_exact.json')); assert _c['part_ii_complete']\n"
           "            big.append((RigorousBoundCurve(f\"[[{_c['n']},{_c['k']},8]]_7 (paired cap, proved)\", 7, _c['n'], _c['k'], 8, float(_c['L8_found']), int(_c['A8_found']),\n"
           "                        heavier=heavier_bounds(343 - _c['k'] // 2, _c['dprime_words'])), g6, 'rigorous'))\n")
src = src.replace(old_big, new_big).replace("for p in (7, 11, 17):", "for p in (7,):")
src = src.replace("[(1e-3, 1e-15), (1e-3, 1e-30)]", "[(1e-3, 1e-15), (1e-3, 1e-30), (1e-2, 1e-30)]")
exec(compile("import os\n" + src, "total_cost.py", "exec"), {"__name__": "__main__"})
