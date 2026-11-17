"""total_cost.py at p = 7 with the single-species paired-cap codes added (leading-order surrogates from
census_paired_w8.py, L_8 = 599 / 252) next to the proved sphere curves (total_cost_w8.py).  Gates per round for the paired
codes: the same one-shot grid encoder model as the sphere codes, stage_gates(7, 4, 6) (an over-count: only the slices
x1 = 5, 6 are used).  All raw magic states counted as 1 (for the paired codes this is also the U-only count).
usage: PYTHONPATH=.:scratchpad python scratchpad/total_cost_paired.py OUT.json"""
import sys
src = open("scratchpad/total_cost.py").read()
old_big = src[src.index("        big = [(LeadingOrderCurve(\"[[1746"):src.index("    if p == 11:")]
new_big = ("        from qmsd.rigorous_curve import from_census\n        from w8_proved import upgrade\n"
           "        g6 = stage_gates(7, 4, 6)[0]\n"
           "        print('gates per round: sphere', g7, ' paired (grid model)', g6, ' paired (two 3D encoders)', 2 * stage_gates(7, 3, 6)[0] + 2 * 343, flush=True)\n"
           "        big = [(upgrade(from_census('runs/weighted/p7_m4_sphere_1746_270_7_census.json'), 402)[0], g7, 'rigorous'),\n"
           "               (LeadingOrderCurve('[[596,90,8]]_7 (paired cap, surrogate)', 7, 596, 90, 8, 599.07), g6, 'surrogate'),\n"
           "               (LeadingOrderCurve('[[606,80,8]]_7 (paired cap, surrogate)', 7, 606, 80, 8, 252.0), g6, 'surrogate')]\n")
src = src.replace(old_big, new_big).replace("for p in (7, 11, 17):", "for p in (7,):")
src = src.replace("[(1e-3, 1e-15), (1e-3, 1e-30)]", "[(1e-3, 1e-15), (1e-3, 1e-30), (1e-2, 1e-30)]")
exec(compile(src, "total_cost.py", "exec"), {"__name__": "__main__"})
