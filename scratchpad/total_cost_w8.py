"""total_cost.py at p = 7 with the PROVED curves of the sphere codes: complete weight-7 census + the
rigorous weight-8 bound (scratchpad/w8_proved.py).  usage: total_cost_w8.py OUT.json"""
import sys
src = open("scratchpad/total_cost.py").read()
old_big = src[src.index("        big = [(LeadingOrderCurve(\"[[1746"):src.index("    if p == 11:")]
new_big = ("        from qmsd.rigorous_curve import from_census\n        from w8_proved import upgrade\n"
           "        big = [(upgrade(from_census('runs/weighted/p7_m4_sphere_1746_270_7_census.json'), 402)[0], g7, 'rigorous'),\n"
           "               (upgrade(from_census('runs/weighted/p7_m4_sphere_1756_260_7_census.json'), 412)[0], g7, 'rigorous')]\n")
src = src.replace(old_big, new_big).replace("for p in (7, 11, 17):", "for p in (7,):")
src = src.replace("[(1e-3, 1e-15), (1e-3, 1e-30)]", "[(1e-3, 1e-15), (1e-3, 1e-30), (1e-2, 1e-30)]")
exec(compile(src, "total_cost.py", "exec"), {"__name__": "__main__"})
