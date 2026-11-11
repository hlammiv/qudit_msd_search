"""species_sensitivity.py with the rigorous weight-8 bound for [[1746,270,7]]_7.  Writes runs/weighted/species_sensitivity_w8.json."""
src = open("scratchpad/species_sensitivity.py").read().replace("runs/weighted/species_sensitivity.json", "runs/weighted/species_sensitivity_w8.json")
src = src.replace('c = from_census("runs/weighted/p7_m4_sphere_1746_270_7_census.json")',
                  'c = __import__("w8_proved").upgrade(from_census("runs/weighted/p7_m4_sphere_1746_270_7_census.json"), 402)[0]')
assert "upgrade" in src
exec(compile(src, "species_sensitivity.py", "exec"), {"__name__": "__main__"})
