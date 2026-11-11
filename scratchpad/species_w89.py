"""species_cost.py with the rigorous weight-8 and weight-9 bounds (scratchpad/w8_proved.py) for the d = 7 sphere codes.
Writes runs/weighted/species_cost_w89.json."""
import sys
sys.argv = ["species_cost.py"]
src = open("scratchpad/species_cost.py").read().replace("runs/weighted/species_cost.json", "runs/weighted/species_cost_w89.json")
src = src.replace("if __name__ == \"__main__\":\n    main()", "")
g = {"__name__": "species_w8"}
exec(compile(src, "species_cost.py", "exec"), g)
from w8_proved import upgrade
_orig = g["from_census"]
BLOCK = {"[[1746,270,7]]_7": 402, "[[1756,260,7]]_7": 412}
def patched(path, name=None):
    c = _orig(path, name)
    return upgrade(c, BLOCK[c.name])[0] if c.name in BLOCK and c.d == 7 else c
g["from_census"] = patched
g["main"]()
