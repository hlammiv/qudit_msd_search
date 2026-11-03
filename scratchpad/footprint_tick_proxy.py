"""Footprint with idle noise at every tick (proxy until the tick CSUM runs land): the single-patch fits of
runs/surface/zd_d4_tick{f}.jsonl replace HOST_Z4, and the CSUM / two-patch fits are rescaled by the same factors per basis
(A ratio, p_th ratio).  python scratchpad/footprint_tick_proxy.py f [cult=a,p_L] [footprint args...]
(f = 0: no host patch; cult=a,p_L replaces the 1e-3 cultivation input, e.g. the end-to-end cultivation -> escape numbers;
 e2e5=attempts_per_kept,p_out replaces cultivation and escape at 5e-4 by the end-to-end numbers)"""
import sys

from primepower_msd import footprint as FP, host_census as HC

TICK = {  # f: (both-basis A, p_th), X basis (A, p_th) from surface_fit 4_tick{f}; untick X basis = (0.049, 7.2e-3) (zd_d4.jsonl)
    1.0: ((0.088, 5.42e-3), (0.0921, 5.46e-3)),
    0.1: ((0.054, 7.12e-3), (0.0591, 7.07e-3)),
}
f = float(sys.argv.pop(1))
if len(sys.argv) > 1 and sys.argv[1].startswith("cult="):
    a_c, p_c = map(float, sys.argv.pop(1)[5:].split(","))
    FP.CULT_D4[FP.CULT_KEY][1e-3] = (a_c, p_c)
if len(sys.argv) > 1 and sys.argv[1].startswith("e2e5="):  # e2e5=attempts,p_out: measured cultivation -> escape at 5e-4
    att, p_c = map(float, sys.argv.pop(1)[5:].split(","))   # (escape discard and error are inside both numbers)
    FP.CULT_D4[FP.CULT_KEY][5e-4] = (1 / att, p_c)
    FP.ESCAPE_Z4_GAP[5e-4] = (0.0, 0.0)
(Ab, pb), (Ax, px) = TICK.get(f, ((FP.HOST_Z4["A"], FP.HOST_Z4["p_th"]), (0.049, 7.2e-3)))
ra_b, rp_b = Ab / FP.HOST_Z4["A"], pb / FP.HOST_Z4["p_th"]
ra_x, rp_x = Ax / 0.049, px / 7.2e-3
FP.HOST_Z4.update(A=Ab, p_th=pb)
for D_ in (FP.CSUM_Z4, FP.BASE2_Z4):
    D_.update(A=D_["A"] * ra_b, p_th=D_["p_th"] * rp_b)
for D_ in (HC.BASE2_X, HC.CSUM_X):
    D_.update(A=D_["A"] * ra_x, p_th=D_["p_th"] * rp_x)
print(f"tick f={f}: HOST {FP.HOST_Z4}, CSUM {FP.CSUM_Z4}, BASE2_X {HC.BASE2_X}, CSUM_X {HC.CSUM_X}")
import textwrap
src = open(FP.__file__).read().split('if __name__ == "__main__":\n', 1)[1]
exec(textwrap.dedent(src), FP.__dict__)                 # footprint's CLI block, run against the patched constants
