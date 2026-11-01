"""Numbers for NOTES.md: (1) the check space of the k=2290 output set (rank of RM(13,4) on O,
dim I(O)_{<=13} = 2320 - 2290 = 30 < 35 = dim RM(3,4), so O lies on no hypersurface of degree
<= 10); (2) the random-code expectation of kernel words of weight w for the 30-check code;
(3) the design table for O = Z(g): cap on k, theorem bound, gamma, for deg g = 1..9."""
import json, math
import numpy as np
from qmsd.reedmuller import rm_generator, d_rm
from hyp_common import rank_modp

p, m, r = 11, 4, 13
d = json.load(open("runs/rm_m4/p11_r13_k2290_seed1.json"))
O = np.array(d["puncture_columns_1indexed"]) - 1
G = np.asarray(rm_generator(r, m, p), dtype=np.int64)
K = G.shape[0]
h = rank_modp(G[:, O], p)
print(f"dim RM(13,4) = {K}; rank on O (k={len(O)}) = {h}; checks = {K - h}; dim RM(3,4) = {rm_generator(3, m, p).shape[0]}, dim RM(2,4) = {rm_generator(2, m, p).shape[0]}")

n = p ** m - len(O); rho = K - len(O)
print("random-code model E[#projective kernel words of weight w] = C(n,w)(p-1)^(w-1)/p^rho, n =", n, "rho =", rho)
for w in range(3, 10):
    E = math.comb(n, w) * (p - 1) ** (w - 1) / p ** rho
    print(f"  w={w}: {E:.3e}")

print("\nDesign table for O = Z(g), deg g = delta (unweighted RM_11(13,4) on F_11^4):")
rt = m * (p - 1) - r - 1
for delta in range(1, 10):
    cap = K - rm_generator(r - delta, m, p).shape[0]
    b = d_rm(rt + delta, m, p)
    gam = math.log((p ** m - cap) / cap) / math.log(b)
    print(f"  delta={delta}: cap k <= {cap}, d >= d(RM({rt+delta},4)) = {b}, gamma at k=cap: {gam:.3f}, "
          f"typical |Z| ~ {p**(m-1)} -> gamma {math.log((p**m - p**(m-1)) / p**(m-1)) / math.log(b):.3f}")
