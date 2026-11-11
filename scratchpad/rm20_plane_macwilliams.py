"""Independent check of A_8, A_9 of RM_7(20,4) (rm20_w8_count.py, rm20_wW_count.py found every weight-8 and
weight-9 word planar): the weight-w words supported in one affine plane are the weight-w words of the dual of the
plane-cubic evaluation code, a [49, 10] code over F_7 (bivariate cubics restricted from 4 variables are all
bivariate cubics).  Enumerate all 7^10 codewords of that code, MacWilliams-transform to the dual, and compare
139 650 * B_w with the frame counts (a weight >= 8 word cannot lie on a line, so it lies in exactly one plane).
Writes runs/weighted/rm20_plane_macwilliams.json.  ~10 min, < 1 GB."""
import json, time, itertools
from math import comb
import numpy as np

p, q = 7, 7
pts = np.array([(x, y) for x in range(p) for y in range(p)], dtype=np.int64)
mons = [(a, b) for a in range(4) for b in range(4 - a)]
G = np.array([[pow(int(x), a, p) * pow(int(y), b, p) % p for x, y in pts] for a, b in mons], dtype=np.int64)
assert G.shape == (10, 49)
from qmsd.structured_ad import _rank
assert _rank(G % p, p) == 10
T0 = time.time()
L = np.array(list(itertools.product(range(p), repeat=5)), dtype=np.int64)
CL = (L @ G[:5]) % p; CR = (L @ G[5:]) % p                      # 16807 x 49 each
hist = np.zeros(50, dtype=np.int64)
for i in range(len(CL)):
    W = ((CR + CL[i]) % p != 0).sum(axis=1)
    hist += np.bincount(W, minlength=50)
    if i % 2000 == 0: print(f"[{time.time()-T0:5.0f}s] {i}/{len(CL)}", flush=True)
A = [int(v) for v in hist]; n = 49
assert sum(A) == p ** 10 and A[0] == 1
def K(w, j):
    return sum((-1) ** i * (q - 1) ** (w - i) * comb(j, i) * comb(n - j, w - i) for i in range(0, w + 1))
B = {}
for w in range(0, 12):
    s = sum(A[j] * K(w, j) for j in range(n + 1))
    assert s % (p ** 10) == 0
    B[w] = s // p ** 10
planes = 139650
frame = {8: json.load(open("runs/weighted/rm20_w8_count.json"))["A8_RM20"], 9: json.load(open("runs/weighted/rm20_w9_count_v2.json"))["A9_RM20"]}
out = {"code_weight_distribution": A, "dual_low_weights": B, "planar_A8": planes * B[8], "planar_A9": planes * B[9],
       "frame_A8": frame[8], "frame_A9": frame[9], "agree": planes * B[8] == frame[8] and planes * B[9] == frame[9]}
# weights 5-7 of the plane dual are line words: each line of the plane (56) carries 21 / 14 / 22 projective words
out["line_check"] = {w: (B[w], 56 * c * 6) for w, c in ((5, 21), (6, 14), (7, 22))}
json.dump(out, open("runs/weighted/rm20_plane_macwilliams.json", "w"), indent=1)
print("dual weights 0..11:", B)
print("A_8: planar", planes * B[8], "frame", frame[8], "| A_9: planar", planes * B[9], "frame", frame[9], "| agree", out["agree"])
print("lines (B_w, 56 lines x words):", out["line_check"])
