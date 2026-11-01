"""Build the validation instances for the ISD route and save them as .npz (matrix H = G0 or the
class block Gc, p, provenance).  Mirrors scratchpad/certify_candidate.py (p = 7 class block) and
scratchpad/rm_m4_scan.py (p = 11 G0).
usage: build_instances.py OUTDIR  [p7 draws + the p11 k=2290 candidate]
       build_instances.py OUTDIR CAND.json   [one rm_m4 candidate json -> OUTDIR/<stem>.npz]"""
import sys, os, json, glob, time, numpy as np
from qmsd.reedmuller import points, rm_generator
from qmsd.weighted import weighted_support, quadric, class_points
from qmsd.structured_ad import _rank
from qmsd.triorthogonal import shorten_matrix

outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)


def p11_from_candidate(cand_json, out):
    c = json.load(open(cand_json)); p, m, r = c["p"], c["m"], c["r"]
    cols0 = [x - 1 for x in c["puncture_columns_1indexed"]]
    G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
    G0 = np.ascontiguousarray(np.asarray(shorten_matrix(G, cols0, p), dtype=np.int64) % p)
    G0 = G0[[t for t in range(G0.shape[0]) if G0[t].any()]]
    rk = _rank(G0.copy(), p)
    np.savez_compressed(out, H=G0, p=p, rank=rk, source=os.path.abspath(cand_json), k=c["k"])
    print(f"{out}: G0 {G0.shape} rank {rk} (dim_A {c['dim_A']})", flush=True)


def p7_class_block(k, seed, files, out):
    p, m, r = 7, 4, 7
    rec = None
    for f in files:
        for l in open(f):
            if not l.strip(): continue
            rr = json.loads(l)
            if rr["k"] == k and rr["seed"] == seed and "cols_in_X" in rr: rec = rr
    assert rec is not None
    cols = rec["cols_in_X"]
    g = quadric(p, m, (1, 1, 1, 1)); GX, X, mu = weighted_support(p, m, g, r)
    cset = set(class_points(mu, p, 1))
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
    keep = [l for l in range(GX.shape[1]) if l not in set(cols)]
    n = G0.shape[1]
    is_class = np.array([keep[j] in cset for j in range(n)])
    C_cols = np.nonzero(is_class)[0]
    Gc = np.ascontiguousarray(G0[:, C_cols]); Gc = Gc[[t for t in range(Gc.shape[0]) if Gc[t].any()]]
    rk = _rank(Gc.copy(), p)
    np.savez_compressed(out, H=Gc, p=p, rank=rk, k=k, seed=seed, d_plane=rec["d_plane"], n_full=n,
                        class_cols=C_cols, G0_rows=G0.shape[0])
    print(f"{out}: n_full {n}, class block {Gc.shape} rank {rk}, d_plane {rec['d_plane']}", flush=True)


t0 = time.time()
if len(sys.argv) > 2:
    p11_from_candidate(sys.argv[2], os.path.join(outdir, os.path.basename(sys.argv[2]).replace(".json", ".npz")))
else:
    p7_class_block(272, 3, ["runs/weighted/orbit_sweep/orbit_shift_k272.jsonl"], os.path.join(outdir, "p7_k272_seed3_block.npz"))
    p7_class_block(270, 19, glob.glob("runs/weighted/plane_sweep/plane_sweep_k*.jsonl"), os.path.join(outdir, "p7_k270_seed19_block.npz"))
    p11_from_candidate("runs/rm_m4/p11_r13_k2290_seed1.json", os.path.join(outdir, "p11_r13_k2290_seed1.npz"))
print(f"done ({time.time()-t0:.0f}s)")
