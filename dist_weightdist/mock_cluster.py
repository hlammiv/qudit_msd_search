"""Local two-mock-node simulation of the 2-machine job (design sec.8).

HARD CONSTRAINT honored: this NEVER ssh's / connects to lenore_remote or any host.
The two-machine topology is simulated as two LOCAL node runs over a shared directory;
"transfer" is a filesystem copy.  Run from /home/hlamm/Desktop/QC/prime_msd:

    python -m dist_weightdist.mock_cluster            # default oracle [[206,37,4]]_3

Test plan (all asserted):
  1. SPLIT: partition.make_jobs splits the blocks across mock-node-A and -B
     (assert_cover: disjoint + exhaustive).
  2. INDEPENDENT ENUMERATION: each mock node runs dwd_core on its chunks -> .partial
     files; per-chunk + per-node checksums verified.
  3. TRANSFER + MERGE: collect both nodes' partials, merge_partials(expected=3^K)
     succeeds, merged A == in-process brute enumeration (vector equality).
  4. CHECKPOINT/RESUME: drop node-B's last partial + its ledger entry (simulating a
     mid-run crash); re-run node-B; assert the resumed merged A is IDENTICAL and no
     block is double-counted.
  5. END TO END: merged A -> full MacWilliams -> (d, A_d) matches qmsd on the oracle.
"""
from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import numpy as np

from qmsd import weightdist as qwd
from qmsd.oracle import load_oracle
from qmsd.triorthogonal import build_triorthogonal_code

from . import merge as mg
from . import partition as pt
from .checkpoint import Ledger
from .harness import core_binary, run_node
from .io_g0 import write_g0


def _build_g0(label):
    oc = next(o for o in load_oracle() if o.label == label)
    b = build_triorthogonal_code(oc.p, oc.m, oc.r_max, oc.puncture_columns_1indexed)
    return np.asarray(b["X_stab"]), oc.p


def run(label="[[206,37,4]]_3", exp_d=4, exp_Ad=880, t=2):
    core = core_binary()
    G0, q = _build_g0(label)
    assert q == 3
    K, n = G0.shape
    print(f"[mock] {label}  K={K} n={n}  t={t}  3^t={3**t} blocks  "
          f"3^K={3**K} messages")

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        g0path = td / "G0.g0"
        write_g0(g0path, G0, q)

        # --- 1. SPLIT across two mock nodes (unequal weights to exercise balance) ---
        plan = pt.make_jobs(K, t, ["nodeA", "nodeB"], [2.0, 1.0], chunk_blocks=max(1, 3 ** t // 4))
        all_chunks = [(c["start"], c["count"])
                      for jb in plan["nodes"].values() for c in jb["chunks"]]
        pt.assert_cover(all_chunks, 3 ** t)
        print(f"[ok] 1. SPLIT disjoint+exhaustive: "
              f"A={plan['nodes']['nodeA']['block_count']} blocks, "
              f"B={plan['nodes']['nodeB']['block_count']} blocks")

        # --- 2. INDEPENDENT ENUMERATION: each node runs its chunks ---
        outA, outB = td / "outA", td / "outB"
        ledA, ledB = td / "ledgerA.jsonl", td / "ledgerB.jsonl"
        sA = run_node(plan, "nodeA", str(g0path), str(outA), str(ledA), core=str(core))
        sB = run_node(plan, "nodeB", str(g0path), str(outB), str(ledB), core=str(core))
        mpb = 3 ** (K - t)
        assert sA["messages"] == plan["nodes"]["nodeA"]["block_count"] * mpb
        assert sB["messages"] == plan["nodes"]["nodeB"]["block_count"] * mpb
        print(f"[ok] 2. ENUMERATION: A ran {sA['chunks_run']} chunks "
              f"({sA['messages']} msgs), B ran {sB['chunks_run']} chunks "
              f"({sB['messages']} msgs)")

        # --- 3. TRANSFER + MERGE (filesystem copy) ---
        collected = td / "collected"
        collected.mkdir()
        for d in (outA, outB):
            for f in d.glob("*.partial"):
                shutil.copy(f, collected / f.name)
        partials = sorted(str(p) for p in collected.glob("*.partial"))
        A_merged, nn, KK = mg.merge_to_A(partials)        # expected_total=3^K enforced
        A_brute = _brute(G0, q)
        assert A_merged == A_brute, "merged A != in-process brute"
        print(f"[ok] 3. TRANSFER+MERGE: {len(partials)} partials merged, "
              f"sum=3^{K}={sum(A_merged)}, == brute enumeration")

        # --- 4. CHECKPOINT/RESUME: simulate node-B crash, resume ---
        bpart = sorted(outB.glob("*.partial"))
        assert bpart, "node B produced no partials"
        victim = bpart[-1]
        # drop the last partial file AND its ledger entry (as if crash before fsync+record)
        led = Ledger(ledB)
        vmissing_id = [cid for cid, rec in led.done.items()
                       if rec["partial_file"] == str(victim)][0]
        victim.unlink()
        # rewrite ledger without the victim entry (atomic-ish for the test)
        keep = [rec for cid, rec in led.done.items() if cid != vmissing_id]
        ledB.write_text("".join(__import__("json").dumps(r) + "\n" for r in keep))
        # resume node-B: only the missing chunk should be recomputed
        sB2 = run_node(plan, "nodeB", str(g0path), str(outB), str(ledB), core=str(core))
        assert sB2["chunks_run"] == 1, f"resume recomputed {sB2['chunks_run']} chunks (want 1)"
        # re-collect and re-merge -> identical A, no double count
        for f in outB.glob("*.partial"):
            shutil.copy(f, collected / f.name)
        partials2 = sorted(str(p) for p in collected.glob("*.partial"))
        A_resumed, _, _ = mg.merge_to_A(partials2)
        assert A_resumed == A_merged, "resume changed merged A (double count / loss!)"
        print(f"[ok] 4. CHECKPOINT/RESUME: dropped 1 partial, resume recomputed exactly "
              f"1 chunk, merged A IDENTICAL (no double-count)")

        # --- 5. END TO END certification ---
        res = mg.certify(partials2, search_kmax=exp_d + 2)
        ref = qwd.exact_distance_and_Ad(G0, q, max_words=q ** K + 1)
        assert res["d"] == ref["distance"] == exp_d, (res["d"], ref["distance"], exp_d)
        assert res["A_d"] == ref["B_d"] == exp_Ad, (res["A_d"], ref["B_d"], exp_Ad)
        print(f"[ok] 5. END-TO-END: d={res['d']} A_d={res['A_d']} "
              f"(== qmsd & published {exp_Ad}); dual invariants pass")

    print("ALL MOCK-CLUSTER TESTS PASSED (2 nodes simulated locally, no SSH)")


def _brute(G0, q):
    M = (np.asarray(G0).astype(np.int64)) % q
    K, n = M.shape
    total = q ** K
    hist = np.zeros(n + 1, dtype=np.int64)
    qpows = q ** np.arange(K, dtype=np.int64)
    chunk = max(1, 4_000_000 // (n + 1))
    start = 0
    while start < total:
        end = min(start + chunk, total)
        idx = np.arange(start, end, dtype=np.int64)
        msgs = (idx[:, None] // qpows[None, :]) % q
        cws = (msgs @ M) % q
        hist += np.bincount(np.count_nonzero(cws, axis=1), minlength=n + 1).astype(np.int64)
        start = end
    return [int(v) for v in hist]


if __name__ == "__main__":
    run()
