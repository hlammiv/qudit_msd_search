# Cap qutrit `[[1968,219]]_3` — full deep distributed validation runbook

Two-machine deep pass for **local (20 cores)** + **lenore (`ssh -p 60022`, ~32 cores)**.

> **Heads-up before you run anything.** The verdict is already settled, *rigorously and
> exactly*, on the local box: `G0` has 83 zero columns => a weight-1 codeword of `G0^perp`
> exists => `d = 1 < 10` => `gamma >= 1`, **refuted** (see `CAP_VALIDATION.md`). The deep
> distributed pass below does **not** change this; it independently re-derives it and maps
> the full low-weight spectrum. Run it for confirmation/coverage, not to decide the verdict.

All paths absolute. Code lives in `/home/hlamm/Desktop/QC/prime_msd/cap_validate/`; run
Python from `/home/hlamm/Desktop/QC/prime_msd`. No credentials are stored in the repo — you
fill the `EDIT ME` block.

---

## 0. Self-tests first (no network — proves the tools before trusting them)

```bash
cd /home/hlamm/Desktop/QC/prime_msd
python cap_validate/test_stern_isd.py            # M1 vs qmsd exact_distance_and_Ad on known-d codes
python -m cap_validate.stern_distribute mock     # M1 two-mock-node seed split + idempotent resume
python -m cap_validate.structured.selftest       # M2 build/structure/count/kernel==brute/minima
python -m cap_validate.structured.mock_cluster   # M2 two-mock-node block split == single-process
```
All must pass before deploying.

## 0b. One-command local refutation (30 s, no cluster needed)

```bash
cd /home/hlamm/Desktop/QC/prime_msd
python cap_validate/run_cap_stern.py --minutes 5 --threads 20   # zero-column scan refutes; Stern confirms witness
python -m cap_validate.structured.run_all                       # M2 rigorous min over weight<=45 = 10 (necessary, NOT sufficient)
```
`run_cap_stern.py` reports `trivial_min_weight = 1` (the 83 zero columns) and a probabilistic
weight-2 witness within ~12 s. Either already proves `d < 10`.

---

## 1. Fill in your parameters (no secrets in the repo)

```bash
# ---- EDIT ME ----
REMOTE_USER=youruser
REMOTE_HOST=lenore
REMOTE_PORT=60022                  # default
REMOTE_KEY=$HOME/.ssh/id_ed25519   # path to your PRIVATE key
REMOTE_DIR=/home/$REMOTE_USER/cap_validate_run
LOCAL_CORES=20
REMOTE_CORES=32
REPO=/home/hlamm/Desktop/QC/prime_msd
SHARED=$REPO/cap_validate/_run     # local collection dir
SHARED_M2=$REPO/cap_validate/structured/_run
# -----------------
mkdir -p "$SHARED" "$SHARED_M2"
```

Sanity-check connectivity (optional):
```bash
ssh -p $REMOTE_PORT -i $REMOTE_KEY $REMOTE_USER@$REMOTE_HOST 'python -c "import numpy, numba; print(\"ok\", numba.__version__)"'
```

---

## 2. Deploy `cap_validate/` to lenore (one rsync)

```bash
ssh -p $REMOTE_PORT -i $REMOTE_KEY $REMOTE_USER@$REMOTE_HOST "mkdir -p $REMOTE_DIR/_run $REMOTE_DIR/structured/_run"
rsync -az -e "ssh -p $REMOTE_PORT -i $REMOTE_KEY" \
      $REPO/qmsd $REPO/cap_validate $REPO/cap_qutrit_code.json $REPO/pyproject.toml \
      $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/
```
numba JIT compiles on first call on lenore (one-time ~30-60 s).

---

## 3. Method 1 (Stern / ISD) — split disjoint random-seed ranges

Build a core-proportional plan (disjoint seed ranges; merge = global MIN weight):

```bash
cd $REPO
TOTAL=200000                       # total ISD iterations (deep). Scale up freely.
python - <<PY
from cap_validate.stern_distribute import split_seed_ranges
for name,(s,it) in zip(["local","lenore"],
        split_seed_ranges($TOTAL,[float($LOCAL_CORES),float($REMOTE_CORES)])):
    open(f"stern_plan_{name}.txt","w").write(f"{s} {it}")
    print(name, "seed_start", s, "iters", it)
PY
```

Local node:
```bash
read S IT < stern_plan_local.txt
NUMBA_NUM_THREADS=$LOCAL_CORES python -m cap_validate.stern_distribute local \
  --H cap --seed-start $S --iters $IT --out $SHARED --threads $LOCAL_CORES
```

lenore node (run concurrently):
```bash
read S IT < stern_plan_lenore.txt
ssh -p $REMOTE_PORT -i $REMOTE_KEY $REMOTE_USER@$REMOTE_HOST \
  "cd $REMOTE_DIR && NUMBA_NUM_THREADS=$REMOTE_CORES python -m cap_validate.stern_distribute \
   lenore --H cap --seed-start $S --iters $IT --out $REMOTE_DIR/_run --threads $REMOTE_CORES"
```

Collect lenore's partials and merge (global minimum weight + witness):
```bash
rsync -az -e "ssh -p $REMOTE_PORT -i $REMOTE_KEY" \
      $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/_run/ $SHARED/
python - <<PY
from cap_validate.stern_distribute import merge_partials
print(merge_partials("$SHARED"))
PY
```

**Read the M1 verdict.** The merge prints the global lowest weight + witness support. Verify any
witness directly: `(G0 @ witness) % 3 == 0`. Lowest weight `< 10` => `d < 10` => **refuted**.
(For this code it will report `1` from the exact zero-column scan and `2` from probabilistic Stern.)

To go deeper: raise `TOTAL`, and sweep higher `p` with a randomized `info_width` (the harness
keeps `l < r`; do not set `l >= r=55` — that degenerate window silently under-reports).

---

## 4. Method 2 (structured) — split disjoint direction blocks

Per family (`w18` then `w27`); merge = MAX of per-block `max|supp∩S|` -> MIN punctured weight.
The split commutes with the merge (asserted by `mock_cluster`), so distributed == single-machine.

```bash
cd $REPO
export LOCAL_CORES REMOTE_CORES
for FAMILY in w18 w27; do
  export FAMILY
  python - <<'PY'
import os
from cap_validate.structured import enum_families as ef, distribute as dist
fam = os.environ["FAMILY"]
n = ef.FAMILIES[fam].n_directions
plan = dist.split_blocks(fam, n,
        node_weights=[float(os.environ["LOCAL_CORES"]), float(os.environ["REMOTE_CORES"])],
        chunk=12000)
dist.assert_cover(plan, n)
for ni, blocks in plan.items():
    open(f"m2_plan_{fam}_node{ni}.txt","w").write(
        " ".join(f"{b.start} {b.count}" for b in blocks))
    print(fam, "node", ni, "dirs", sum(b.count for b in blocks), "blocks", len(blocks))
PY
done
```

Run each family on both nodes (local example for `w18`; repeat with `w27`, and on lenore over SSH
pointing at `$REMOTE_DIR/structured/_run`):
```bash
# local
python -m cap_validate.structured.distribute run --family w18 \
  --plan $REPO/m2_plan_w18_node0.txt --out $SHARED_M2 --threads $LOCAL_CORES
# lenore
ssh -p $REMOTE_PORT -i $REMOTE_KEY $REMOTE_USER@$REMOTE_HOST \
  "cd $REMOTE_DIR && python -m cap_validate.structured.distribute run --family w18 \
   --plan $REMOTE_DIR/m2_plan_w18_node1.txt --out $REMOTE_DIR/structured/_run --threads $REMOTE_CORES"
```
(If your `distribute` CLI differs, consult `cap_validate/structured/RUNBOOK.md`, the authoritative
M2 command reference — it carries the same host/port/key parameters.)

Collect + merge:
```bash
rsync -az -e "ssh -p $REMOTE_PORT -i $REMOTE_KEY" \
      $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/structured/_run/ $SHARED_M2/
python -m cap_validate.structured.distribute merge --out $SHARED_M2
```

**Read the M2 verdict.** It prints the rigorous min `|supp\S|` over the enumerated families
(expected: `w18 -> 10`, `w27 -> 18`). This is a lower bound on `d` **only within those families** —
it is NECESSARY but NOT SUFFICIENT for `d >= 10` (it cannot see RM-weight `> 45`, where the actual
refuting codewords live).

---

## 5. Resume / fault tolerance

Both methods write idempotent, self-describing `.partial` JSON. Re-running a node skips completed
slices/ranges — just re-run that node's command and re-merge. The two `mock` self-tests exercise the
wipe-one-node-and-rerun path.

---

## 6. Reading the combined verdict

| M1 lowest weight | M2 structured min | meaning for `gamma < 1` |
|---|---|---|
| `< 10` (witness verified) | (irrelevant) | **REFUTED.** `d < 10`, `gamma >= 1`. Proof = the witness. **<- this code: M1 = 1.** |
| `>= 10` after a deep sweep | `>= 10` | **Supported-so-far.** Rigorous for RM-weight<=45 (M2) + strong probabilistic evidence beyond (M1). High-confidence, **not a proof** of `d >= 10`. |

For `[[1968,219]]_3` the top row applies and is *exact* (not probabilistic): the 83 zero columns of
`G0` are weight-1 codewords of `G0^perp`. `d = 1`. `gamma >= 1`. The `gamma < 1` claim is refuted.

A genuine *proof* of `d >= 10` (not relevant here) would need either a complete rigorous structured
enumeration across all RM weights or an exact global method — both infeasible at `3^55` / `C(1968,5)`.
