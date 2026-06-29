# Method 2 (structured_enum) — runbook

Rigorous structured enumeration of low-weight RM_3(9,7) codewords for the m=7 cap
qutrit code `[[1968,219]]_3`, distributed over two machines.

`d = min_{c in RM_3(9,7), c!=0} |supp(c) \ S|`, S = the 219 cap points.
Method 2 lower-bounds `d` by minimising `|supp(c) \ S|` over the *structured* families
(unions of `j` parallel `w`-flats — genuine RM_3(9,7) codewords). For each family the
computation is the exact additive reduction "per direction V, top-j coset occupancy of
S", run as a numba `prange` hot loop over the `gaussian_binomial(7,w)` directions.

## What runs

| family | structure | weight | directions | full support count |
|--------|-----------|--------|-----------:|-------------------:|
| `w18`  | 2 parallel 2-flats (the complete min-weight class) | 18 | 99,463 | 2,924,510,589 |
| `w27`  | affine 3-flat (degree-8 indicator) | 27 | 925,771 | 74,987,451 |

Result (single machine, already verified): `w18 -> min 10`, `w27 -> min 18`, so the
rigorous min over enumerated families is **10** — consistent with `d >= 10`.

> The weight-18 class alone is exact and complete. Non-flat Leducq "second/third weight"
> classes (18 < wt <= 45 with non-flat support) are NOT enumerated here and are covered
> probabilistically by Method 1 (Stern). The m=4 cautionary case (a *higher*-weight
> codeword drops below the class minimum) is exactly why Method 1 is also required.

## Self-tests (run first, no network)

```bash
cd /home/hlamm/Desktop/QC/prime_msd
python -m cap_validate.structured.selftest        # build, structure, count, kernel==brute, minima
python -m cap_validate.structured.mock_cluster    # LOCAL two-node split == single-process (+resume)
python -m cap_validate.structured.run_all         # single-machine Method-2 verdict
```

## Two-machine distribution

The directions of each family are split into disjoint slices; each node runs
`max_topj` over its slice and writes `.partial` JSON files; the merge takes the MAX of
per-block `max|supp∩S|` -> MIN punctured weight. The split commutes with the merge, so
the distributed result is bit-identical to the single-machine result (asserted by
`mock_cluster`).

### Parameters (fill these in — no credentials are stored in the repo)

```bash
# ---- EDIT ME ----
REMOTE_USER=youruser
REMOTE_HOST=lenore
REMOTE_PORT=60022                 # default
REMOTE_KEY=$HOME/.ssh/id_ed25519  # path to your private key
REMOTE_DIR=/home/$REMOTE_USER/cap_validate_run
LOCAL_CORES=20
REMOTE_CORES=32
FAMILY=w18                        # then repeat for w27
REPO=/home/hlamm/Desktop/QC/prime_msd
SHARED=$REPO/cap_validate/structured/_run                # local collection dir
# -----------------
```

### 1. Build the plan (proportional to core counts) on the local box

```bash
cd $REPO
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
    args = " ".join(f"{b.start} {b.count}" for b in blocks)
    open(f"plan_node{ni}.txt","w").write(args)
    print("node", ni, "dirs", sum(b.count for b in blocks), "blocks", len(blocks))
PY
```

`plan_node0.txt` -> local node, `plan_node1.txt` -> remote (lenore).

### 2. Stage the remote node (one rsync; numba JIT compiles on first call there)

```bash
ssh -p $REMOTE_PORT -i $REMOTE_KEY $REMOTE_USER@$REMOTE_HOST "mkdir -p $REMOTE_DIR"
rsync -az -e "ssh -p $REMOTE_PORT -i $REMOTE_KEY" \
      $REPO/qmsd $REPO/cap_validate $REPO/cap_qutrit_code.json $REPO/pyproject.toml \
      $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/
```

### 3. Launch both nodes

Local:
```bash
mkdir -p $SHARED
cd $REPO
NUMBA_NUM_THREADS=$LOCAL_CORES python -m cap_validate.structured.distribute \
    local $FAMILY $(cat plan_node0.txt) $SHARED
```

Remote (lenore):
```bash
ssh -p $REMOTE_PORT -i $REMOTE_KEY $REMOTE_USER@$REMOTE_HOST \
  "cd $REMOTE_DIR && NUMBA_NUM_THREADS=$REMOTE_CORES python -m cap_validate.structured.distribute \
   lenore $FAMILY $(cat plan_node1.txt) $REMOTE_DIR/_run"
```

### 4. Collect remote partials and merge

```bash
rsync -az -e "ssh -p $REMOTE_PORT -i $REMOTE_KEY" \
      $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/_run/ $SHARED/
cd $REPO
python - <<PY
from cap_validate.structured import enum_families as ef, distribute as dist
fam="$FAMILY"; n=ef.FAMILIES[fam].n_directions
print(dist.merge_partials("$SHARED", fam, n))
PY
```

Repeat steps 1–4 with `FAMILY=w27`. The Method-2 verdict is `min` of the per-family
`min_punct` (and `run_all.py` reproduces it single-machine for cross-check).

### Resume

`.partial` files are idempotent and self-describing; re-running a node skips slices whose
`.partial` already exists and validates. After any crash, just re-run step 3 for that node
and re-merge — `mock_cluster` exercises exactly this path.
