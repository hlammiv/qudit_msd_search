# Cap qutrit code `[[1968,219]]_3` — two-machine validation runbook

Validates whether the m=7 cap qutrit code has minimum distance `d >= 10`
(=> `gamma = log(1968/219)/log(10) = 0.954 < 1`) or `d < 10` (=> refuted).

`d = min_{c in RM_3(9,7), c != 0} |supp(c) \ S|`, with `S` = the 219 cap puncture
points; `G0 = X_stab` (55 x 1968 over F_3) and `G0^perp = ker(G0) =` punctured RM_3(9,7).

Two complementary validators, each distributable across **local (20 cores)** and
**lenore (SSH port 60022, ~32 cores)**:

| method | what | guarantee |
|---|---|---|
| **M1 Stern/ISD** (`cap_validate/stern_isd.py`, `stern_distribute.py`) | probabilistic low-weight-codeword search on `G0^perp`; **plus an O(n) exact zero-/proportional-column scan** (`trivial_low_weight`) | finding `w<10` PROVES `d<10` (witness); finding nothing is probabilistic evidence |
| **M2 structured enum** (`cap_validate/structured/`) | rigorous min `|supp\S|` over RM_3(9,7) codewords of weight 18..45 | a true LOWER BOUND on `d` *restricted to RM-weight<=45* |

> **RESULT (already established locally, rigorously): `d = 1`, code REFUTED.**
> `G0` has **83 zero columns**; each is a genuine weight-1 codeword of punctured
> RM_3(9,7) (verified: `e_j in rowspan(G9_punctured)`, `rank=1913`). They are the
> punctures of RM_3(9,7) codewords of full weight **77, 103, 112, 145, ...** whose
> support is almost entirely inside the cap `S` (e.g. 76 of 77 points in `S`). These
> high-weight codewords are exactly what M2's weight<=45 window and the prior
> "cap meets any 2-flat in <=4 points" argument (which only bounds the weight-18 class)
> never covered, and what the infeasible 3^55 MITM could not reach. `gamma >= 1`.

## Self-tests (run first, no network)

```bash
cd /home/hlamm/Desktop/QC/prime_msd
python cap_validate/test_stern_isd.py               # M1 vs exact engine on known-d codes
python -m cap_validate.stern_distribute mock        # M1 two-mock-node split + resume
python -m cap_validate.structured.selftest          # M2 build/structure/count/kernel/minima
python -m cap_validate.structured.mock_cluster      # M2 two-mock-node split + resume
```

## Single-machine runs (reproduce the verdict)

```bash
cd /home/hlamm/Desktop/QC/prime_msd
# M1: the zero-column scan alone refutes; Stern confirms with a witness
python cap_validate/run_cap_stern.py --minutes 5 --threads 20
# M2: rigorous min over weight<=45 structured families (=> 10, NECESSARY but NOT sufficient)
python -m cap_validate.structured.run_all
```

## Parameterized two-machine deployment (fill these in — no credentials in repo)

```bash
# ---- EDIT ME ----
REMOTE_USER=youruser
REMOTE_HOST=lenore
REMOTE_PORT=60022                  # default
REMOTE_KEY=$HOME/.ssh/id_ed25519   # private key path
REMOTE_DIR=/home/$REMOTE_USER/cap_validate_run
LOCAL_CORES=20
REMOTE_CORES=32
REPO=/home/hlamm/Desktop/QC/prime_msd
SHARED=$REPO/cap_validate/_run     # local collection dir
# -----------------
mkdir -p "$SHARED"
```

### Stage lenore (one rsync; numba JIT compiles on first call there)

```bash
ssh -p $REMOTE_PORT -i $REMOTE_KEY $REMOTE_USER@$REMOTE_HOST "mkdir -p $REMOTE_DIR"
rsync -az -e "ssh -p $REMOTE_PORT -i $REMOTE_KEY" \
      $REPO/qmsd $REPO/cap_validate $REPO/cap_qutrit_code.json $REPO/pyproject.toml \
      $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/
```

### Method 1 (Stern) — split random-seed ranges across nodes

Plan (proportional to cores) — disjoint seed ranges, merge = global MIN weight:

```bash
cd $REPO
TOTAL=200000
python - <<PY
from cap_validate.stern_distribute import split_seed_ranges
for name,(s,it) in zip(["local","lenore"],
        split_seed_ranges($TOTAL,[float($LOCAL_CORES),float($REMOTE_CORES)])):
    open(f"stern_plan_{name}.txt","w").write(f"{s} {it}")
    print(name, s, it)
PY
```

Local node:
```bash
read S IT < stern_plan_local.txt
NUMBA_NUM_THREADS=$LOCAL_CORES python -m cap_validate.stern_distribute local \
  --H cap --seed-start $S --iters $IT --out $SHARED --threads $LOCAL_CORES
```

lenore:
```bash
read S IT < stern_plan_lenore.txt
ssh -p $REMOTE_PORT -i $REMOTE_KEY $REMOTE_USER@$REMOTE_HOST \
  "cd $REMOTE_DIR && NUMBA_NUM_THREADS=$REMOTE_CORES python -m cap_validate.stern_distribute \
   lenore --H cap --seed-start $S --iters $IT --out $REMOTE_DIR/_run --threads $REMOTE_CORES"
```

Collect + merge:
```bash
rsync -az -e "ssh -p $REMOTE_PORT -i $REMOTE_KEY" \
      $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR/_run/ $SHARED/
python - <<PY
from cap_validate.stern_distribute import merge_partials
print(merge_partials("$SHARED"))
PY
```

### Method 2 (structured) — split direction blocks across nodes

See `cap_validate/structured/RUNBOOK.md` (same host/port/key parameters; family `w18`
then `w27`; merge = MAX of block `max|supp∩S|` -> MIN punctured weight).

### Resume

Both methods write idempotent, self-describing `.partial` JSON. Re-running a node skips
already-completed slices/ranges; just re-run that node's command and re-merge. The two
`mock` self-tests exercise the wipe-one-node-and-rerun path.
```
```
