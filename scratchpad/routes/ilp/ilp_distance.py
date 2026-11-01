"""Route: integer / constraint programming / SAT for the minimum-distance certificate.

Decide whether ker(G) (G an r x n matrix over F_p, a parity-check matrix) contains a nonzero
vector of weight exactly w (``--exact``, assumes weights < w were excluded, like the MITM
scans) or of weight in [1, w].  Three encodings:

  cpsat-onehot   booleans x[j][v] = [e_j = v], v = 1..p-1 (at most one per column); for every
                 row t of the row-reduced matrix sum_j,v ((G[t,j] v) mod p) x[j][v] = p q_t with
                 q_t an integer in [0, floor(w (p-1)/p)];  sum_j s_j = w (or 1 <= sum <= w).
  cpsat-sys      information-set form: row-reduce G to [I | A] (pivot set I, |I| = rank), free
                 part e_J one-hot as above, y_t = (A e_J) mod p via one-hot residues per pivot
                 row, z_t = [y_t != 0];  wt(e_J) + sum_t z_t = w.  Fewer free columns, no q.
  sat            CNF for CaDiCaL (pysat): one-hot values, one residue automaton per row
                 (states Z_p, p*(p) transition clauses per column), sequential-counter
                 cardinality on the support indicators.

Symmetry breaking (all encodings): scalar multiples -- the FIRST nonzero coordinate equals 1
(prefix booleans b_j = OR_{i<=j} s_i; x[j][v] for v >= 2 implies b_{j-1}).  Support ordering is
implicit (the support is a set of column indicators, no permutation symmetry).  Column-group
symmetry: ``--orbit-reps`` takes a list of column orbit representatives of a known automorphism
group and adds "support meets the representatives"; none of the test instances has a known
nontrivial group (random puncture sets), so it is unused in the measurements.

usage:
  ilp_distance.py --instance {A,B,C} [--subset N --seed S] --w W [--exact] \
                  [--encoding cpsat-onehot|cpsat-sys|sat|mitm] [--time SEC]
  ilp_distance.py --matrix FILE.npy --p P --w W ...
Prints a JSON line with status (INFEASIBLE / FEASIBLE with the word / TIMEOUT), timing and solver
statistics.  Always run with one thread (OMP_NUM_THREADS=1; CP-SAT num_workers = 1).
"""
from __future__ import annotations
import sys, os, json, time, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))


# ----------------------------------------------------------------------------- F_p linear algebra
def rref_mod_p(G, p):
    """Row-reduced echelon form over F_p; returns (R, pivots) with R of shape (rank, n)."""
    A = np.array(G, dtype=np.int64) % p
    r, n = A.shape
    piv = []
    row = 0
    for col in range(n):
        if row >= r:
            break
        nz = np.nonzero(A[row:, col])[0]
        if len(nz) == 0:
            continue
        k = row + nz[0]
        if k != row:
            A[[row, k]] = A[[k, row]]
        inv = pow(int(A[row, col]), -1, p)
        A[row] = (A[row] * inv) % p
        others = np.nonzero(A[:, col])[0]
        for t in others:
            if t != row:
                A[t] = (A[t] - A[t, col] * A[row]) % p
        piv.append(col)
        row += 1
    return A[:row], piv


def verify_word(G, p, cols, vals):
    e = np.zeros(G.shape[1], dtype=np.int64)
    for c, v in zip(cols, vals):
        e[c] = v
    return bool(((np.asarray(G, dtype=np.int64) @ e) % p == 0).all())


# ----------------------------------------------------------------------------- instances
def _load_cached(name, builder):
    os.makedirs(CACHE, exist_ok=True)
    f = os.path.join(CACHE, name + ".npy")
    if os.path.exists(f):
        return np.load(f)
    M = builder()
    np.save(f, M)
    return M


def _class_block(k, seed, draw_file):
    """Class block Gc of a p = 7 sphere-charge draw, as in scratchpad/certify_candidate.py."""
    sys.path.insert(0, ROOT)
    from qmsd.reedmuller import points
    from qmsd.weighted import weighted_support, quadric, class_points
    from qmsd.triorthogonal import shorten_matrix
    p, m, r = 7, 4, 7
    rec = None
    for l in open(os.path.join(ROOT, draw_file)):
        if not l.strip():
            continue
        rr = json.loads(l)
        if rr["k"] == k and rr["seed"] == seed and "cols_in_X" in rr:
            rec = rr
    assert rec is not None, "draw not found"
    cols = rec["cols_in_X"]
    g = quadric(p, m, (1, 1, 1, 1))
    GX, X, mu = weighted_support(p, m, g, r)
    cand = class_points(mu, p, 1); cset = set(cand)
    G0 = np.asarray(shorten_matrix(GX, cols, p), dtype=np.int64) % p
    keep = [l for l in range(GX.shape[1]) if l not in set(cols)]
    n = G0.shape[1]
    is_class = np.array([keep[j] in cset for j in range(n)])
    C_cols = np.nonzero(is_class)[0]
    Gc = np.ascontiguousarray(G0[:, C_cols])
    Gc = Gc[[t for t in range(Gc.shape[0]) if Gc[t].any()]]
    return Gc


def _p11_full():
    sys.path.insert(0, ROOT)
    from qmsd.reedmuller import rm_generator
    from qmsd.triorthogonal import shorten_matrix
    c = json.load(open(os.path.join(ROOT, "runs/rm_m4/p11_r13_k2290_seed1.json")))
    p, m, r = c["p"], c["m"], c["r"]
    cols0 = [x - 1 for x in c["puncture_columns_1indexed"]]
    G = np.asarray(rm_generator(r, m, p), dtype=np.int64) % p
    G0 = np.ascontiguousarray(np.asarray(shorten_matrix(G, cols0, p), dtype=np.int64) % p)
    G0 = G0[[t for t in range(G0.shape[0]) if G0[t].any()]]
    return G0


INSTANCES = {
    "A": (7, lambda: _class_block(270, 19, "runs/weighted/plane_sweep/plane_sweep_k270.jsonl")),
    "B": (7, lambda: _class_block(272, 3, "runs/weighted/orbit_sweep/orbit_shift_k272.jsonl")),
    "C": (11, _p11_full),
}


def load_instance(name, subset=None, seed=0):
    p, builder = INSTANCES[name]
    G = _load_cached(name, builder)
    cols = None
    if subset is not None and subset < G.shape[1]:
        rng = np.random.default_rng(seed)
        cols = np.sort(rng.choice(G.shape[1], subset, replace=False))
        G = np.ascontiguousarray(G[:, cols])
        G = G[[t for t in range(G.shape[0]) if G[t].any()]]
    return G, p, cols


# ----------------------------------------------------------------------------- CP-SAT encodings
def _cpsat_params(model_solver, time_limit, workers=1, seed=0):
    model_solver.parameters.max_time_in_seconds = float(time_limit)
    model_solver.parameters.num_workers = int(workers)
    model_solver.parameters.random_seed = int(seed)
    model_solver.parameters.log_search_progress = False


def _first_nonzero_is_one(m, x, s, p):
    """CP-SAT: b_j = OR_{i<=j} s_i; a value >= 2 in column j requires an earlier nonzero."""
    b_prev = None
    for j in range(len(s)):
        if b_prev is None:
            for v in range(2, p):
                m.Add(x[j][v - 1] == 0)
            b = s[j]
        else:
            b = m.NewBoolVar(f"b{j}")
            m.AddBoolOr([b_prev, s[j]]).OnlyEnforceIf(b)
            m.AddImplication(b_prev, b); m.AddImplication(s[j], b)
            for v in range(2, p):
                m.AddImplication(x[j][v - 1], b_prev)
        b_prev = b


def solve_cpsat_onehot(G, p, w, exact=True, time_limit=600, orbit_reps=None, workers=1,
                       linearization=None, symmetry=True):
    from ortools.sat.python import cp_model
    R, piv = rref_mod_p(G, p)
    rho, n = R.shape
    m = cp_model.CpModel()
    x = [[m.NewBoolVar(f"x{j}_{v}") for v in range(1, p)] for j in range(n)]
    s = [m.NewBoolVar(f"s{j}") for j in range(n)]
    for j in range(n):
        m.Add(sum(x[j]) == s[j])
    if exact:
        m.Add(sum(s) == w)
    else:
        m.Add(sum(s) <= w); m.Add(sum(s) >= 1)
    qmax = (w * (p - 1)) // p
    for t in range(rho):
        q = m.NewIntVar(0, qmax, f"q{t}")
        terms = []
        for j in range(n):
            c = int(R[t, j])
            if c == 0:
                continue
            for v in range(1, p):
                res = (c * v) % p
                if res:
                    terms.append(res * x[j][v - 1])
        m.Add(sum(terms) == p * q)
    if symmetry:
        _first_nonzero_is_one(m, x, s, p)
    if orbit_reps:
        m.AddBoolOr([s[j] for j in orbit_reps])
    solver = cp_model.CpSolver()
    _cpsat_params(solver, time_limit, workers)
    if linearization is not None:
        solver.parameters.linearization_level = linearization
    t0 = time.time()
    st = solver.Solve(m)
    dt = time.time() - t0
    out = {"solver": "cpsat-onehot", "rows": rho, "n": n, "w": w, "exact": exact, "seconds": round(dt, 2),
           "status": solver.StatusName(st), "conflicts": solver.NumConflicts(), "branches": solver.NumBranches()}
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        cols, vals = [], []
        for j in range(n):
            for v in range(1, p):
                if solver.Value(x[j][v - 1]):
                    cols.append(j); vals.append(v)
        out["word_cols"] = cols; out["word_vals"] = vals; out["verified"] = verify_word(G, p, cols, vals)
    return out


def solve_cpsat_sys(G, p, w, exact=True, time_limit=600, orbit_reps=None, workers=1, symmetry=True):
    """Information-set form.  e = (e_I, e_J) with pivot set I, e_I = -A e_J.  Free one-hot e_J;
    y_t = residue of (A e_J) mod p (sum = p q_t + y_t, q_t in [0, qmax]); z_t = [y_t != 0].
    Weight = wt(e_J) + sum z_t."""
    from ortools.sat.python import cp_model
    R, piv = rref_mod_p(G, p)
    rho, n = R.shape
    pivset = set(piv)
    J = [j for j in range(n) if j not in pivset]
    A = R[:, J]                              # R e = 0  <=>  e_I = -A e_J
    nJ = len(J)
    m = cp_model.CpModel()
    x = [[m.NewBoolVar(f"x{j}_{v}") for v in range(1, p)] for j in range(nJ)]
    s = [m.NewBoolVar(f"s{j}") for j in range(nJ)]
    for j in range(nJ):
        m.Add(sum(x[j]) == s[j])
    z = [m.NewBoolVar(f"z{t}") for t in range(rho)]
    qmax = (w * (p - 1)) // p
    for t in range(rho):
        y = m.NewIntVar(0, p - 1, f"y{t}")
        q = m.NewIntVar(0, qmax, f"q{t}")
        terms = []
        for j in range(nJ):
            c = int(A[t, j])
            if c == 0:
                continue
            for v in range(1, p):
                res = (c * v) % p
                if res:
                    terms.append(res * x[j][v - 1])
        m.Add(sum(terms) == p * q + y)
        m.Add(y <= (p - 1) * z[t]); m.Add(y >= z[t])
    if exact:
        m.Add(sum(s) + sum(z) == w)
    else:
        m.Add(sum(s) + sum(z) <= w)
    m.Add(sum(s) >= 1)                        # e_J = 0 forces e = 0
    if symmetry:
        _first_nonzero_is_one(m, x, s, p)     # scaling fixes the first nonzero FREE coordinate
    if orbit_reps:
        pos = {c: i for i, c in enumerate(piv)}
        posJ = {c: i for i, c in enumerate(J)}
        m.AddBoolOr([z[pos[c]] if c in pos else s[posJ[c]] for c in orbit_reps])
    solver = cp_model.CpSolver()
    _cpsat_params(solver, time_limit, workers)
    t0 = time.time()
    st = solver.Solve(m)
    dt = time.time() - t0
    out = {"solver": "cpsat-sys", "rows": rho, "n": n, "free": nJ, "w": w, "exact": exact, "seconds": round(dt, 2),
           "status": solver.StatusName(st), "conflicts": solver.NumConflicts(), "branches": solver.NumBranches()}
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        e = np.zeros(n, dtype=np.int64)
        for j in range(nJ):
            for v in range(1, p):
                if solver.Value(x[j][v - 1]):
                    e[J[j]] = v
        eI = (-(A @ e[J])) % p
        for t, c in enumerate(piv):
            e[c] = eI[t]
        cols = [int(c) for c in np.nonzero(e)[0]]; vals = [int(e[c]) for c in cols]
        out["word_cols"] = cols; out["word_vals"] = vals; out["verified"] = verify_word(G, p, cols, vals)
    return out


# ----------------------------------------------------------------------------- CNF / CaDiCaL
def solve_sat(G, p, w, exact=True, time_limit=600, orbit_reps=None, solver_name="cadical195", symmetry=True):
    from pysat.formula import CNF, IDPool
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Solver
    import threading
    R, piv = rref_mod_p(G, p)
    rho, n = R.shape
    pool = IDPool()
    X = lambda j, v: pool.id(("x", j, v))
    S = lambda j: pool.id(("s", j))
    cnf = CNF()
    for j in range(n):
        vs = [X(j, v) for v in range(1, p)]
        cnf.append([-S(j)] + vs)
        for v in vs:
            cnf.append([-v, S(j)])
        for a in range(len(vs)):                      # at most one value
            for b in range(a + 1, len(vs)):
                cnf.append([-vs[a], -vs[b]])
    # residue automaton per row: state st[t][j][r] = residue after columns 0..j-1
    for t in range(rho):
        ST = lambda j, r: pool.id(("st", t, j, r))
        cnf.append([ST(0, 0)])
        for r in range(1, p):
            cnf.append([-ST(0, r)])
        for j in range(n):
            c = int(R[t, j])
            for r in range(p):
                if c == 0:
                    cnf.append([-ST(j, r), ST(j + 1, r)])
                    continue
                cnf.append([-ST(j, r), S(j), ST(j + 1, r)])
                for v in range(1, p):
                    cnf.append([-ST(j, r), -X(j, v), ST(j + 1, (r + c * v) % p)])
            sts = [ST(j + 1, r) for r in range(p)]   # at most one state: functional chain
            for a in range(p):
                for b in range(a + 1, p):
                    cnf.append([-sts[a], -sts[b]])
        cnf.append([ST(n, 0)])
    lits = [S(j) for j in range(n)]
    if exact:
        card = CardEnc.equals(lits=lits, bound=w, vpool=pool, encoding=EncType.seqcounter)
    else:
        card = CardEnc.atmost(lits=lits, bound=w, vpool=pool, encoding=EncType.seqcounter)
        cnf.append(lits)
    cnf.extend(card.clauses)
    if symmetry:
        B = lambda j: pool.id(("b", j))
        for j in range(n):
            if j == 0:
                for v in range(2, p):
                    cnf.append([-X(j, v)])
                cnf.append([-B(0), S(0)]); cnf.append([B(0), -S(0)])
            else:
                cnf.append([-B(j), B(j - 1), S(j)]); cnf.append([B(j), -B(j - 1)]); cnf.append([B(j), -S(j)])
                for v in range(2, p):
                    cnf.append([-X(j, v), B(j - 1)])
    if orbit_reps:
        cnf.append([S(j) for j in orbit_reps])
    nv, nc = pool.top, len(cnf.clauses)
    t0 = time.time()
    with Solver(name=solver_name, bootstrap_with=cnf.clauses) as sol:
        # time limit by chunked conflict budgets (pysat's interrupt() does not stop CaDiCaL reliably)
        res = None
        while True:
            sol.conf_budget(5000)
            res = sol.solve_limited()
            if res is not None or time.time() - t0 > time_limit:
                break
        dt = time.time() - t0
        out = {"solver": "sat-" + solver_name, "rows": rho, "n": n, "w": w, "exact": exact, "vars": nv,
               "clauses": nc, "seconds": round(dt, 2),
               "status": {True: "FEASIBLE", False: "INFEASIBLE", None: "TIMEOUT"}[res]}
        try:
            out["conflicts"] = sol.accum_stats().get("conflicts")
        except Exception:
            pass
        if res:
            model = set(l for l in sol.get_model() if l > 0)
            cols, vals = [], []
            for j in range(n):
                for v in range(1, p):
                    if X(j, v) in model:
                        cols.append(j); vals.append(v)
            out["word_cols"] = cols; out["word_vals"] = vals; out["verified"] = verify_word(G, p, cols, vals)
    return out


# ----------------------------------------------------------------------------- MITM ground truth
def solve_mitm(G, p, w, exact=True):
    sys.path.insert(0, ROOT)
    from qmsd.mindist_nb import weight_d_exists_fast, weight6_witness_3x3
    G = np.ascontiguousarray(np.asarray(G, dtype=np.int64) % p)
    t0 = time.time()
    ds = [w] if exact else list(range(2, w + 1))
    found = None
    for d in ds:
        if d == 6:
            wit = weight6_witness_3x3(G, p)
            if wit is not None:
                found = d; break
        else:
            if weight_d_exists_fast(G, p, d):
                found = d; break
    dt = time.time() - t0
    return {"solver": "mitm", "n": G.shape[1], "w": w, "exact": exact, "seconds": round(dt, 2),
            "status": "FEASIBLE" if found else "INFEASIBLE", "found_weight": found}


# ----------------------------------------------------------------------------- CP-SAT ordered-support model
def solve_cpsat_elem(G, p, w, exact=True, time_limit=600, orbit_reps=None, workers=1, symmetry=True,
                     fixed_search=False):
    """The natural CP model: w ordered support positions.  For position i: column c_i (IntVar in
    [0, n)), value v_i in [1, p-1], flat index idx_i = (p-1) c_i + v_i - 1, and for each row t a
    residue r_{t,i} = (G[t, c_i] v_i) mod p obtained with an element constraint over a constant
    table of size n (p-1).  Row t: sum_i r_{t,i} = p q_t, q_t in [0, floor(w (p-1)/p)].
    Symmetry: c_1 < c_2 < ... < c_w (lexicographic support, kills position permutations) and
    v_1 = 1 (scalar multiples).  Pruning happens at depth w-1: the residues of the last
    position are forced by every row and the element domains intersect to empty.  ``exact``
    means exactly w positions; the non-exact variant loops d = 1..w."""
    from ortools.sat.python import cp_model
    R, piv = rref_mod_p(G, p)
    rho, n = R.shape
    ds = [w] if exact else list(range(1, w + 1))
    t_all = time.time()
    total_conf = total_br = 0
    for d in ds:
        m = cp_model.CpModel()
        c = [m.NewIntVar(0, n - 1, f"c{i}") for i in range(d)]
        v = [m.NewIntVar(1, p - 1, f"v{i}") for i in range(d)]
        idx = [m.NewIntVar(0, n * (p - 1) - 1, f"i{i}") for i in range(d)]
        for i in range(d):
            m.Add(idx[i] == (p - 1) * c[i] + v[i] - 1)
            if i:
                m.Add(c[i] >= c[i - 1] + 1)
        if symmetry:
            m.Add(v[0] == 1)
        qmax = (d * (p - 1)) // p
        for t in range(rho):
            table = [int((int(R[t, j]) * vv) % p) for j in range(n) for vv in range(1, p)]
            rs = []
            for i in range(d):
                r_ = m.NewIntVar(0, p - 1, f"r{t}_{i}")
                m.AddElement(idx[i], table, r_)
                rs.append(r_)
            q = m.NewIntVar(0, qmax, f"q{t}")
            m.Add(sum(rs) == p * q)
        if orbit_reps:
            lits = []
            for i in range(d):
                for j in orbit_reps:
                    b = m.NewBoolVar(f"o{i}_{j}")
                    m.Add(c[i] == j).OnlyEnforceIf(b)
                    lits.append(b)
            m.AddBoolOr(lits)
        solver = cp_model.CpSolver()
        _cpsat_params(solver, max(1.0, time_limit - (time.time() - t_all)), workers)
        if fixed_search:
            # plain DFS in support order: c_1, v_1, c_2, v_2, ... (the enumeration whose leaves sit at depth w-1)
            order = []
            for i in range(d):
                order += [c[i], v[i]]
            m.AddDecisionStrategy(order, cp_model.CHOOSE_FIRST, cp_model.SELECT_MIN_VALUE)
            solver.parameters.search_branching = cp_model.FIXED_SEARCH
        st = solver.Solve(m)
        total_conf += solver.NumConflicts(); total_br += solver.NumBranches()
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            cols = [solver.Value(ci) for ci in c]; vals = [solver.Value(vi) for vi in v]
            return {"solver": "cpsat-elem", "rows": rho, "n": n, "w": w, "exact": exact,
                    "seconds": round(time.time() - t_all, 2), "status": "FEASIBLE", "found_weight": d,
                    "conflicts": total_conf, "branches": total_br, "word_cols": cols, "word_vals": vals,
                    "verified": verify_word(G, p, cols, vals)}
        if st != cp_model.INFEASIBLE:
            return {"solver": "cpsat-elem", "rows": rho, "n": n, "w": w, "exact": exact,
                    "seconds": round(time.time() - t_all, 2), "status": solver.StatusName(st), "stuck_at_weight": d,
                    "conflicts": total_conf, "branches": total_br}
    return {"solver": "cpsat-elem", "rows": rho, "n": n, "w": w, "exact": exact,
            "seconds": round(time.time() - t_all, 2), "status": "INFEASIBLE", "conflicts": total_conf, "branches": total_br}


SOLVERS = {"cpsat-onehot": solve_cpsat_onehot, "cpsat-sys": solve_cpsat_sys, "sat": solve_sat, "mitm": solve_mitm,
           "cpsat-elem": solve_cpsat_elem,
           "cpsat-elem-dfs": lambda *a, **k: dict(solve_cpsat_elem(*a, fixed_search=True, **k), solver="cpsat-elem-dfs")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instance", choices=list(INSTANCES))
    ap.add_argument("--matrix"); ap.add_argument("--p", type=int)
    ap.add_argument("--subset", type=int); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--w", type=int, required=True)
    ap.add_argument("--exact", action="store_true", help="weight exactly w (assumes < w excluded)")
    ap.add_argument("--encoding", default="cpsat-onehot", choices=list(SOLVERS))
    ap.add_argument("--time", type=float, default=600)
    ap.add_argument("--no-symmetry", action="store_true")
    ap.add_argument("--orbit-reps", help="comma-separated column orbit representatives")
    ap.add_argument("--cols", help="restrict to these comma-separated columns (model validation)")
    ap.add_argument("--lin", type=int, help="CP-SAT linearization_level (one-hot model only)")
    ap.add_argument("--tag", default="")
    ap.add_argument("--log", help="append the JSON result line to this file")
    a = ap.parse_args()
    if a.matrix:
        G = np.load(a.matrix); p = a.p; cols = None
        if a.subset:
            rng = np.random.default_rng(a.seed); cols = np.sort(rng.choice(G.shape[1], a.subset, replace=False))
            G = G[:, cols]
    else:
        G, p, cols = load_instance(a.instance, a.subset, a.seed)
    if a.cols:
        cols = [int(v) for v in a.cols.split(",")]; G = np.ascontiguousarray(G[:, cols])
        G = G[[t for t in range(G.shape[0]) if G[t].any()]]
    kw = {}
    if a.lin is not None and a.encoding == "cpsat-onehot":
        kw["linearization"] = a.lin
    if a.encoding != "mitm":
        kw["time_limit"] = a.time
        kw["symmetry"] = not a.no_symmetry
        if a.orbit_reps:
            kw["orbit_reps"] = [int(v) for v in a.orbit_reps.split(",")]
    res = SOLVERS[a.encoding](G, p, a.w, exact=a.exact, **kw)
    res.update({"instance": a.instance or a.matrix, "subset": a.subset, "seed": a.seed, "p": p, "tag": a.tag,
                "G_shape": list(G.shape)})
    line = json.dumps(res)
    print(line, flush=True)
    if a.log:
        with open(a.log, "a") as f:
            f.write(line + "\n")


if __name__ == "__main__":
    main()
