# Plan: a Fortran DAQP package (fpm)

A plan for a standalone fpm package: a modern Fortran port of
[DAQP](https://github.com/darnstrom/daqp), a dual active-set solver for
dense convex quadratic programs, with a small object-oriented interface. It
is written to be read on its own: section 10 says what
[sqpopt](https://github.com/jacobwilliams/sqpopt) needs from it, as its
first user, but the package must not depend on sqpopt.

The upstream C code is kept beside the port as a **git submodule**, built
by a helper script (never by fpm's default build), and used to check the
port against the original for both **accuracy** and **speed** (section 5).

Facts about the upstream C code (file and function names, structures,
flags, exit codes, storage order) are as I remember them from DAQP's
repository; check each against the pinned upstream commit before porting
(stage 0 does that, and records the differences here).

---

## 1. What DAQP is, and why port it

- **What it solves.** Dense convex QPs
  \[ \min_x \; \tfrac12 x^T H x + f^T x \quad \text{s.t.} \quad
     b_l \le \begin{bmatrix} x_{1:m_s} \\ A x \end{bmatrix} \le b_u \]
  with \(H\) symmetric positive definite (or semidefinite, through a
  proximal-point outer loop), two-sided bounds on every row, the first
  \(m_s\) rows being simple bounds on variables, and per-row flags for
  equality, "immutable" (never leaves the working set), and soft
  constraints.
- **How.** It turns the QP into a *least-distance problem* (LDP): with the
  Cholesky factor \(H = R^TR\), the variable \(u = Rx + R^{-T}f\) and the
  rows \(M = AR^{-1}\), the QP becomes \(\min \tfrac12\lVert u \rVert^2\)
  subject to \(Mu \le d\). That is solved by a *dual* active-set method,
  which keeps the \(LDL^T\) factors of \(M_{\mathcal{W}} M_{\mathcal{W}}^T\)
  (the working set's rows) and **updates** them when a row is added or
  removed, in \(O(n\,|\mathcal{W}|)\) operations instead of refactoring.
  Simple bounds are handled without forming their rows. The dual method
  starts feasible for the dual (the unconstrained minimizer), so it needs
  no phase 1.
- **Why it suits small and medium dense QPs.** One Cholesky factorization of
  \(H\) per QP (\(n^3/3\)), then cheap updates per iteration, warm starts
  from a given working set, and very little overhead. Its authors report it
  as one of the fastest solvers for small dense QPs (model predictive
  control, where it was developed).
- **Why a Fortran port.** sqpopt's dense QP solver
  (`src/sqpopt_qp_dense_module.f90`, the default for up to 200 variables)
  refactors everything at each active-set iteration: a Householder QR of the
  working set with the full \(n \times n\) \(Q\), the reduced Hessian
  \(Z^THZ\), and its Cholesky factor, \(O(n^3)\) each time (roadmap F30).
  DAQP's updates are the algorithmic fix. A Fortran port gives that with
  nothing but fpm: no C compiler setup, any real kind, copyable types.
- **What it isn't.** A solver for nonconvex QPs (it needs \(H \succ 0\), or
  \(H \succeq 0\) with the proximal loop), for large sparse QPs (\(H\) and
  \(A\) are dense), or for an \(\ell_1\) elastic mode as sqpopt's dense QP
  does it (see section 10).

**Upstream:** <https://github.com/darnstrom/daqp> (Daniel Arnström; MIT
licence, to confirm). The method: D. Arnström, A. Bemporad, D. Axehill,
*A dual active-set solver for embedded quadratic programming using
recursive LDLᵀ updates*, IEEE Transactions on Automatic Control 67(8),
2022. Upstream also has a branch-and-bound layer for binary variables
(BnB-DAQP) and hierarchical (lexicographic) QPs; both are out of scope here
(section 4, stage 6).

---

## 2. Upstream C interface (to port first)

Verified against the pinned commit (v0.10.3,
`29a8987d8b9a596c1ab6acf84c17d38b068b73bc`, stage 0). Upstream has grown
since the first draft of this plan: the differences are marked **(changed)**.

| C entity | What it is |
|---|---|
| `DAQPProblem` | `n`, `m` (rows including the simple bounds), `ms` (simple bounds, the first rows), `H` (\(n \times n\), **row-major**, read as \(\tfrac12(H+H^T)\) by the factorization), `f`, `A` (\((m-m_s) \times n\), **row-major**: confirmed), `bupper`, `blower`, `sense` (an `int` per row), plus `break_points`, `nh` (hierarchical QPs) and `problem_type` (1 = AVI, 2 = prefactorized `H`) **(changed)** |
| `DAQPSettings` | `primal_tol` (1e-6), `dual_tol` (1e-12), `zero_tol` (1e-11), `pivot_tol` (1e-8), `progress_tol` (1e-14), `cycle_tol` (10), `iter_limit` (10000), `fval_bound` (1e30), `eps_prox` (**-1e-6**: negative = automatic, only when `H` is singular), `eta_prox` (**-1**: automatic), `rho_soft` (1e-6), `rel_subopt`, `abs_subopt` (0, BnB), `sing_tol` (3.7e-11), `refactor_tol` (1e-9), `time_limit` (0), `w_soft` (0), `eq_reduction` (AUTO) **(changed)**; `daqp_default_settings(&settings)` fills the defaults |
| `DAQPResult` | `x`, `lam`, `fval`, `soft_slack`, `exitflag`, `iter`, `nodes`, `solve_time`, `setup_time` |
| `daqp_quadprog(result, problem, settings)` | the one-call interface: setup (with the unconstrained-optimum shortcut and, with `eq_reduction` AUTO, the equality elimination), solve, free |
| `setup_daqp(problem, work, &setup_time)`, `daqp_solve(result, work)`, `free_daqp_workspace(work)` + `free_daqp_ldp(work)` | the split interface; the workspace keeps a pointer to the problem |
| `daqp_update_ldp(mask, work, problem)` **(changed name)** | recompute what the mask says: `DAQP_UPDATE_Rinv` 1, `_M` 2, `_v` 4, `_d` 8, `_sense` 16, `_hierarchy` 32, `_unconstrained` 64, `_eliminate` 128 |
| sense flags | `ACTIVE` 1, `LOWER` 2, `IMMUTABLE` 4, `SOFT` 8, `BINARY` 16, `SLACK_FIXED` 32, `SET_ASIDE` 64, `AUTO_EQUALITY` 128; an equality is `ACTIVE+IMMUTABLE` (5), or equal bounds (detected at setup) |
| exit flags | optimal 1, soft-optimal 2, **optimal-inexact 4** (after the noise floor), infeasible -1, cycling -2, unbounded -3, iteration limit -4, nonconvex -5, overdetermined initial working set -6, **time limit -7, unsupported -8** |

Source files: `daqp.c` (the LDP loop `daqp_ldp`), `factorization.c` (the
LDLᵀ updates), `auxiliary.c` (CSP, add/remove, blocking, refinement, soft
constraints, noise floor), `utils.c` (LDP transformation, the factorization of
`H` with deferred inverse and automatic regularization, the unconstrained
check), `api.c` (setup, solve, result), `daqp_prox.c` (proximal and
semi-proximal outer loop, with accelerated steps), and, not ported: `bnb.c`,
`hierarchical.c`, `avi.c`, `eq_elim.c` (1272 lines: elimination of
equalities), `codegen/`. Upstream's own tests are in the interfaces (Julia
`core_tests.jl`, Python, MATLAB, C++/Eigen); there are no plain C tests, and
the build is CMake (not installed here: `tools/build_upstream.sh` compiles the
C files directly).

**Inputs the C code assumes:** `H` dense, row-major, symmetrized by the
factorization (the proximal steps read the raw rows); `A` row-major;
infinite bounds as \(\pm 10^{30}\) (`DAQP_INF`); `sense` zero for a cold
start, or with `ACTIVE` (and `LOWER`) flags for a warm start.

**Multiplier sign:** \(Hx + f + A_{all}^T\lambda = 0\), \(\lambda > 0\) at an
upper bound, \(< 0\) at a lower bound. **Soft penalty:** `w_soft*s +
s^2/(2*rho_soft)` per violated side (normalized rows): quadratic by default.

### Status of the port (2026-10-04)

* Stages 0–4 are done, stage 5 except the release tag: `src/daqp_core.f90`
  (the port), `src/daqp_module.f90` (`daqp_type`), six test programs (all
  passing in REAL32, REAL64, REAL128), two examples, the comparison harness
  (`compare/`, `tools/run_compare.sh`) and `compare/RESULTS.md`.
* With `-ffp-contract=off` in both builds, the port is **bit-for-bit
  identical** to the C code on all 756 comparison problems (same exit flags,
  iterations, working sets, and values); with the optimized builds (-O3), the
  same flags, iterations, and working sets, and values equal to round-off.
* Speed (Apple M5, gfortran 15 `-O3 -funroll-loops`, clang 21 `-O3`): median
  ratio port/C 1.06 (setup), 0.98 (cold solve), 0.97 (warm solve); at most
  1.12 for n ≥ 20.
* Choices: `daqp_inf` = 1e30 as upstream; dependent rows of a warm-start
  working set are dropped (as upstream's activation does); the default
  tolerances are floored at multiples of `epsilon` (only changes single
  precision); `fval` with neither `f` nor a proximal term is `0.5*||u||^2`
  (upstream leaves it unset); the core setup also accepts `A` transposed.
* Not done: set 5 (QPs captured from sqpopt) and set 6 (Maros–Mészáros), which
  need sqpopt's capture option; the `v0.1.0` tag.

---

## 3. Package layout

```
daqp-fortran/
  fpm.toml
  LICENSE            MIT (see section 9), with upstream's copyright line
  NOTICE             credits DAQP's author and the paper
  README.md
  ford.md
  upstream/daqp/     git submodule: the C repository, pinned to a release tag
  src/
    daqp_kinds.F90       real and integer kinds (preprocessor-selected)
    daqp_core.f90        the one-to-one port: LDP setup, the dual active-set loop,
                         the LDLᵀ updates, the proximal loop
    daqp_module.f90      the object-oriented interface (daqp_type) and the public API
  test/
    test_core_*.f90      ports of upstream's tests
    test_kkt_random.f90  random feasible QPs, checked by their KKT conditions
    test_degenerate.f90  dependent rows, ties, cycling guard
    test_warm_start.f90  warm and hot starts give the same solution, in fewer iterations
    test_infeasible.f90  infeasible, unbounded, and semidefinite cases
    test_kinds.f90       run in each real kind by CI
  compare/               NOT part of the package build (section 5)
    fpm.toml             a separate fpm project: depends on the package by path
    c_binding.f90        bind(c) interfaces to the upstream library
    compare.f90          the comparison driver
    qp_sets/             generators and captured QPs
  tools/
    build_upstream.sh    builds upstream/daqp into build/upstream/libdaqp.a
  example/
    example_simple.f90   a small QP, one call
    example_mpc.f90      repeated solves with warm starts
```

Conventions (from the author's other Fortran libraries, and as in
qdldl-fortran):

- `implicit none` everywhere; reals `real(wp)`, literals `1.0_wp`.
- Kind selection: `-DREAL32`, `-DREAL64` (default), `-DREAL128`.
- **No local variable initialized in its declaration** (implicit `save`).
- **Nothing sized by the problem on the stack:** no automatic arrays, no
  array temporaries of size `n` or `m`; the workspace is allocatable
  components of the type, allocated in `setup`, so `solve` allocates
  nothing.
- Allocations use `stat=` and report out of memory through the status.
- FORD comments: `!>` module headers, `!!` for every procedure, type,
  component, and dummy argument.
- **Column-major storage** in the Fortran interface. If upstream's \(A\) is
  row-major, the port stores the rows of \(A\) (and of \(M\)) as columns
  of an \(n \times m\) array, so that a row stays contiguous (section 7).
- 1-based indexing throughout.

---

## 4. Stages

### Stage 0: upstream as a submodule, built separately

- `git submodule add https://github.com/darnstrom/daqp upstream/daqp`, then
  check out a **release tag** and commit the pinned commit. The README says
  `git clone --recursive` (or `git submodule update --init`) is only needed
  for the comparison, never for using the package.
- fpm must not see the C sources: they live outside `src/`, and the
  package's `fpm.toml` doesn't list them. A user of the package (an fpm
  dependency by git URL) never compiles C.
- `tools/build_upstream.sh` builds the C library with upstream's own build
  (CMake, as I recall; else compile the `.c` files directly with `cc -O3`)
  into `build/upstream/libdaqp.a`, in **double precision** and with the same
  optimization level used for the Fortran comparison builds (record both in
  the comparison output).
- Run upstream's own tests once, to know the C library is sound before
  comparing against it.
- Read the pinned source and correct section 2 of this plan (names,
  storage order, flags, exit codes, tolerances, defaults).

**Done when:** `tools/build_upstream.sh` produces the library from a clean
clone, upstream's tests pass, and section 2 matches the pinned commit.

### Stage 1: the one-to-one port (`daqp_core`)

Port the solver as plain procedures and a workspace type, keeping
upstream's algorithm, order of operations, and tolerances, so that the two
can be compared line by line and iteration by iteration:

- the LDP setup: the Cholesky factor of \(H\) (and its inverse, as upstream
  keeps it), \(M = AR^{-1}\), \(v = R^{-T}f\), the scaled bounds; the
  special treatment of simple bounds;
- the dual active-set loop: the step direction, the blocking constraint,
  adding and removing rows, the primal feasibility check at a dual-feasible
  point, the termination tests;
- the \(LDL^T\) updates of the working set's Gram matrix (add a row: one new
  row of \(L\) and a pivot; remove a row: the downdate and the shifted
  rows), with the `pivot_tol` and `zero_tol` checks;
- the cycling guard and the iteration limit;
- the proximal outer loop for semidefinite \(H\) (`eps_prox`).

Not in stage 1: branch and bound, hierarchical QPs, timing.

- **Indexing.** Convert to 1-based; check every loop bound, every `-1`
  marker, every `< 0` test, and the bit flags of `sense` (use `iand`/`ior`
  with named constants, or separate logical arrays: decide once).
- **Storage.** Decide the Fortran layout of \(A\) and \(M\) (rows as
  columns, section 7) and port the indexing with it consistently.
- **Exit codes** as named integer constants with upstream's meaning.

**Done when:** upstream's test problems give the same exit flags, the same
iteration counts, and solutions equal to round-off (checked in stage 4's
harness, which can already run here on a few problems).

### Stage 2: the object-oriented interface (`daqp_type`)

```fortran
type(daqp_type) :: qp
call qp%setup(H, f, A, bupper, blower, istat [, ms] [, sense])  ! once per QP (factor H, form the LDP)
call qp%solve(x, lam, istat [, working_set])                     ! solve, from a cold or warm start
call qp%update(istat [, H] [, f] [, A] [, bupper] [, blower])    ! new data, same sizes: recompute only what changed
call qp%get_working_set(active, at_lower)                         ! for a later warm start
call qp%info(iter=..., fval=..., ...)
call qp%destroy()
```

- **Data in Fortran layout:** `H(n,n)`, `A(m-ms, n)` column-major (the
  interface copies into the internal layout once, in `setup`), bounds of
  size `m`. Accept an absent `A` (bounds only).
- **Warm start:** a working set given as the indices of the active rows and
  which bound each is at (`at_lower`); checked for consistency (no
  duplicates, no more active rows than `n` unless dependent rows are
  handled, section 7).
- **Options as components with defaults** (upstream's settings), and a
  `set_defaults` for clarity.
- **Status codes:** success; soft-optimal; infeasible; unbounded;
  nonconvex (Cholesky failed and no proximal term); cycling; iteration
  limit; invalid input; out of memory; not set up.
- **Multipliers:** document the sign convention exactly (which sign at a
  lower bound, which at an upper bound), and give a helper that converts to
  the common "\(\lambda \ge 0\) at a lower bound" convention if upstream's
  differs.
- **Copyable:** allocatable components only.

**Tests:** the stage-1 problems through the type; `update` with new \(f\)
and bounds only (no refactorization) gives the same as a fresh setup; a
warm start from the optimal working set finishes in 0 or 1 iterations with
the same solution; a copied object solves independently.

### Stage 3: the package's own tests

- **KKT checks** on random feasible QPs of sizes 2 to 500: stationarity
  \(Hx + f + A^T\lambda = 0\) (with the package's sign convention), primal
  feasibility, dual feasibility, complementarity, each to a tolerance
  scaled by \(\epsilon\) and the data's norms.
- **Degenerate problems:** duplicated rows, linearly dependent active
  rows, more active constraints than variables at the solution, ties in the
  ratio test.
- **Infeasible and unbounded** problems, with the right exit codes.
- **Semidefinite \(H\)** with the proximal loop (including \(H = 0\): an LP).
- **Equality and immutable rows; soft constraints** (and what their
  penalty is: verify whether `rho_soft` penalizes the slack linearly or
  quadratically, and document it).
- **Kinds:** every test in single, double, and quadruple precision, with
  tolerances scaled by `epsilon(1.0_wp)`.

### Stage 4: the comparison with upstream (accuracy and speed)

See section 5 for the harness. Run it on every problem set, record the
results in `compare/RESULTS.md` (with the upstream commit, compilers, flags,
and machine), and fix the port until the accuracy criteria hold. Then tune
the port where it is slower than the C code by more than the speed target.

**Done when:** the criteria of section 5 hold on every set, and the results
are recorded.

### Stage 5: performance, kinds, docs, release

- **Speed:** profile the port on the comparison sets; the usual causes of a
  slow port are strided access (a row of a column-major array in an inner
  loop), array temporaries, and `merge`/array syntax in hot loops; fix
  those first.
- **Stack:** no automatic arrays or problem-sized temporaries
  (`-Warray-temporaries`, and a grep for explicit-shape locals).
- **Docs:** README (what it is, when to use it, the convexity requirement,
  the warm start, the sign convention, the licence, how to run the
  comparison), FORD API docs.
- **CI:** gfortran (and ifx if convenient) with the kind matrix; a separate
  job that checks out the submodule, builds upstream, and runs the
  comparison on a small set (so the port can't drift from upstream
  unnoticed when either changes).
- **Release:** tag `v0.1.0`.

### Stage 6 (optional, later)

- Branch and bound for binary variables (BnB-DAQP), hierarchical QPs: only
  if a user needs them.
- Rank-one updates of \(H\)'s factor (\(H\) changes by a BFGS update between
  SQP iterations), if sqpopt's measurements show the per-QP Cholesky
  matters.

---

## 5. Comparing the port with upstream

### The harness

- `compare/` is a **separate fpm project** (its own `fpm.toml`), depending
  on the package by path. It is built only when wanted, with the upstream
  library linked: `fpm run --link-flag "-L../build/upstream -ldaqp"`
  (or a `pixi` / shell task wrapping `tools/build_upstream.sh` and that
  command). The package's own build never needs the C library.
- `c_binding.f90`: `bind(c)` derived types mirroring `DAQPProblem`,
  `DAQPSettings`, and `DAQPResult` (`c_double`, `c_int`, `type(c_ptr)` for
  the arrays), and interfaces to `daqp_default_settings`, `daqp_quadprog`,
  and the split interface. Check every field's order and type against the
  pinned headers: a mismatch here produces plausible-looking garbage.
- The driver builds each QP once, in the Fortran layout, converts it to
  upstream's layout for the C call (transposing \(A\) if upstream is
  row-major), and solves it with both, with **the same settings** (copied
  field by field from the C defaults, so a difference in defaults can't
  pass for a bug).

### Problem sets

1. **Upstream's test problems** (the ones its own tests generate or read).
2. **Random QPs:** \(n\) from 2 to 500, \(m\) from 0 to \(3n\), with \(H\)
   of condition number \(10^0\) to \(10^{10}\), a known feasible point, a
   prescribed number of active constraints at the solution, and fixed
   seeds.
3. **Degenerate QPs:** dependent and duplicated rows, more active rows
   than variables, weakly active constraints (zero multiplier at a bound).
4. **Warm-start sequences:** the MPC-like sequence (the same QP with
   shifted data), solved with each solver's working set from the previous
   solve.
5. **QPs captured from sqpopt:** an option in sqpopt's HS harness writes
   every QP subproblem (dense \(H\), \(g\), the Jacobian, the bounds of the
   step) to files; the driver reads and solves them. These are the
   problems that matter for section 10: small, often degenerate, sometimes
   ill-conditioned.
6. **Maros–Mészáros** convex QPs that are small enough to be dense
   (optional; from the CUTEst/Maros–Mészáros collection already used by
   `tools/cutest_benchmark.py`'s setup, if convenient).

### Accuracy criteria (per problem)

| Quantity | Pass when |
|---|---|
| exit flag | equal |
| iterations | equal (a difference is allowed only where a tie in the ratio test or in the choice of the constraint to drop is broken differently by round-off; each such case is listed and explained) |
| solution | \(\lVert x_F - x_C \rVert_\infty \le c\,\epsilon\,\kappa\,(1 + \lVert x_C \rVert_\infty)\), with \(\kappa\) the condition number of \(H\) (or of the working set's Gram matrix) and \(c\) a small constant |
| multipliers | the same relative test on \(\lambda\), after converting to one sign convention |
| objective | \(|f_F - f_C| \le c\,\epsilon\,(1 + |f_C|)\) |
| final working set | equal (same rows, same bounds) |
| KKT residuals | each solver's own residuals, computed by the driver in quadruple precision from its \(x\) and \(\lambda\); the port's must not be worse than \(c\) times the C code's |

Where iterations differ, compare the KKT residuals of both solutions: the
port is acceptable if both are optimal to tolerance; a difference that
changes the exit flag or the working set at the solution is a bug until
explained.

### Speed

- Time **setup** and **solve** separately, each as the median of many
  repetitions (enough for 0.1 s of total time per problem), with the clock
  around the calls only; report the time per QP and the ratio port / C.
- Sizes \(n\) = 5, 10, 20, 50, 100, 200, 500, with \(m = n\) and \(m = 3n\);
  cold and warm starts.
- Builds: gfortran and gcc at `-O3 -march=native` (and `-O2`), the same
  machine; record compiler versions. Optionally ifx/icx.
- **Target:** the port within **1.5×** of the C code on the median of each
  size, and no worse than 2× on any problem of size 20 or more. (Below
  that, call overhead dominates both, and the ratio means little.)
- Also time sqpopt's current dense QP on set 5, for section 10.

### Output

`compare/RESULTS.md`: one table per set (problems, pass counts, the
largest deviations, the iteration mismatches with their explanations), the
speed table, and the environment. Regenerate it whenever the port or the
pinned upstream changes.

---

## 6. Public API (target)

```fortran
module daqp_module
    integer, parameter, public :: daqp_wp   ! real kind

    ! exit statuses (named constants; 0 or positive = solved)
    integer, parameter, public :: daqp_optimal, daqp_soft_optimal, daqp_infeasible, &
                                  daqp_unbounded, daqp_nonconvex, daqp_cycling, &
                                  daqp_iteration_limit, daqp_invalid_input, &
                                  daqp_out_of_memory, daqp_not_setup

    type, public :: daqp_type
        ! options (upstream's settings, with its defaults)
        real(daqp_wp) :: primal_tol, dual_tol, zero_tol, pivot_tol, progress_tol, cycle_tol
        integer       :: iter_limit
        real(daqp_wp) :: eps_prox, eta_prox, rho_soft
        ! results
        integer       :: iter
        real(daqp_wp) :: fval
        ! private: R (or its inverse), M, v, d, the working set, L and D of its Gram matrix,
        !          the work arrays
    contains
        procedure :: setup, update, solve, get_working_set, info, destroy
    end type
end module
```

The low-level procedures of `daqp_core` stay public too, for callers that
manage their own storage.

---

## 7. Pitfalls to watch

- **Storage order.** If upstream's \(A\) and \(M\) are row-major, the hot
  loops walk a row; in Fortran that row must be contiguous (store
  \(M^T\), \(n \times m\)), or the port will be several times slower.
- **Bit flags** in `sense`: port them with named constants and `iand`, or
  split them into logical arrays; never as magic numbers.
- **Triangles of \(H\):** which one upstream reads, and whether the
  interface should symmetrize.
- **Simple bounds** (the first `ms` rows) take a different code path
  upstream; test with and without them.
- **Infinite bounds:** upstream's threshold for "no bound"; map the
  caller's large values to it.
- **Dependent rows** in a warm-start working set: upstream reports an
  overdetermined initial working set; the Fortran interface should either
  report it the same way or drop the dependent rows (decide, document).
- **Tolerances** are absolute in places upstream; keep upstream's meaning
  in the port (for the comparison), and offer scaled ones only as an
  option.
- **Ties and round-off:** the iteration-by-iteration comparison will show
  differences where two candidates tie; don't "fix" those by changing the
  algorithm, explain them (section 5).
- **Timing inside the solver** (upstream's `solve_time`, `setup_time`):
  keep it out of the core, or behind an option, so the comparison times
  the same work.

---

## 8. Validation matrix (what "done" means)

| Check | Problems | Pass when |
|---|---|---|
| Upstream tests | upstream's cases | same exit flags and solutions |
| Port vs C, accuracy | section 5, sets 1–5 | the criteria of section 5 |
| Port vs C, speed | section 5, sizes 5–500 | within 1.5× (median), 2× (any, \(n \ge 20\)) |
| KKT | random and degenerate | residuals at round-off level |
| Warm start | sequences | same solutions, fewer iterations |
| Edge cases | infeasible, unbounded, semidefinite, LP | right statuses |
| Kinds | REAL32, REAL64, REAL128 | all tests pass |
| Memory | huge `n` | out-of-memory status, no crash |
| Stack | the whole library | no automatic arrays or problem-sized temporaries |

---

## 9. Licensing

- DAQP is, as I recall, MIT-licensed: confirm in the pinned commit. A port
  is a derivative work: keep upstream's copyright notice in the `LICENSE`
  (as the original author, alongside the port's), credit the paper in
  `NOTICE` and the README, and mark the ported files as changed (a line at
  the top of each).
- The submodule is upstream's code under its own licence; it isn't part of
  the released package's sources (fpm doesn't build it), but the repository
  that contains it should say so in the README.
- MIT is compatible with sqpopt's MIT licence. (Not legal advice; check if
  it matters.)

---

## 10. The first user: sqpopt

sqpopt's dense QP solver, `solve_dense_qp` in
`src/sqpopt_qp_dense_module.f90`, is selected by `options%qp_solver_mode`
(`sqpopt_qp_dense`, or `sqpopt_qp_auto` for \(n \le\)
`qp_solver%auto_dense_max_n` = 200, so it solves every QP of the HS
suite). Its QP is

\[ \min_p \; g^Tp + \tfrac12 p^THp \quad \text{s.t.} \quad
   c_l - c \le Jp \le c_u - c, \quad x_l - x \le p \le x_u - x \]

which maps directly onto DAQP's form: the variable bounds are DAQP's simple
bounds (`ms = n`), the linearized constraints its general rows, two-sided
in both. What the dense solver does beyond that, and how DAQP would cover
it:

| sqpopt's dense QP does | With DAQP |
|---|---|
| **\(\ell_1\) elastic mode**: each linearized constraint violated at \(p = 0\) gets a slack with a linear penalty \(\rho\), raised if needed, so an inconsistent QP still has a solution (`sqpopt_infeasible` if slacks remain at the largest weight) | Not native. Either explicit slack variables (\(n\) grows by the number of elastic rows, \(H\) becomes semidefinite in them, so the proximal loop is needed), or DAQP's soft constraints if their penalty can be made linear (check: if it is quadratic, the multipliers and the infeasibility test differ from sqpopt's). Or: try DAQP first, and fall back to the current solver when DAQP reports infeasible. |
| **Forced elastic mode** (`elastic_multiplier_limit`) | the same question |
| **Indefinite \(H\)** (exact Hessian, SR1): follows negative curvature on the face, and reports it (the outer loop then shifts the Hessian) | DAQP's Cholesky of \(H\) fails: report `negative_curvature` (or a QP failure) at once, and let the outer loop shift \(H\), or fall back to the current solver. With the default L-BFGS, \(H\) is positive definite, so this doesn't arise. |
| **Warm start** from the previous QP's working set | DAQP's warm start (`sense` flags), with the working set mapped between the two numbering schemes |
| **Multipliers** with sqpopt's sign (\(\lambda \ge 0\) at a lower bound) | convert from DAQP's convention (section 4, stage 2) |
| **Dependent rows** never enter the working set (an independence check) | DAQP's handling of dependent rows (its overdetermined-working-set exit, and its pivot tolerance): verify on the captured QPs (section 5, set 5) |
| **Out of memory** reported as `sqpopt_out_of_memory` | the package's out-of-memory status |
| **Any real kind** | the package's kinds |

How it would be introduced in sqpopt (a roadmap item, not part of this
package):

1. Add the package as an fpm dependency.
2. A new QP mode, or a backend inside `sqpopt_qp_dense` chosen by an option
   (for example `qp_solver%dense_qp%method`), that forms the dense \(H\)
   (from `hessian_dense`, as now) and the dense Jacobian, calls DAQP, and
   falls back to the current active-set code when DAQP reports infeasible,
   nonconvex, cycling, or an iteration limit. The fallback keeps every
   current behaviour (elastic mode, negative curvature) available.
3. Measure on the HS suite, release build (the protocol of CLAUDE.md): the
   default configuration, `--hessian=exact`, and `--hessian=sr1`, comparing
   solved/local/failed, `fc`, and time against the current dense QP, and
   counting how often the fallback runs. Also the dense QP's time alone,
   from the captured QPs (section 5, set 5), at \(n\) up to 200 and beyond:
   if DAQP is much faster, re-measure `auto_dense_max_n`.
4. Adopt it as the default only if the HS results don't regress; then
   update the guide (QP solver section, "Options that form dense
   matrices"), the Performance page, and the roadmap (F30).

---

## References

- D. Arnström, A. Bemporad, D. Axehill, *A dual active-set solver for
  embedded quadratic programming using recursive LDLᵀ updates*, IEEE
  Transactions on Automatic Control 67(8), 4362–4369, 2022 (check the
  pages).
- D. Arnström, A. Bemporad, D. Axehill, *BnB-DAQP: a mixed-integer QP
  solver for embedded applications*, IFAC World Congress, 2023 (the
  branch-and-bound layer; out of scope).
- D. Goldfarb, A. Idnani, *A numerically stable dual method for solving
  strictly convex quadratic programs*, Mathematical Programming 27, 1–33,
  1983 (the classical dual active-set method DAQP builds on).
- DAQP: <https://github.com/darnstrom/daqp>.
