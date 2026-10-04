# daqp-fortran

A modern Fortran port of [DAQP](https://github.com/darnstrom/daqp), the dual active-set
solver for dense convex quadratic programs by Daniel Arnström, with an object-oriented
interface. It is built with [fpm](https://fpm.fortran-lang.org) and needs no C compiler.

It solves

```
min  0.5 x'Hx + f'x
s.t. blower(1:ms)  <= x(1:ms) <= bupper(1:ms)      (simple bounds)
     blower(ms+1:) <= A x     <= bupper(ms+1:)     (general constraints)
```

with every row two-sided (`daqp_inf` for no bound), and optional per-row flags for
equality, immutable, and soft constraints.

## When to use it

* **Small and medium dense QPs** (up to a few hundred variables): the Hessian is factored
  once per QP, and the working set's factorization is updated, not recomputed, at each
  iteration.
* **Sequences of similar QPs** (model predictive control, SQP): new data can be given
  without refactoring, and each solve starts from the previous working set.
* **Convex problems only.** `H` must be symmetric positive definite. A positive
  semidefinite `H` (including `H = 0`, or no `H`: an LP) is handled by a proximal-point
  outer loop. An indefinite `H` is reported as `daqp_nonconvex`. Large sparse problems
  are better served by a sparse solver.

## Usage

```fortran
use daqp_module

type(daqp_type) :: qp
real(daqp_wp) :: x(n), lam(m)
integer :: istat

call qp%setup(H, f, A, bupper, blower, istat)   ! factor H, form the problem
call qp%solve(x, lam, istat)                    ! istat > 0: solved

call qp%update(istat, f=f2, bupper=bu2)         ! new data, same sizes
call qp%solve(x, lam, istat)                    ! hot start from the previous working set
```

* `H(n,n)`, `f(n)`, `A(m-ms,n)` in the usual Fortran (column-major) layout, and bounds
  `bupper(m)`, `blower(m)`. The number of simple bounds is `ms = m - size(A,1)`: the first
  `ms` entries of the bounds apply to `x(1:ms)`. `H`, `f`, and `A` are optional (use
  keyword arguments after an omitted one): no `H` is an LP, no `A` means only simple bounds.
* Options are components of `daqp_type` (upstream's settings, with its defaults):
  `primal_tol`, `dual_tol`, `zero_tol`, `pivot_tol`, `progress_tol`, `cycle_tol`,
  `iter_limit`, `fval_bound`, `eps_prox`, `eta_prox`, `rho_soft`, `w_soft`, `sing_tol`,
  `refactor_tol`. `set_defaults` restores them.
* Results of the latest solve: `qp%status`, `qp%iter`, `qp%fval`, `qp%soft_slack`,
  `qp%outer_iter` (or `qp%info(...)`).
* Constraint flags: `sense(m)` in `setup`, with `daqp_equality` (or simply equal bounds),
  `daqp_immutable`, `daqp_soft`.
* `daqp_type` is copyable (allocatable components only), and allocates nothing during
  `solve`.

See [`example/example_simple.f90`](example/example_simple.f90) and
[`example/example_mpc.f90`](example/example_mpc.f90).

### Multipliers

`lam` satisfies \( Hx + f + A_{all}^T\lambda = 0 \) with \( A_{all} = [I_{ms}\;0;\,A] \):
`lam(i) > 0` when the upper bound of row `i` is active, `lam(i) < 0` when its lower bound
is active, 0 when inactive. `daqp_lower_positive(lam)` converts to the convention
\( \lambda \ge 0 \) at a lower bound (\( Hx + f = A_{all}^T\lambda \)).

### Warm and hot starts

* By default, `solve` starts from the working set of the previous solve (a *hot start*);
  after `setup`, from the equality constraints only.
* `qp%get_working_set(active, at_lower)` returns the working set of a solution;
  `qp%solve(x, lam, istat, active=active, at_lower=at_lower)` (or `set_working_set`)
  starts from a given one (a *warm start*). Rows that are linearly dependent on the others
  are dropped.
* `qp%solve(..., cold=.true.)` starts from the empty working set.
* `update` with new `f` or bounds keeps the factorization of `H` and the working set; with
  a new `A`, the problem is formed again but `H` is not refactored; with a new `H`,
  everything is recomputed.

### Soft constraints

A row flagged `daqp_soft` may be violated: a violation `s` (of the normalized row) adds
`w_soft*s + s**2/(2*rho_soft)` to the objective (`w_soft = 0` by default: a quadratic
penalty; a large enough `w_soft` gives an exact penalty). The exit status is then
`daqp_soft_optimal` and `qp%soft_slack` is the largest violation.

### Statuses

| status | value | meaning |
|---|---|---|
| `daqp_success` | 0 | `setup`/`update` succeeded |
| `daqp_optimal` | 1 | solved |
| `daqp_soft_optimal` | 2 | solved, a soft constraint is violated |
| `daqp_optimal_inexact` | 4 | solved after persistent cycling; violates a constraint by more than `primal_tol` |
| `daqp_infeasible` | -1 | infeasible |
| `daqp_cycling` | -2 | cycling |
| `daqp_unbounded` | -3 | unbounded (LP) |
| `daqp_iteration_limit` | -4 | iteration limit |
| `daqp_nonconvex` | -5 | `H` is not positive (semi)definite |
| `daqp_overdetermined` | -6 | inconsistent equalities in the initial working set |
| `daqp_invalid_input` | -101 | invalid sizes or indices |
| `daqp_out_of_memory` | -102 | an allocation failed (or the problem is too large) |
| `daqp_not_setup` | -103 | `setup` was not called, or failed |

### Low-level interface

`daqp_core` holds the one-to-one port (a workspace type and plain procedures:
`daqp_setup`, `daqp_solve`, `daqp_update_ldp`, `daqp_quadprog`, `daqp_ldp`, ...), with the
same algorithm, order of operations, and tolerances as the C code. It also accepts `A` in
its internal layout (`At(n,m-ms)`, a row of `A` per column), which saves a transposed copy.

## Real kinds

The real kind is chosen by a preprocessor flag: `-DREAL32`, `-DREAL64` (the default), or
`-DREAL128`, e.g. `fpm test --flag "-DREAL128"`. In double and quadruple precision the
default tolerances are upstream's; in single precision, those below a few `epsilon` are
raised to a multiple of `epsilon`.

## Install

```toml
[dependencies]
daqp-fortran = { git = "https://github.com/jacobwilliams/daqp-fortran" }
```

To build and test: `fpm test`. To build the documentation: `ford ford.md`.

## Comparison with the C code

The upstream C repository is a git submodule in `upstream/daqp`, pinned to the release the
port follows (v0.10.3). It is only needed for the comparison, not to use the package
(`git clone --recursive`, or `git submodule update --init`).

```
tools/run_compare.sh          # writes compare/RESULTS.md
tools/run_compare.sh quick    # a short check (CI)
```

`tools/build_upstream.sh` builds the C library into `build/upstream/libdaqp.a`;
`compare/` is a separate fpm project that solves the same problems with both (upstream's
test problems, random QPs of sizes 2 to 500 and condition numbers up to 1e10, degenerate,
infeasible, semidefinite, soft, LP, and warm-start sequences) and compares exit flags,
iteration counts, working sets, solutions, multipliers, objectives, KKT residuals, and
speed. **With fused multiply-adds disabled in both builds, the port gives bit-for-bit the
same results as the C code on every problem; the optimized builds agree to round-off, and
run at the same speed** (see [`compare/RESULTS.md`](compare/RESULTS.md)).

### Not ported

Branch and bound (binary constraints), hierarchical QPs, affine variational inequalities,
the elimination of equality constraints (an optional reduction for problems with many
equalities; the comparison disables it in the C code), individual weights of soft
constraints, code generation, the time limit, and timing.

## Licence

MIT (see [LICENSE](LICENSE)). DAQP is MIT-licensed, Copyright (c) 2022 Daniel Arnström;
this package is a derivative work (see [NOTICE](NOTICE)).

## Reference

D. Arnström, A. Bemporad, D. Axehill, *A dual active-set solver for embedded quadratic
programming using recursive LDLᵀ updates*, IEEE Transactions on Automatic Control 67(8),
4362–4369, 2022.
